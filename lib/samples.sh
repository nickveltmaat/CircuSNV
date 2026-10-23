# shellcheck shell=bash
# Per-sample processing and the loop over all samples.
# Sourced by CircuSNV.sh; uses its global settings.

# Run the whole pipeline for one tumor sample.
# Usage: run_tumor_sample <sample.bam>
run_tumor_sample() {
  local bam_file=$1
  local samplename sample_dir
  samplename=$(basename "$bam_file" .bam)
  sample_dir="$OutputDir/$samplename"
  local pid_vardict pid_mutect2 pid_lofreq pid_sinvict
  local status_vardict status_mutect2 status_lofreq status_sinvict
  local MUTS matched_control annotate_input

  echo -e "\nProcessing sample $samplename\n"
  mkdir -p "$sample_dir"
  mkdir -p "$sample_dir/mnvcalling"

  ensure_bam_index "$bam_file" "$MODULE_SAMTOOLS"

  echo "Variant calling in tumor sample..."

  # Launch all callers in the background; each retries internally
  run_caller_with_retries Vardict "$bam_file" "$sample_dir/" 3 10 & pid_vardict=$!
  run_caller_with_retries Mutect2 "$bam_file" "$sample_dir/" 3 10 & pid_mutect2=$!
  run_caller_with_retries Lofreq "$bam_file" "$sample_dir/" 3 10 & pid_lofreq=$!
  run_caller_with_retries Sinvict "$bam_file" "$sample_dir/" 3 10 & pid_sinvict=$!

  # Wait for all to finish and capture exit statuses
  # ('wait' returns the job's exit status; capture it instead of letting set -e stop here)
  status_vardict=0; wait "$pid_vardict" || status_vardict=$?
  status_mutect2=0; wait "$pid_mutect2" || status_mutect2=$?
  status_lofreq=0; wait "$pid_lofreq" || status_lofreq=$?
  status_sinvict=0; wait "$pid_sinvict" || status_sinvict=$?

  echo "Vardict: $status_vardict, Mutect2: $status_mutect2, LoFreq: $status_lofreq, SiNVICT: $status_sinvict"
  if (( status_vardict || status_mutect2 || status_lofreq || status_sinvict )); then
    echo "ERROR: variant calling failed for $samplename after all retries (exit codes above); stopping this sample" >&2
    return 1
  fi

  # Merge outputs from all tools, zip and index
  (
    load_modules

    # Variants called by at least $Calls tools (file order SiNVICT, Mutect2, LoFreq, VarDict
    # gives isec/0000-0003 and the SV-M2-LF-VD tool tag)
    bcftools isec -p "$sample_dir/isec" -O z -n "+$Calls" "$sample_dir/"*_SiNVICT.vcf.gz "$sample_dir/"*_Mutect2.vcf.gz "$sample_dir/"*_LoFreq.vcf.gz "$sample_dir/"*_VarDict.vcf.gz
    # Convert to vcf, zip and index
    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$sample_dir/isec/sites.txt" "$sample_dir/isec/sites.vcf" 'FALSE'
    bgzip -c "$sample_dir/isec/sites.vcf" > "$sample_dir/isec/sites.vcf.gz"
    bcftools index -f "$sample_dir/isec/sites.vcf.gz"
    # PoN filtering (or not)
    if [ -z "$PoN" ]; then
      echo -e "No PoN argument given, proceeding with intersect results..."
      echo "Adding mutation metrics back to merged mutation list"
      pythonscript "$CIRCUSNV_HOME/scripts/merge_variants.py" "$samplename" "no_pon" "$OutputDir"
    else
      # Filter the combined list with the combined blacklist
      echo -e "PoN argument given, filtering intersect results with PoN..."
      bcftools isec -p "$sample_dir/isec/PoN-filtered" -O v -n~10 "$sample_dir/isec/sites.vcf.gz" "$PON_DIR/BLACKLIST_combined.vcf.gz"

      # Add VAF/RD/MRD data to the merged list
      echo "Adding mutation metrics back to pon-filtered merged mutation list"
      pythonscript "$CIRCUSNV_HOME/scripts/merge_variants.py" "$samplename" "pon" "$OutputDir"
    fi
  )

  # Filtering on UNMET score
  if [[ "$UNMET" == "true" ]]; then
    unmet_filter "$samplename" "$sample_dir/merged_vcfs.txt" "$UNMET_BED" "$UNMET_THRESHOLD" "UNMET filtering"
  fi

  # MNV calling
  if [ "$MNV" == "true" ]; then
    echo -e "\nMNV process start: -H argument is set to 'true', proceeding with creating input list for MNV calling"
    mnv_calling "$samplename" "$bam_file" "$CrosslinkTable"
  elif [ "$MNV" == "false" ]; then
    echo "-H argument is set to 'false'; not proceeding with MNV calling"
    rm -rf "$sample_dir/mnvcalling/"
  else
    echo "-H argument not supplied or is not 'true' nor 'false'; not proceeding with MNV calling"
    rm -rf "$sample_dir/mnvcalling/"
  fi

  # First annotation step: gathers the data post_filtering.py filters on
  echo "First step of annotating mutations..."
  awk -v samplename="$samplename" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" samplename "\t" $5}' "$sample_dir/merged_vcfs.txt" > "$sample_dir/merged_vcfs2.txt"
  annotate "$sample_dir/merged_vcfs2.txt" "merged_vcfs_annotated" excel

  # Filtering annotated output
  echo "Post filtering based on functional effect and public SNP databases..."
  pythonscript "$CIRCUSNV_HOME/scripts/post_filtering.py" "$samplename" "$OutputDir" "$UNMET"

  # Second annotation step, if there are mutations left
  MUTS=$(awk '{ FS="\t" } END { print $2 }' "$sample_dir/filtering_info.txt")
  echo "$MUTS mutations left after filtering..."
  if [[ $MUTS -eq 0 ]]; then
    echo 'NO mutations left after filtering...'
    echo "After post-filtering, no mutations were left to annotate or plot, therefore, these steps are skipped for this sample" > "$sample_dir/NO_MUTS.txt"
  else
    echo 'At least one mutation remains, on to annotating final list...'

    if [[ "${Matched,,}" == "true" ]]; then
      # Find the matched control (strip CR)
      matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable")

      if [[ -n "$matched_control" ]]; then
        (
          load_module "$MODULE_SAMTOOLS"
          add_matched_control_metrics "$sample_dir/merged_vcfs_annotated-filtered.txt" "$samplename"
        )
        annotate_input="$sample_dir/merged_vcfs_annotated-filtered.with_control.txt"
      else
        echo "[MATCHED CTRL] No matched control for $samplename — skipping matched metrics"
        annotate_input="$sample_dir/merged_vcfs_annotated-filtered.txt"
      fi

      annotate "$annotate_input" "Final_list_$samplename" excel vcf
    else
      annotate "$sample_dir/merged_vcfs_annotated-filtered.txt" "Final_list_$samplename" excel vcf
    fi
  fi

  # Sort the results into subdirectories
  mkdir -p "$sample_dir/raw-data/" "$sample_dir/filtering-info/" "$sample_dir/intermediate-files/" "$sample_dir/final-list/"
  mv "$sample_dir/filtering_info.txt" "$sample_dir/filtering-info/"
  mv "$sample_dir/removed_variants.xlsx" "$sample_dir/filtering-info/"
  mv "$sample_dir/"merged* "$sample_dir/intermediate-files/"
  # No final list exists when post-filtering left no mutations
  if compgen -G "$sample_dir/Final_list*" >/dev/null; then
    mv "$sample_dir/"Final_list* "$sample_dir/final-list/"
  fi

  if [[ "$CLEANUP" == "true" ]]; then
    cleanup_sample_intermediates "$sample_dir"
  fi
  mark_sample_done "$samplename"

  echo -e "DONE analyzing $samplename!!\n\n\n"
}

# --cleanup: remove intermediate files that later steps and samples no longer need. Kept:
# the per-tool VCFs, intermediate-files/, filtering-info/, final-list/ and mnvcalling/results/
# (follow-up samples read their baseline's results.vcf).
cleanup_sample_intermediates() {
  local sample_dir=$1
  echo "Removing intermediate files (--cleanup): $sample_dir/isec, $sample_dir/mnvcalling (except results/)"
  rm -rf "$sample_dir/isec"
  if [[ -d "$sample_dir/mnvcalling" ]]; then
    find "$sample_dir/mnvcalling" -mindepth 1 -maxdepth 1 ! -name results -exec rm -rf {} +
  fi
}

# A sample counts as finished once this marker exists; --resume skips finished samples.
sample_done_marker() {
  echo "$(run_state_dir)/done/$1"
}

mark_sample_done() {
  mkdir -p "$(run_state_dir)/done"
  date '+%Y-%m-%d %H:%M:%S' > "$(sample_done_marker "$1")"
}

# Print the state of a sample's output: 'new', 'done' or 'incomplete'.
sample_state() {
  local name=$1
  if [[ -f "$(sample_done_marker "$name")" ]]; then
    echo "done"
  elif [[ -d "$OutputDir/$name" ]]; then
    echo "incomplete"
  else
    echo "new"
  fi
}

# Without --resume, refuse to write into sample directories left by an earlier run: mixing
# old and new files would give misleading results.
check_existing_output() {
  local bam name state
  local found=()
  [[ "$RESUME" == "true" ]] && return 0
  while IFS= read -r bam; do
    name=$(basename "$bam" .bam)
    state=$(sample_state "$name")
    [[ "$state" == "new" ]] || found+=("$name ($state)")
  done < <(list_sample_bams)
  if (( ${#found[@]} )); then
    echo "ERROR: the output directory already has results for: ${found[*]}" >&2
    echo "Use --resume to skip finished samples and redo unfinished ones, or choose another -O/--output." >&2
    exit 1
  fi
}

# --resume: move the output of an unfinished sample out of the way before it is run again.
set_aside_incomplete_sample() {
  local name=$1
  local target
  target="$(run_state_dir)/incomplete/${name}_$(date '+%Y%m%d_%H%M%S')"
  mkdir -p "$(run_state_dir)/incomplete"
  mv "$OutputDir/$name" "$target"
  echo "Unfinished output of $name from an earlier run moved to $target"
}

# Print the sample names from the crosslink table in processing order: the baselines first,
# then the follow-ups (each deduplicated, first occurrence kept; controls are not processed).
crosslink_sample_order() {
  awk -F $'\t' 'NR>1 && $2=="baseline" && !seen[$1]++{print $1}' "$CrosslinkTable"
  awk -F $'\t' 'NR>1 && $2=="follow-up" && !seen[$1]++{print $1}' "$CrosslinkTable"
}

# Print the BAM files to process, one per line, in processing order: a single BAM, all BAMs
# in a directory, or (with a crosslink table) the baselines first and then the follow-ups.
list_sample_bams() {
  local s bam
  if [[ -d "$Inputbam" ]]; then
    if [[ -n "${CrosslinkTable:-}" && -f "$CrosslinkTable" ]]; then
      while IFS= read -r s; do
        echo "$Inputbam/${s}.bam"
      done < <(crosslink_sample_order)
    else
      for bam in "$Inputbam"/*.bam; do
        echo "$bam"
      done
    fi
  elif [[ -f "$Inputbam" ]]; then
    echo "$Inputbam"
  fi
}

# Process every sample in order.
run_samples() {
  local bam
  local bams=()

  mkdir -p "$OutputDir"

  if [[ -d "$Inputbam" ]]; then
    if [[ -n "${CrosslinkTable:-}" && -f "$CrosslinkTable" ]]; then
      echo "Using CrosslinkTable: processing baselines --> follow-ups (controls ignored)"
    else
      echo "No CrosslinkTable: processing all BAMs in '$Inputbam'"
    fi
  elif [[ ! -f "$Inputbam" ]]; then
    echo "Invalid Inputbam: $Inputbam" >&2
    usage
  fi

  mapfile -t bams < <(list_sample_bams)

  if [[ ${#bams[@]} -eq 0 ]]; then
    echo "Warning: no samples to process"
    return 0
  fi
  if [[ -n "${CrosslinkTable:-}" && -f "$CrosslinkTable" && -d "$Inputbam" ]]; then
    echo -e "\nSamples order ="
    crosslink_sample_order
  fi

  local name
  for bam in "${bams[@]}"; do
    if [[ "$RESUME" == "true" ]]; then
      name=$(basename "$bam" .bam)
      case "$(sample_state "$name")" in
        done)
          echo "Skipping $name: finished in an earlier run (--resume)"
          continue
          ;;
        incomplete)
          set_aside_incomplete_sample "$name"
          ;;
      esac
    fi
    run_tumor_sample "$bam"
  done
}

# Describe what a run would do, without writing anything (--dry-run).
print_run_plan() {
  local bams=() bam n=0 steps

  mapfile -t bams < <(list_sample_bams)

  echo -e "\nDry run: nothing is written and no tool is started.\n"
  echo "Output directory: $OutputDir"
  if [[ -z "$PoN" ]]; then
    echo "Panel of normals: none"
  elif [[ -d "$PoN" ]]; then
    echo "Panel of normals: generate from $(compgen -G "$PoN/*.bam" | wc -l) normal BAM(s) in $PoN, work directory $PON_DIR"
  else
    echo "Panel of normals: pre-made blacklist $PoN, work directory $PON_DIR"
  fi

  echo -e "\nSamples (${#bams[@]}), in processing order:"
  local name state note
  for bam in ${bams[@]+"${bams[@]}"}; do
    n=$((n + 1))
    name=$(basename "$bam" .bam)
    state=$(sample_state "$name")
    note=""
    if [[ "$RESUME" == "true" && "$state" == "done" ]]; then
      note="  [finished earlier: skipped]"
    elif [[ "$RESUME" == "true" && "$state" == "incomplete" ]]; then
      note="  [unfinished earlier: set aside and run again]"
    elif [[ "$state" != "new" ]]; then
      note="  [output exists: the real run stops unless --resume is given]"
    fi
    echo "  $n. $name  ($bam)$note"
  done
  [[ "$CLEANUP" == "true" ]] && echo -e "\nIntermediate files are removed after each sample (--cleanup)."

  steps="VarDict + Mutect2 + LoFreq + SiNVICT in parallel ($THREADS threads each)"
  steps+=" -> merge (called by >= $Calls tool(s))"
  [[ -n "$PoN" ]] && steps+=" -> PoN filter"
  [[ "$UNMET" == "true" ]] && steps+=" -> UNMET filter (< $UNMET_THRESHOLD)"
  [[ "$MNV" == "true" ]] && steps+=" -> MNV calling (MNVista)"
  steps+=" -> annotation (OpenCRAVAT) -> post-filtering"
  [[ "${Matched,,}" == "true" ]] && steps+=" -> matched-control read counts"
  steps+=" -> final annotation"
  echo -e "\nSteps per sample:\n  $steps"
}
