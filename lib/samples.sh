# shellcheck shell=bash
# Per-sample processing and the loop over all samples.
# Sourced by CircuSNV.sh; uses its global settings.

# Run the whole pipeline for one tumor sample.
# Usage: run_tumor_sample <sample.bam>
run_tumor_sample() {
  local bam_file=$1
  local samplename bam_dir sample_dir
  samplename=$(basename "$bam_file" .bam)
  bam_dir=$(dirname "$bam_file")
  sample_dir="$OutputDir/$samplename"
  local bam_index=("${bam_dir}/${samplename}"*.bai)
  local pid_vardict pid_mutect2 pid_lofreq pid_sinvict
  local status_vardict status_mutect2 status_lofreq status_sinvict
  local MUTS matched_control annotate_input

  echo -e "\nProcessing sample $samplename\n"
  mkdir -p "$sample_dir"
  mkdir -p "$sample_dir/mnvcalling"

  # Create the BAM index if it does not exist
  if [ ! -e "${bam_index[0]}" ]; then
    echo "BAM index file for ${samplename} not found. Creating BAM index..."
    (
      module load "$MODULE_SAMTOOLS"
      samtools index "$bam_file"
    )
  else
    echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
    touch "${bam_index[0]}"
  fi

  echo "Variant calling in tumor sample..."

  # Launch all callers in the background; each retries internally
  run_caller_with_retries Vardict "$bam_file" "$sample_dir/" 3 10 & pid_vardict=$!
  run_caller_with_retries Mutect2 "$bam_file" "$sample_dir/" 3 10 & pid_mutect2=$!
  run_caller_with_retries Lofreq "$bam_file" "$sample_dir/" 3 10 & pid_lofreq=$!
  run_caller_with_retries Sinvict "$bam_file" "$sample_dir/" 3 10 & pid_sinvict=$!

  # Wait for all to finish and capture exit statuses
  wait "$pid_vardict"; status_vardict=$?
  wait "$pid_mutect2"; status_mutect2=$?
  wait "$pid_lofreq"; status_lofreq=$?
  wait "$pid_sinvict"; status_sinvict=$?

  echo "Vardict: $status_vardict, Mutect2: $status_mutect2, LoFreq: $status_lofreq, SiNVICT: $status_sinvict"

  # Merge outputs from all tools, zip and index
  (
    load_modules

    # Variants called by at least $Calls tools (file order SiNVICT, Mutect2, LoFreq, VarDict
    # gives isec/0000-0003 and the SV-M2-LF-VD tool tag)
    bcftools isec -p "$sample_dir/isec" -O z -n "+$Calls" "$sample_dir/"*_SiNVICT.vcf.gz "$sample_dir/"*_Mutect2.vcf.gz "$sample_dir/"*_LoFreq.vcf.gz "$sample_dir/"*_VarDict.vcf.gz
    # Convert to vcf, zip and index
    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$sample_dir/isec/sites.txt" "$sample_dir/isec/sites.vcf" 'FALSE'
    bgzip -c "$sample_dir/isec/sites.vcf" > "$sample_dir/isec/sites.vcf.gz" && bcftools index "$sample_dir/isec/sites.vcf.gz"
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
          module load "$MODULE_SAMTOOLS"
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
  mkdir "$sample_dir/raw-data/" && mkdir "$sample_dir/filtering-info/" && mkdir "$sample_dir/intermediate-files/" && mkdir "$sample_dir/final-list/"
  mv "$sample_dir/filtering_info.txt" "$sample_dir/filtering-info/" && mv "$sample_dir/removed_variants.xlsx" "$sample_dir/filtering-info/"
  mv "$sample_dir/"merged* "$sample_dir/intermediate-files/" && mv "$sample_dir/"Final_list* "$sample_dir/final-list/"

  echo -e "DONE analyzing $samplename!!\n\n\n"
}

# Process every sample: a single BAM, all BAMs in a directory, or (with a crosslink table)
# the baselines first and then the follow-ups.
run_samples() {
  local s bam
  local baselines=() follows=() samples=()

  mkdir -p "$OutputDir"

  if [[ -d "$Inputbam" ]]; then
    if [[ -n "${CrosslinkTable:-}" && -f "$CrosslinkTable" ]]; then
      echo "Using CrosslinkTable: processing baselines --> follow-ups (controls ignored)"

      # Baselines (col1 where col2=="baseline"), deduplicated, first occurrence kept
      mapfile -t baselines < <(
        awk -F $'\t' 'NR>1 && $2=="baseline" && !seen[$1]++{print $1}' "$CrosslinkTable"
      )

      # Follow-ups (col1 where col2=="follow-up"), deduplicated, first occurrence kept
      mapfile -t follows < <(
        awk -F $'\t' 'NR>1 && $2=="follow-up" && !seen[$1]++{print $1}' "$CrosslinkTable"
      )

      samples=( "${baselines[@]}" "${follows[@]}" )

      if [[ ${#samples[@]} -eq 0 ]]; then
        echo "Warning: no baseline or follow-up samples found in '$CrosslinkTable'"
      else
        echo -e "\nSamples order ="
        printf '%s\n' "${samples[@]}"
        for s in "${samples[@]}"; do
          run_tumor_sample "$Inputbam/${s}.bam"
        done
      fi

    else
      # No table: process all BAMs in the directory
      echo "No CrosslinkTable: processing all BAMs in '$Inputbam'"
      for bam in "$Inputbam"/*.bam; do
        run_tumor_sample "$bam"
      done
    fi

  elif [[ -f "$Inputbam" ]]; then
    # Input is a single BAM file
    run_tumor_sample "$Inputbam"
  else
    echo "Invalid Inputbam: $Inputbam" >&2
    usage
  fi
}
