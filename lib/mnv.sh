# shellcheck shell=bash
# Multi-nucleotide variant (MNV) input generation and MNVista calls.
# Sourced by CircuSNV.sh; uses its global settings.

# Normalize and atomize one tool's unfiltered VCF for MNV calling; keeps the plain VCF and
# writes a compressed, indexed copy next to it.
# Usage: normalize_for_mnv <output.vcf> <input.vcf.gz>
normalize_for_mnv() {
  local output_vcf=$1
  shift
  bcftools norm -O v -m -any -f "$Reference" -a --atom-overlaps . -c w "$@" 2>/dev/null | bcftools norm -d all > "$output_vcf" 2>/dev/null
  bgzip -c "$output_vcf" > "$output_vcf.gz"
  bcftools index -f "$output_vcf.gz"
}

# Build mnvcalling/sites.vcf for one sample from the unfiltered per-tool VCFs.
# Usage: generate_mnv_input <samplename>
generate_mnv_input() {
  local samplename=$1
  (
    load_modules

    local mnvpath="$OutputDir/$samplename/mnvcalling"

    echo 'MNV input: normalizing, atomizing, zipping and indexing 4 inputs for MNV analysis --> merging into one list'
    bgzip -c "$mnvpath/"*_M2_F2stats.vcf > "$mnvpath/${samplename}_M2_F2stats.vcf.gz"
    bcftools index -f "$mnvpath/${samplename}_M2_F2stats.vcf.gz"

    # The four tools are independent: normalize, compress and index them in parallel
    local pids=() pid
    normalize_for_mnv "$mnvpath/VarDict_MNV.vcf" "$mnvpath/"*VarDict_inter.vcf.gz & pids+=($!)
    normalize_for_mnv "$mnvpath/Mutect2_MNV.vcf" "$mnvpath/"*M2_F2stats.vcf.gz & pids+=($!)
    normalize_for_mnv "$mnvpath/LoFreq_MNV.vcf" "$mnvpath/"*LoFreq_inter.vcf.gz & pids+=($!)
    normalize_for_mnv "$mnvpath/SiNVICT_MNV.vcf" "$mnvpath/"*SiNVICT_inter_RAW.vcf.gz & pids+=($!)
    for pid in "${pids[@]}"; do
      wait "$pid"
    done

    # isec to merge into one list (variants called by at least 1 tool)
    bcftools isec -p "$mnvpath" -O z -n +1 "$mnvpath/SiNVICT_MNV.vcf.gz" "$mnvpath/Mutect2_MNV.vcf.gz" "$mnvpath/LoFreq_MNV.vcf.gz" "$mnvpath/VarDict_MNV.vcf.gz"
    # Merge sites with data: --> sites2.txt
    pythonscript "$CIRCUSNV_HOME/scripts/merge_variants.py" "$samplename" "phased" --infile "sites.txt" --outfile "sites2.txt" "$OutputDir"

    # MNV filtering on UNMET score
    if [[ "$UNMET" == "true" ]]; then
      unmet_filter "$samplename" "$mnvpath/sites2.txt" "$UNMET_BED" "$UNMET_THRESHOLD" "MNV input -- UNMET filtering"
    fi

    # Convert to vcf
    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$mnvpath/sites2.txt" "$mnvpath/sites.vcf" 'TRUE'
    bgzip -c "$mnvpath/sites.vcf" > "$mnvpath/sites.vcf.gz"
    bcftools index -f "$mnvpath/sites.vcf.gz"
    echo "MNV input: merging 4 tools to 1 list DONE --> sites.vcf"
  )
}

# Run MNVista for one sample; baselines first build their own input, follow-ups reuse the
# MNVs found in their baseline.
# Usage: mnv_calling <samplename> <bam_file> [crosslink_table]
mnv_calling() {
  local samplename=$1
  local bam_file=$2
  local crosslink_table=${3:-}
  local sample_info timepoint corresponding_baseline

  local mnvista_args=("$bam_file" "$OutputDir/$samplename/mnvcalling/sites.vcf" "$OutputDir/$samplename/mnvcalling/results" -T 15 -M 6 -R 150 -F 0 -J 0)
  if [[ -n "$PoN" ]]; then
    mnvista_args+=(-C "$MNV_PON")
  fi

  if [ -n "$crosslink_table" ]; then
    echo "MNV input: Crosslink table provided: $crosslink_table"

    # Extract the timepoint and corresponding baseline for the sample
    sample_info=$(awk -v sample="$samplename" '$1 == sample {print $2, $3}' "$crosslink_table")
    timepoint=$(echo "$sample_info" | awk '{print $1}')
    corresponding_baseline=$(echo "$sample_info" | awk '{print $2}' | tr -d '\r')

    if [ "$timepoint" == "follow-up" ]; then
      echo "MNV input: Sample $samplename is a follow-up sample."
      echo -e "MNV input: corresponding baseline sample = $corresponding_baseline"
      echo -e "MNV input: Copying results.vcf from $corresponding_baseline and analyzing that with MNV algorithm..."

      # The SNVs that make up the phased MNVs in the baseline are the (pre-filtered) input here
      cp "$OutputDir/$corresponding_baseline/mnvcalling/results/results.vcf" "$OutputDir/$samplename/mnvcalling/sites.vcf"
      (
        load_module "$MODULE_HTSLIB"
        load_module "$MODULE_GCC_MNVISTA"
        "$MNVISTA" "${mnvista_args[@]}" -S 1 -N 1 -B 0.00000000001 -P 1
        echo -e "MNV process for sample $samplename DONE ...\n"
      )

    elif [ "$timepoint" == "baseline" ]; then
      echo "MNV input: Sample $samplename is a baseline sample."
      generate_mnv_input "$samplename"
      (
        load_module "$MODULE_HTSLIB"
        load_module "$MODULE_GCC_MNVISTA"
        "$MNVISTA" "${mnvista_args[@]}" -S 10 -N 3 -B 0.95 -P "$PRIOR" -K true
        echo -e "MNV process for baseline sample $samplename DONE ...\n"

        if [[ "${Matched,,}" == "true" ]]; then
          # Matched mode: also run MNVista on the matched control and store it in its own folder
          local matched_control ctrl matched_output
          matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable")
          if [[ -n "$matched_control" ]]; then
            echo "MNV input: Running MNVista on baseline tumor MNVs in matched control ${matched_control}"
            ctrl="${Inputbam%/}/${matched_control}.bam"
            matched_output="$OutputDir/$samplename/mnvcalling/results/matched_control"
            mkdir -p "$matched_output"
            if [[ -f "$ctrl" ]]; then
              "$MNVISTA" "$ctrl" "$OutputDir/$samplename/mnvcalling/results/results.vcf" "$matched_output" -T 12 -M 6 -R 100 -F 0 -J 0 -S 10 -N 3 -B 0.95 -P "$PRIOR" -K true
            else
              echo "WARNING: matched control BAM not found: $ctrl — skipping MNVista"
            fi
          else
            echo "No matched control for sample $samplename — skipping matched MNV analysis"
          fi
        fi
      )

    else
      echo -e "MNV input: Crosslink table is provided but sample $samplename is not annotated as either a 'baseline' or 'follow-up' sample. \nMNV input will be generated for $samplename but not be analyzed"
      generate_mnv_input "$samplename"
    fi

  else
    # Treat as baseline sample
    echo "MNV input: Crosslink table not provided. Treating sample as baseline."
    generate_mnv_input "$samplename"
    (
      load_module "$MODULE_HTSLIB"
      load_module "$MODULE_GCC_MNVISTA"
      "$MNVISTA" "${mnvista_args[@]}" -S 10 -N 3 -B 0.95 -P "$PRIOR" -K true
      echo -e "MNV process for sample $samplename DONE ...\n"
    )
  fi
}
