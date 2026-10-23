# shellcheck shell=bash
# Multi-nucleotide variant (MNV) input generation and MNVista calls.
# Sourced by CircuSNV.sh; uses its global settings.

generate_mnv_input() {
  (
    load_modules
    
    local mnvpath="$OutputDir/$samplename/mnvcalling"
    
    echo 'MNV input: normalizing, atomizing, zipping and indexing 4 inputs for MNV analysis --> merging into one list'
    bgzip -c $mnvpath/*_M2_F2stats.vcf > $mnvpath/${samplename}_M2_F2stats.vcf.gz && bcftools index $mnvpath/${samplename}_M2_F2stats.vcf.gz 
    
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*VarDict_inter.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/VarDict_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*M2_F2stats.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/Mutect2_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*LoFreq_inter.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/LoFreq_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*SiNVICT_inter_RAW.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/SiNVICT_MNV.vcf 2>/dev/null
    
    bgzip -c $mnvpath/VarDict_MNV.vcf > $mnvpath/VarDict_MNV.vcf.gz && bcftools index $mnvpath/VarDict_MNV.vcf.gz
    bgzip -c $mnvpath/Mutect2_MNV.vcf > $mnvpath/Mutect2_MNV.vcf.gz && bcftools index $mnvpath/Mutect2_MNV.vcf.gz
    bgzip -c $mnvpath/LoFreq_MNV.vcf > $mnvpath/LoFreq_MNV.vcf.gz && bcftools index $mnvpath/LoFreq_MNV.vcf.gz
    bgzip -c $mnvpath/SiNVICT_MNV.vcf > $mnvpath/SiNVICT_MNV.vcf.gz && bcftools index $mnvpath/SiNVICT_MNV.vcf.gz
    
    
    #isec to merge into one list --> now with 1 tools. 
    bcftools isec -p $mnvpath -O z -n +1 $mnvpath/SiNVICT_MNV.vcf.gz $mnvpath/Mutect2_MNV.vcf.gz $mnvpath/LoFreq_MNV.vcf.gz $mnvpath/VarDict_MNV.vcf.gz
    #Merge sites with data: --> sites2.txt
    pythonscript "$CIRCUSNV_HOME"/scripts/merge_variants.py $samplename "phased" --infile "sites.txt" --outfile "sites2.txt" $OutputDir
    
    #### MNV Filtering on UNMET score ###
    if [[ "$UNMET" == "true" ]]; then
      unmet_filter "$samplename" "$mnvpath/sites2.txt" "$UNMET_BED" "$UNMET_THRESHOLD" "MNV input -- UNMET filtering"
    fi  
    
    #convert to vcf
    pythonscript "$CIRCUSNV_HOME"/scripts/isec_to_vcf.py $mnvpath/sites2.txt $mnvpath/sites.vcf 'TRUE'
    bgzip -c $mnvpath/sites.vcf > $mnvpath/sites.vcf.gz && bcftools index $mnvpath/sites.vcf.gz
    echo "MNV input: merging 4 tools to 1 list DONE --> sites.vcf" 
  )
}


mnv_calling() {
  local samplename=$1
  local crosslink_table=$2
  
  mnvista_cmd="$MNVISTA $bam_file $OutputDir/$samplename/mnvcalling/sites.vcf $OutputDir/$samplename/mnvcalling/results -T 15 -M 6 -R 150 -F 0 -J 0"
  [ -n "$PoN" ] && mnvista_cmd+=" -C '$MNV_PON'"
  
  if [ -n "$crosslink_table" ]; then
    echo "MNV input: Crosslink table provided: $crosslink_table"
    
    # Extract the timepoint and corresponding baseline for the sample
    sample_info=$(awk -v sample="$samplename" '$1 == sample {print $2, $3}' "$crosslink_table")
    timepoint=$(echo "$sample_info" | awk '{print $1}')
    corresponding_baseline=$(echo "$sample_info" | awk '{print $2}' | tr -d '\r')
    
    if [ "$timepoint" == "follow-up" ]; then
      # Sample is follow-up
      echo "MNV input: Sample $samplename is a follow-up sample."
      echo -e "MNV input: corresponding baseline sample = $corresponding_baseline"
      echo -e "MNV input: Copying results.vcf from $corresponding_baseline and analyzing that with MNV algorithm..."      

      #fetch list of SNVs (new vcf) that make up the phased MNVs in the baseline sample (so it is pre-filtered already and dont need to look at all variants)
      cp $OutputDir/$corresponding_baseline/mnvcalling/results/results.vcf "$OutputDir/$samplename/mnvcalling/sites.vcf"
      (
        module load $MODULE_HTSLIB
        module load $MODULE_GCC_MNVISTA
        mnvista_cmd+=" -S 1 -N 1 -B 0.00000000001 -P 1"
        eval "$mnvista_cmd"
        echo -e "MNV process for sample $samplename DONE ...\n"
      )
    
    elif [ "$timepoint" == "baseline" ]; then
      # Sample is baseline
      echo "MNV input: Sample $samplename is a baseline sample."
      generate_mnv_input
      (
        module load $MODULE_HTSLIB
        module load $MODULE_GCC_MNVISTA
        mnvista_cmd+=" -S 10 -N 3 -B 0.95 -P $PRIOR -K true"
        eval "$mnvista_cmd"
        echo -e "MNV process for baseline sample $samplename DONE ...\n"
        
        if [[ "${Matched,,}" == "true" ]]; then
          # If matched mode is on, additionally run MNVista on matched control sample and store in new folder
          matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable")
          if [[ -n "$matched_control" ]]; then
            echo "MNV input: Running MNVista on baseline tumor MNVs in matched control ${matched_control}"
            ctrl="${Inputbam%/}/${matched_control}.bam"
            matched_output="$OutputDir/$samplename/mnvcalling/results/matched_control"
            mkdir -p "$matched_output"
            if [[ -f "$ctrl" ]]; then
              "$MNVISTA" "$ctrl" "$OutputDir/$samplename/mnvcalling/results/results.vcf" "$matched_output" -T 12 -M 6 -R 100 -F 0 -J 0 -S 10 -N 3 -B 0.95 -P $PRIOR -K true
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
      generate_mnv_input
    fi
  
  else
    #Treat as baseline sample
    echo "MNV input: Crosslink table not provided. Treating sample as baseline."
    generate_mnv_input
    (
      module load $MODULE_HTSLIB
      module load $MODULE_GCC_MNVISTA
      mnvista_cmd+=" -S 10 -N 3 -B 0.95 -P $PRIOR -K true"
      eval "$mnvista_cmd"
      echo -e "MNV process for sample $samplename DONE ...\n"
    )
  fi
}
