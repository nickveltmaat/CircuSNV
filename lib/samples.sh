# shellcheck shell=bash
# Per-sample processing and the loop over all samples.
# Sourced by CircuSNV.sh; uses its global settings.

###################################################
# Function for running a tumor-sample
run_tumor_sample() {
  local bam_file=$1
  local samplename=$(basename "$bam_file" .bam)
  local bam_dir=$(dirname "$bam_file")
  local bam_index=(${bam_dir}/${samplename}*.bai)

  echo -e "\nProcessing sample $samplename\n" 
  mkdir -p $OutputDir/$samplename
  mkdir -p $OutputDir/$samplename/mnvcalling
  
  # Check if BAM index file exists, if not create it
  if [ ! -e "${bam_index[0]}" ]; then
    echo "BAM index file for ${samplename} not found. Creating BAM index..."
    (
      module load SAMtools/1.19.2-GCCcore-11.3.0
      samtools index "$bam_file"
    )
  else
    echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
    touch $bam_index
  fi
    
  echo "Variant calling in tumor sample..."
  
  
  #NEW
  # launch all wrappers in background; they retry internally
  run_caller_with_retries Vardict "$bam_file" "$OutputDir/$samplename/" 3 10 & pid_vardict=$!
  run_caller_with_retries Mutect2 "$bam_file" "$OutputDir/$samplename/" 3 10 & pid_mutect2=$!
  run_caller_with_retries Lofreq "$bam_file" "$OutputDir/$samplename/" 3 10 & pid_lofreq=$!
  run_caller_with_retries Sinvict "$bam_file" "$OutputDir/$samplename/" 3 10 & pid_sinvict=$!
  
  # wait for all to finish and capture exit statuses
  wait "$pid_vardict"; status_vardict=$?
  wait "$pid_mutect2"; status_mutect2=$?
  wait "$pid_lofreq"; status_lofreq=$?
  wait "$pid_sinvict"; status_sinvict=$?
  
  echo "Vardict: $status_vardict, Mutect2: $status_mutect2, LoFreq: $status_lofreq, SiNVICT: $status_sinvict"
  # handle failures if any non-zero
  
  ### Merge outputs from all tools, zip and index:
  (
    load_modules
    
    #Combine 4 tools' output into one textfile using intersect
    bcftools isec -p $OutputDir/$samplename/isec -O z -n +$Calls $OutputDir/$samplename/*_SiNVICT.vcf.gz $OutputDir/$samplename/*_Mutect2.vcf.gz $OutputDir/$samplename/*_LoFreq.vcf.gz $OutputDir/$samplename/*_VarDict.vcf.gz
    #Convert to vcf, zip and index
    pythonscript "$CIRCUSNV_HOME"/scripts/isec_to_vcf.py $OutputDir/$samplename/isec/sites.txt $OutputDir/$samplename/isec/sites.vcf 'FALSE'
    bgzip -c $OutputDir/$samplename/isec/sites.vcf > $OutputDir/$samplename/isec/sites.vcf.gz && bcftools index $OutputDir/$samplename/isec/sites.vcf.gz
    #PoN filtering (or not)
    if [ -z "$PoN" ]; then
      echo -e "No PoN argument given, proceeding with intersect results..."
      echo "Adding mutation metrics back to merged mutation list"
      pythonscript "$CIRCUSNV_HOME"/scripts/merge_variants.py $samplename "no_pon" $OutputDir
    else
      #Filter combined list with combined blacklist:
      echo -e "PoN argument given, filtering intersect results with PoN..."
      bcftools isec -p $OutputDir/$samplename/isec/PoN-filtered -O v -n~10 $OutputDir/$samplename/isec/sites.vcf.gz $PON_DIR/BLACKLIST_combined.vcf.gz
      
      # Check if variants are left after PoN Filtering ##########################################      
      
      #Add VAF/RD/MRD data to the merged list:
      echo "Adding mutation metrics back to pon-filtered merged mutation list"
      pythonscript "$CIRCUSNV_HOME"/scripts/merge_variants.py $samplename "pon" $OutputDir
    fi
  )  
  
  #### Filtering on UNMET score ###
  if [[ "$UNMET" == "true" ]]; then
    unmet_filter "$samplename" "$OutputDir/$samplename/merged_vcfs.txt" "$UNMET_BED" 0.90 "UNMET filtering"
  fi

  #MNVCalling
  if [ "$MNV" == "true" ]; then
    #If phased argument is true: 
    echo -e "\nMNV process start: -H argument is set to 'true', proceeding with creating input list for MNV calling"
    mnv_calling $samplename $CrosslinkTable
  elif [ "$MNV" == "false" ]; then
    echo "-H argument is set to 'false'; not proceeding with MNV calling"
    rm -rf $OutputDir/$samplename/mnvcalling/
  else
    echo "-H argument not supplied or is not 'true' nor 'false'; not proceeding with MNV calling"
    rm -rf $OutputDir/$samplename/mnvcalling/
  fi
  
  #1st annotation step
  echo "First step of annotating mutations..."
  awk -v samplename="$samplename" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" samplename "\t" $5}' "$OutputDir/$samplename/merged_vcfs.txt" > $OutputDir/$samplename/merged_vcfs2.txt
  annotate $OutputDir/$samplename/merged_vcfs2.txt "merged_vcfs_annotated" excel
    
  #Filtering annotated output
  echo "Post filtering based on functional effect and public SNP databases..." 
  pythonscript "$CIRCUSNV_HOME"/scripts/post_filtering.py $samplename $OutputDir $UNMET
    
  #Annotating again if there are mutations left:
  MUTS=$(awk '{ FS="\t" } END { print $2 }' $OutputDir/$samplename/filtering_info.txt)
  echo "$MUTS mutations left after filtering..."
  if [[ $MUTS -eq 0 ]]; then #If 0 mutations left:
    echo 'NO mutations left after filtering...' 
    echo "After post-filtering, no mutations were left to annotate or plot, therefore, these steps are skipped for this sample" > $OutputDir/$samplename/NO_MUTS.txt
  else #If mutations are found:
    echo 'At least one mutation remains, on to annotating final list...' 
    
    if [[ "${Matched,,}" == "true" ]]; then
      # find matched control (strip CR)
      matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable")
    
      if [[ -n "$matched_control" ]]; then
        (
          module load SAMtools/1.19.2-GCCcore-11.3.0
          add_matched_control_metrics "$OutputDir/$samplename/merged_vcfs_annotated-filtered.txt" "$samplename"
        )
        annotate_input="$OutputDir/$samplename/merged_vcfs_annotated-filtered.with_control.txt"
      else
        echo "[MATCHED CTRL] No matched control for $samplename — skipping matched metrics"
        annotate_input="$OutputDir/$samplename/merged_vcfs_annotated-filtered.txt"
      fi
    
      annotate "$annotate_input" "Final_list_$samplename" excel vcf
    else
      annotate "$OutputDir/$samplename/merged_vcfs_annotated-filtered.txt" "Final_list_$samplename" excel vcf
      echo 'tmp'
    fi
  fi
  
  #Wait for MNV algorithm to finish
  wait
  
  #Cleaning:
  mkdir $OutputDir/$samplename/raw-data/ && mkdir $OutputDir/$samplename/filtering-info/ && mkdir $OutputDir/$samplename/intermediate-files/ && mkdir $OutputDir/$samplename/final-list/
  mv $OutputDir/$samplename/filtering_info.txt $OutputDir/$samplename/filtering-info/ && mv $OutputDir/$samplename/removed_variants.xlsx $OutputDir/$samplename/filtering-info/
  mv $OutputDir/$samplename/merged* $OutputDir/$samplename/intermediate-files/ && mv $OutputDir/$samplename/Final_list* $OutputDir/$samplename/final-list/
  
  echo -e "DONE analyzing $samplename!!\n\n\n" 
}

run_samples() {
###################################################
# Running Tumor samples from directory or as single file
mkdir -p "$OutputDir"

if [[ -d "$Inputbam" ]]; then
  if [[ -n "${CrosslinkTable:-}" && -f "$CrosslinkTable" ]]; then
    echo "Using CrosslinkTable: processing baselines --> follow-ups (controls ignored)"

    # Collect baselines (col1 where col2=="baseline"), deduplicate preserving first occurrence
    mapfile -t baselines < <(
      awk -F $'\t' 'NR>1 && $2=="baseline" && !seen[$1]++{print $1}' "$CrosslinkTable"
    )

    # Collect follow-ups (col1 where col2=="follow-up"), deduplicate preserving first occurrence
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
