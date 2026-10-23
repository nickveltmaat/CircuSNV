# shellcheck shell=bash
# Panel of Normals: generate from normal BAMs or prepare a pre-made blacklist.
# Sourced by CircuSNV.sh; uses its global settings.
#
# The result is $PON_DIR/BLACKLIST_combined.vcf.gz, used to filter the tumor samples.

prepare_pon() {
  local healthycontrol samplename bam_dir
  local bam_index=()

  if [ -z "$PoN" ]; then
    echo -e "\nNo PoN argument given\n"
    return 0
  fi

  echo -e "\nPoN argument is given!!\n"
  mkdir -p "$PON_DIR/"

  # Directory with normal BAMs: call variants in every normal and build the blacklist
  if [[ -d "$PoN" ]]; then
    echo 'PoN is directory'
    echo -e "\nPoN included: source for PoN .bam files = $PoN"
    for healthycontrol in "$PoN"/*.bam
    do
      samplename=$(basename "$healthycontrol" .bam)
      bam_dir=$(dirname "$healthycontrol")
      bam_index=("${bam_dir}/${samplename}"*.bai)

      echo -e "\n\nGenerating VCF for normal control: $samplename"
      mkdir "$PON_DIR/$samplename/"
      mkdir -p "$PON_DIR/$samplename/mnvcalling"

      # Create the BAM index if it does not exist
      if [ ! -e "${bam_index[0]}" ]; then
        echo "BAM index file for ${samplename} not found. Creating BAM index..."
        (
          module load "$MODULE_SAMTOOLS_PON"
          samtools index "$healthycontrol"
        )
      else
        echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
        touch "${bam_index[0]}"
      fi

      Vardict "$healthycontrol" "$PON_DIR/$samplename/" & Mutect2 "$healthycontrol" "$PON_DIR/$samplename/" & Lofreq "$healthycontrol" "$PON_DIR/$samplename/" & Sinvict "$healthycontrol" "$PON_DIR/$samplename/"
      wait
    done

    echo -e "\nMerging outputs for PoN..."
    # Merge per tool -> 4 separate blacklists as .vcf file
    ls "$PON_DIR/"*/*Mutect2.vcf.gz | sed 's/^/ /' > "$PON_DIR/Mutect2-normals.dat" & ls "$PON_DIR/"*/*VarDict.vcf.gz | sed 's/^/ /' > "$PON_DIR/VarDict-normals.dat" & \
    ls "$PON_DIR/"*/*LoFreq.vcf.gz | sed 's/^/ /' > "$PON_DIR/LoFreq-normals.dat" & ls "$PON_DIR/"*/*SiNVICT.vcf.gz | sed 's/^/ /' > "$PON_DIR/SiNVICT-normals.dat"
    (
      module load "$MODULE_BCFTOOLS"
      xargs -a "$PON_DIR/Mutect2-normals.dat" bcftools isec -o "$PON_DIR/PoN_Mutect2.txt" -O v -n +1 && xargs -a "$PON_DIR/VarDict-normals.dat" bcftools isec -o "$PON_DIR/PoN_Vardict.txt" -O v -n +1 && \
      xargs -a "$PON_DIR/LoFreq-normals.dat" bcftools isec -o "$PON_DIR/PoN_Lofreq.txt" -O v -n +1 && xargs -a "$PON_DIR/SiNVICT-normals.dat" bcftools isec -o "$PON_DIR/PoN_Sinvict.txt" -O v -n +1
    )

    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/PoN_Mutect2.txt" "$PON_DIR/PoN_Mutect2.vcf" 'FALSE' && pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/PoN_Vardict.txt" "$PON_DIR/PoN_Vardict.vcf" 'FALSE' && \
    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/PoN_Lofreq.txt" "$PON_DIR/PoN_Lofreq.vcf" 'FALSE' && pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/PoN_Sinvict.txt" "$PON_DIR/PoN_Sinvict.vcf" 'FALSE'

    (
      # Compress and index the per-tool blacklists
      load_modules
      bgzip -c "$PON_DIR/PoN_Mutect2.vcf" > "$PON_DIR/PoN_Mutect2.vcf.gz" && bgzip -c "$PON_DIR/PoN_Vardict.vcf" > "$PON_DIR/PoN_Vardict.vcf.gz" && \
      bgzip -c "$PON_DIR/PoN_Lofreq.vcf" > "$PON_DIR/PoN_Lofreq.vcf.gz" && bgzip -c "$PON_DIR/PoN_Sinvict.vcf" > "$PON_DIR/PoN_Sinvict.vcf.gz"

      bcftools index "$PON_DIR/PoN_Mutect2.vcf.gz" && bcftools index "$PON_DIR/PoN_Vardict.vcf.gz" && bcftools index "$PON_DIR/PoN_Lofreq.vcf.gz" && bcftools index "$PON_DIR/PoN_Sinvict.vcf.gz"

      # Variants found by at least 2 tools form the combined blacklist
      bcftools isec -o "$PON_DIR/BLACKLIST_combined.txt" -O v -n +2 "$PON_DIR/PoN_Sinvict.vcf.gz" "$PON_DIR/PoN_Mutect2.vcf.gz" "$PON_DIR/PoN_Lofreq.vcf.gz" "$PON_DIR/PoN_Vardict.vcf.gz"
    )

    # Add variants found by a single tool in at least PON_MIN_NORMALS_PCT % of the normals
    echo -e "Adding mutations to PoN that are found in ${PON_MIN_NORMALS_PCT}% of normals, but only called with one tool \n$PoN \n"
    pythonscript "$CIRCUSNV_HOME/scripts/finish_pon.py" "$PON_MIN_NORMALS_PCT" "$PON_DIR/"

    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/BLACKLIST_combined_complete.txt" "$PON_DIR/BLACKLIST_combined.vcf" 'FALSE' && \
    (
      load_modules
      bgzip -c "$PON_DIR/BLACKLIST_combined.vcf" > "$PON_DIR/BLACKLIST_combined.vcf.gz"
      bcftools index "$PON_DIR/BLACKLIST_combined.vcf.gz"
    )
    echo -e "PoN generation is finished! On to the tumor samples... \n"

  # Single file: a pre-made blacklist
  elif [[ -f "$PoN" ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp "$PoN" "$PON_DIR/BLACKLIST_combined_complete.txt"

    pythonscript "$CIRCUSNV_HOME/scripts/isec_to_vcf.py" "$PON_DIR/BLACKLIST_combined_complete.txt" "$PON_DIR/BLACKLIST_combined.vcf" 'FALSE' && \
    (
      load_modules
      bgzip -c "$PON_DIR/BLACKLIST_combined.vcf" > "$PON_DIR/BLACKLIST_combined.vcf.gz"
      bcftools index "$PON_DIR/BLACKLIST_combined.vcf.gz"
    )
    echo -e "PoN preparation is finished! On to the tumor samples... \n"

  else
    echo "Invalid PoN input. Input folder with .bam files from healthy controls or a pre-generated PoN.txt "
    exit 1
  fi
}
