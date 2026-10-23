# shellcheck shell=bash
# Variant callers (VarDict, Mutect2, LoFreq, SiNVICT), VCF normalisation and the retry wrapper.
# Sourced by CircuSNV.sh; uses its global settings.

# Function to normalize, atomize, compress, and index filtered VCF files
process_vcf() {
  local input_vcf=$1
  local output_vcf=$2

  bcftools norm -O v -m -any -f "$Reference" -a --atom-overlaps . -c w "$input_vcf" 2>/dev/null | bcftools norm -d all > "$output_vcf" 2>/dev/null
  bgzip "$output_vcf"
  bcftools index "${output_vcf}.gz"
}

################################################### 
# Functions to run each tool, filter on Read Dept, VAF and Mutant allele dept & creating .vcf file

Vardict() { #1: Input .bam file   #2 Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with VarDict..."
  (
    module load $MODULE_R
    load_modules
    
    vardict-java -th $THREADS -r 1 -G $Reference -f 0.00000000001 -N $base_name -I 20 -L 20 -Y 20 --nosv --deldupvar -q $Qual -Q $Qual -b $1 -c 1 -S 2 -E 3 -g 4 $Listregions | \
    teststrandbias.R | var2vcf_valid.pl -S -N $base_name -E -f 0.00000000001 > $2${base_name}_VarDict_RAW_inter.vcf
    
    bcftools filter -O v -i "(FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) >= ${SB} & (FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) <= ${SBU} & (INFO/ODDRATIO[0]) > 0 & INFO/PMEAN < 67" \
    $2${base_name}_VarDict_RAW_inter.vcf 2>/dev/null > $2${base_name}_VarDict_RAW.vcf
    
    bgzip -c $2${base_name}_VarDict_RAW.vcf > $2${base_name}_VarDict_inter.vcf.gz && bcftools index $2${base_name}_VarDict_inter.vcf.gz && cp $2${base_name}_VarDict_inter.vc* $2/mnvcalling
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/VD > ${MRD}" $2${base_name}_VarDict_inter.vcf.gz > $2${base_name}_VarDict_inter2.vcf.gz && \
    process_vcf "$2${base_name}_VarDict_inter2.vcf.gz" "$2${base_name}_VarDict.vcf"
    rm -rf $2${base_name}*inter*.*
    echo -e "VarDict DONE for $base_name\n"
  )
}

######
Mutect2() { #1: Input .bam file    #2: Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with Mutect2..."
  (
    module load $MODULE_GATK
    load_modules
    
    gatk --java-options -Xmx$MUTECT2_JAVA_MEM Mutect2 --QUIET --verbosity ERROR --create-output-variant-index false --native-pair-hmm-threads $THREADS --min-base-quality-score $Qual --callable-depth 10 --minimum-allele-fraction 0.00000000001 -R $Reference -I $1 -O $2${base_name}_Mutect2_RAW.vcf -L $Listregions $3 && \
    gatk FilterMutectCalls --QUIET --verbosity ERROR --create-output-variant-index false --min-allele-fraction 0.00000000001 -R $Reference -V $2${base_name}_Mutect2_RAW.vcf -O $2${base_name}_M2_Fstats1.vcf
    bcftools filter -O v -i "(FORMAT/SB[0:2]/FORMAT/SB[0:3]) >= ${SB} & (FORMAT/SB[0:2]/FORMAT/SB[0:3]) <= ${SBU}" $2${base_name}_M2_Fstats1.vcf > $2${base_name}_M2_Fstats.vcf
    
    grep -v '^##contig=<ID=chrUn\|^##contig=<ID=HLA' $2${base_name}_M2_Fstats.vcf | grep -v '_alt' | grep -v '_random' > $2${base_name}_M2_F2stats.vcf && cp $2${base_name}_M2_F2stats.vc* $2/mnvcalling
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/AD[0:1] > ${MRD}" $2${base_name}_M2_F2stats.vcf > $2${base_name}M2_F3stats.vcf.gz && \
    process_vcf "$2${base_name}M2_F3stats.vcf.gz" "$2${base_name}_Mutect2.vcf"
    
    rm -rf $2${base_name}*tats.*
    echo -e "Mutect2 DONE for $base_name\n"
  )
}


######
Lofreq() { #1: Input .bam file    #2: Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with LoFreq..."
  (
    load_modules
    module load $MODULE_SAMTOOLS
    
    # create a temporary BAM of just the regions you care about
    tmpbam=$(mktemp --suffix=.bam "$2/${base_name}_regions_XXXX")
    samtools view -b -M -L "$Listregions" "$1" -o "$tmpbam"
    samtools index "$tmpbam"

    # 1) run one multi-threaded LoFreq on the subset BAM
    lofreq call-parallel \
      --pp-threads $THREADS \
      --no-default-filter \
      --call-indels \
      -C 10 \
      -f "$Reference" \
      -q "$Qual" -Q "$Qual" \
      -o "$2/${base_name}_LoFreq_inter.vcf" \
      "$tmpbam" && \

    # 2) dedupe
    awk '!seen[$0]++' "$2/${base_name}_LoFreq_inter.vcf" \
      > "$2/${base_name}_LoFreq_inter_RAW.vcf" && \

    # 3) apply your DP4‐ratio filter
    bcftools filter -O v \
      -i "(INFO/DP4[2]/INFO/DP4[3]) >= ${SB} & (INFO/DP4[2]/INFO/DP4[3]) <= ${SBU}" \
      "$2/${base_name}_LoFreq_inter_RAW.vcf" \
      > "$2/${base_name}_LoFreq_RAW.vcf" && \

    # 4) compress & index
    bgzip -c "$2/${base_name}_LoFreq_RAW.vcf" \
      > "$2/${base_name}_LoFreq_inter.vcf.gz" && \
    bcftools index "$2/${base_name}_LoFreq_inter.vcf.gz" && \
    cp "$2/${base_name}_LoFreq_inter.vcf."* "$2/mnvcalling" && \

    # 5) final DP+AF+MRD filter
    bcftools filter -O z \
      -i "INFO/DP > ${RDP} & INFO/AF > ${VAF} & INFO/DP4[2]+INFO/DP4[3] > ${MRD}" \
      "$2/${base_name}_LoFreq_inter.vcf.gz" \
      > "$2/${base_name}_LoFreq_inter2.vcf.gz" && \
    process_vcf "$2/${base_name}_LoFreq_inter2.vcf.gz" \
                "$2/${base_name}_LoFreq.vcf" && \

   rm -f "$tmpbam" "$tmpbam.bai"
   rm -rf $2${base_name}_LoFreq_inter*.*
   echo -e "LoFreq DONE for $base_name\n"
  )
}

######
Sinvict() { #1: Input .bam file    #2: Output folder
  InputSinvict=$1
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with bam-readcount && SiNVICT..."
  
  mkdir -p $2/${base_name}_SiNVICT_inter_bamreadcount/
  
  # Define the number of chunks for parallel processing
  N=$THREADS
  
  # Split the BED file into N chunks (excluding the header lines) and add the header lines back to each chunk
  header_lines=$(grep -P "^browser|^track" $Listregions | wc -l)
  tail -n +$((header_lines + 1)) $Listregions | split -l $((($(wc -l < $Listregions) - header_lines) / N + 1)) -d -a 2 --additional-suffix=.bed - $2/bed_chunk_

  for chunk_file in $2/bed_chunk_*.bed; do
    (grep -P "^browser|^track" $Listregions; cat $chunk_file) > ${chunk_file}.tmp && mv ${chunk_file}.tmp $chunk_file
  done
  
  # Function to run bam-readcount for a chunk
  run_bam_readcount() {
    local chunk_file=$1
    local output_file=$2
    load_modules
    bam-readcount -l $chunk_file -w 1 -b $Qual -f $Reference $InputSinvict 2>/dev/null > $output_file
  }
  
  # Run bam-readcount in parallel using background jobs
  (
    for chunk_file in $2/bed_chunk_*.bed; do
      run_bam_readcount $chunk_file $2/${base_name}_SiNVICT_inter_bamreadcount/output_$(basename $chunk_file .bed).bamreadcount &
    done
    wait
    echo 'bamreadcount finished'
  )
  
  # Merge bam-readcount output files
  cat $2/${base_name}_SiNVICT_inter_bamreadcount/output_*.bamreadcount > $2/${base_name}_SiNVICT_inter_bamreadcount/output.bamreadcount
  
  # Run SiNVICT with the merged bam-readcount output
  (
    module load $MODULE_HTSLIB
    sinvict -m 10 -t $2/${base_name}_SiNVICT_inter_bamreadcount/ -o $2
  )
  
  (
    module load $MODULE_SAMTOOLS
    pythonscript "$CIRCUSNV_HOME"/scripts/sinvict_to_vcf2.py $2/calls_level1.sinvict $2${base_name}_SiNVICT_inter.vcf $Reference
  )

  (
    load_modules

    bcftools sort -o $2${base_name}_SiNVICT_inter_RAW1.vcf -O v $2${base_name}_SiNVICT_inter.vcf
    bcftools filter -O v -i "(FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) >= ${SB} & (FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) <= ${SBU}" $2${base_name}_SiNVICT_inter_RAW1.vcf > $2${base_name}_SiNVICT_RAW.vcf
    bgzip -c $2${base_name}_SiNVICT_RAW.vcf > $2${base_name}_SiNVICT_inter_RAW.vcf.gz

    bcftools index $2${base_name}_SiNVICT_inter_RAW.vcf.gz && cp $2${base_name}_SiNVICT_inter_RAW.vc* $2/mnvcalling  
    
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/VAF > ${VAF} & FORMAT/MDP > ${MRD}" $2${base_name}_SiNVICT_inter_RAW.vcf.gz > $2${base_name}_SiNVICT_inter3.vcf.gz && \
    process_vcf "$2${base_name}_SiNVICT_inter3.vcf.gz" "$2${base_name}_SiNVICT.vcf"  

    rm -rf $2*SiNVICT_inter* && rm -rf $2*alls_level* rm -rf $2*bed_chunk_*
    echo -e "SiNVICT DONE for $base_name\n"
  )
}

### NEW FUNCTION: 
# Retry wrapper: run a caller and ensure expected files exist, retry if necessary.
# Usage: run_caller_with_retries <CallerName> <bam> <outdir> [max_attempts] [sleep_seconds]
run_caller_with_retries() {
  local caller="$1"; shift
  local bam="$1"; shift
  local outdir="$1"; shift
  local max_attempts="${1:-3}"; shift || true
  local sleep_seconds="${1:-10}"; shift || true

  local base
  base=$(basename "$bam" .bam)

  # determine expected outputs for each caller
  local expected=()
  case "$caller" in
    Vardict)
      expected+=("$outdir/${base}_VarDict_RAW.vcf" "$outdir/${base}_VarDict.vcf.gz.csi")
      ;;
    Mutect2)
      expected+=("$outdir/${base}_Mutect2_RAW.vcf" "$outdir/${base}_Mutect2.vcf.gz.csi")
      ;;
    Lofreq)
      expected+=("$outdir/${base}_LoFreq_RAW.vcf" "$outdir/${base}_LoFreq.vcf.gz.csi")
      ;;
    Sinvict|SiNVICT|SiNVict)
      expected+=("$outdir/${base}_SiNVICT_RAW.vcf" "$outdir/${base}_SiNVICT.vcf.gz.csi")
      ;;
    *)
      echo "run_caller_with_retries: unknown caller '$caller'" >&2
      return 2
      ;;
  esac

  local attempt=1
  local missing
  while (( attempt <= max_attempts )); do
    echo "[${caller}] attempt ${attempt}/${max_attempts} for sample ${base} (outdir=${outdir})"

    # call the function if present, otherwise call as external command
    if declare -F "$caller" >/dev/null 2>&1; then
      # it's a shell function - call it
      if ! "$caller" "$bam" "$outdir"; then
        echo "[${caller}] returned non-zero exit on attempt $attempt"
      fi
    else
      # try as external command
      if ! command -v "$caller" >/dev/null 2>&1; then
        echo "[${caller}] not found as function or external command" >&2
        return 2
      fi
      if ! "$caller" "$bam" "$outdir"; then
        echo "[${caller}] external command returned non-zero exit on attempt $attempt"
      fi
    fi

    # verify expected files exist and are non-empty
    missing=()
    for f in "${expected[@]}"; do
      if [[ ! -s "$f" ]]; then
        missing+=("$f")
      fi
    done

    if (( ${#missing[@]} == 0 )); then
      echo "[${caller}] success for ${base}"
      return 0
    fi

    echo "[${caller}] missing/empty outputs after attempt $attempt: ${missing[*]}"
    if (( attempt < max_attempts )); then
      echo "[${caller}] retrying in ${sleep_seconds}s..."
      sleep "$sleep_seconds"
    fi
    ((attempt++))
  done

  echo "[${caller}] failed after ${max_attempts} attempts for ${base}. Missing: ${missing[*]}" >&2
  return 1
}
