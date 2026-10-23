#!/bin/bash
#SBATCH --job-name=CircuSNV
#SBATCH --output="CircuSNV.out"
#SBATCH --error="CircuSNV.err"
#SBATCH --time=159:59:00
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --get-user-env
#SBATCH --export=NONE



#TODO :
#Integrate a check on which samples have been processed (see if output file has been generated --> if so: Skip sample). This way, the pipeline can easily continue where it breaks
# Change first annotation step: Less tools to annotate to save time. Only SNPdbs & functional effect (see post filtering steps) that are needed can be kept. 



###  Helper: print usage and exit  ###
usage() {
  cat <<EOF

Usage: bash $0 -I <input_bam_or_dir> -R <reference.fasta> -L <regions.bed> [options]

Required arguments:
  -I   Input BAM file or directory
  -R   Reference FASTA
  -L   BED of target regions

Optional arguments:
  -V   Minimum VAF (default: 0.001)
  -D   Minimum read depth (default: [100])
  -C   Minimum number of variant calls (default: [1])
  -P   Panel-of-Normals directory
  -Q   Minimum base quality (default: [25])
  -M   Minimum reads supporting a mutation (default: [3])
  -S   Strand-bias threshold (default: [0.1])
  -H   Enable MNV calling (true/[false])
  -B   Bayesian Prior for MNVista (default: [0.5])
  -X   Crosslink table TSV
  -T   Use matched control logic (true/[false])
  -U   Use UNMET score to filter error-prone variants (true/[false])
  -O   Main output directory (default: [./output])
  -h   Show this help message and exit

EOF
  exit 1
}

###  Default (optional) parameter values  ###
VAF=0.001            # minimum variant allele frequency
RDP=100              # minimum read depth
Calls=1              # minimum number of variant-calls
PoN=""               # panel-of-normals path (empty = off)
Qual="25"            # minimum base quality
MRD="3"              # minimum reads supporting a mutation
SB="0.9"             # strand-bias threshold
MNV="false"          # enable multi-nucleotide-variant logic?
PRIOR="0.5"          # Bayesian prior for MNVista
CrosslinkTable=""    # TSV file for crosslinking (empty = off)
Matched="false"      # use matched control logic?
UNMET="false"        # use UNMET score filtering
OutputDir="./output"

# number of threads; define or inherit from SLURM_CPUS_PER_TASK
THREADS=${SLURM_CPUS_PER_TASK:-12}

###  Parsing parameters  ###
while getopts "R:L:I:V:D:C:P:Q:M:S:H:B:X:T:U:O:h:" arg; do
  case $arg in
    R) Reference=$OPTARG;;      # "/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly38.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation_hg38.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
    Q) Qual=$OPTARG;;           # BaseQuality (25)
    M) MRD=$OPTARG;;            # Minimum reads supporting a mutation
    S) SB=$OPTARG;;             # Strand bias. Lowest fraction of FW reads vs RV reads. (SB = 0.2 --> FW = 1, RV = 5 and vice versa)
    H) MNV=$OPTARG;;            # if "true", apply MultiNucleotide Variant calling (detect variants in pHase)
    B) PRIOR=$OPTARG;;          # Bayesian prior for calling MNVs with MNVista. Default 0.5 (equal odds)
    X) CrosslinkTable=$OPTARG;; # tab separated table (.tsv file) with 'samplename', 'timepoint' and 'corresponding_baseline' columns
    T) Matched=$OPTARG;;        # if "true", apply matched_control logic
    U) UNMET=$OPTARG;;          # if "true", apply UNMET score to filter hard to sequence regions
    O) OutputDir=$OPTARG;;      # Sets the path to the main directory in which the data will be outputted. 
    h) usage;;
    *) usage;;
  esac
done

#֖ Bail out if any required flags are missing ֖#
check_required_args() {
  local miss=()
  for var in Reference Listregions Inputbam; do
    [[ -z "${!var:-}" ]] && miss+=("$var")
  done
  if (( ${#miss[@]} )); then
    echo "ERROR: Missing required arguments: ${miss[*]}"
    usage
  fi
}

check_required_args

# Initiating settings #
cd /groups/umcg-pmb/tmp02/projects/hematopathology/Nick/CircuSNV/
envpath="./env"
SBU=$(echo "scale=2; 1 / $SB" | bc)

echo -e "\nSettings: \nInput bam(s): $Inputbam\nReference: $Reference\nPanel: $Listregions\nMinimum VAF: $VAF\nMinimum Read Depth: $RDP\nMinimum Base Quality: $Qual\nMinimum Overlapping Calls: $Calls\nPoN: $PoN\nStrand Bias: $SB\nMultiNucleotide Variant Calling: $MNV\nBayesian Prior for MNV calling: $PRIOR\nCrosslink Table: $CrosslinkTable\nMatched Control Logic: $Matched\nUNMET & repeat filtering: $UNMET\nMain output directory: $OutputDir"


# Function to load modules
load_modules() {
  module load zlib/1.2.12-GCCcore-11.3.0
  module load HTSlib/1.19.1-GCCcore-11.3.0
  module load BCFtools/1.19-GCCcore-11.3.0
}

# Function to run Python script
pythonscript() {
  source "$envpath/bin/activate"
  python3 "$@"
  deactivate
}

# Function to annotate with OpenCravat
annotate() { #1 input.txt (tab delim snv data)   #2 Output filenames    #3&#4 output 
  echo -e "annotating $1 with OpenCravat..."
  source "$envpath/bin/activate"
  oc run $1 \
    -l hg38 \
    -n $2 --silent --mp $THREADS \
    -a cgc cadd_exome clinpred cosmic mutation_assessor thousandgenomes gnomad3 repeat \
    -t $3 $4
  deactivate 
}

# Function to normalize, atomize, compress, and index filtered VCF files
process_vcf() {
  local input_vcf=$1
  local output_vcf=$2

  bcftools norm -O v -m -any -f "$Reference" -a --atom-overlaps . -c w "$input_vcf" 2>/dev/null | bcftools norm -d all > "$output_vcf" 2>/dev/null
  bgzip "$output_vcf"
  bcftools index "${output_vcf}.gz"
}

# Function to add UNMET score and filter 
unmet_filter() {
  local samp="$1"                             # sample name, e.g. “SAMPLE123”
  local curr_file="$2"                        # merged_vcfs.txt for that sample
  local bed_file="$3"                         # path to the .bed you want
  local threshold="$4"                        # numeric cutoff, e.g. 0.93
  local label="$5"                            # short label for logging
  local tmp_bed_tdf="${bed_file%.bed}.tdf"
  local out_with="${curr_file%.txt}_UNMET.txt"
  

  # 1) Make sure the UNMET.bed exists
  if [[ ! -f "$bed_file" ]]; then
    echo "UNMET is true and $bed_file not found — downloading + converting…"
    wget -q "http://web16.kazusa.or.jp/data/hg38/$(basename "$tmp_bed_tdf")"
    (
      module load IGV
      igvtools tdftobedgraph "$tmp_bed_tdf" "$bed_file" 2>/dev/null
    )
    echo "UNMET: created $bed_file"
  fi

  # 2) Attach UNMET scores
  echo "[$samp] $label: variants before = $(wc -l < "$curr_file")"
  (
    module load BEDTools
    bedtools intersect -a <(
      awk 'BEGIN{FS=OFS="\t"} { print $1, $2, $2+1, $0 }' "$curr_file"
    ) -b "$bed_file" -loj \
    | awk -v OFS="\t" -v th="$threshold" '
        {
          sc = ($12 == "." ? "NA" : $12)
          printf "%s\t%s\t%s\t%s\t%s_UNMET:%s\n",
                 $4, $5, $6, $7, $8, sc
        }' > "$out_with"
  )
  # 3) Filter by threshold
  awk -F'\t' -v th="$threshold" '
    {
      split($5, a, "_UNMET:")
      sc = a[2]
      if (sc == "NA" || sc+0 < th) print
    }
  ' "$out_with" > "$curr_file"

  echo "[$samp] $label: variants after  = $(wc -l < "$curr_file")"
}

# Function to add reads info from matched control to each mutation found in tumor sample 
add_matched_control_metrics() {
  local mutations="$1"
  local samplename="$2"

  echo "[MATCHED CTRL] Adding matched control metrics for each mutation..."

  # find matched control (strip CR if present)
  matched_control=$(
    awk -F $'\t' -v s="$samplename" '$1==s{gsub(/\r/,"",$4); print $4; exit}' "$CrosslinkTable"
  )
  [[ -z "$matched_control" ]] && { echo "ERROR: no matched_control for $samplename" >&2; return 1; }

  # control BAM path
  ctrl="${Inputbam%/}/${matched_control}.bam"
  [[ ! -f "$ctrl" ]] && { echo "ERROR: control BAM not found: $ctrl" >&2; return 1; }

  # ensure BAM index exists (create if necessary)
  if [[ ! -f "${ctrl}.bai" && ! -f "${ctrl%.bam}.bai" ]]; then
    echo "[MATCHED CTRL] Index for $ctrl not found; creating..."
    samtools index "$ctrl" || { echo "ERROR: samtools index failed for $ctrl" >&2; return 1; }
  fi

  local out="${mutations%.txt}.control_reads.txt"
  local merged="${mutations%.txt}.with_control.txt"

  # make temporary files
  local bedfile
  local pileupfile
  bedfile=$(mktemp "${TMPDIR:-/tmp}/ctrl_pos.XXXXXX.bed")
  pileupfile=$(mktemp "${TMPDIR:-/tmp}/ctrl_pileup.XXXXXX.txt")

  # cleanup on exit
  trap 'rm -f "$bedfile" "$pileupfile"' RETURN INT TERM

  # create BED with 0-based start, 1-based end for each mutation position
  # (also remove potential CRs)
  awk -F $'\t' 'BEGIN{OFS="\t"} {gsub(/\r/,"",$0); if(NF>=2 && $1!="" && $2!="") print $1, $2-1, $2}' "$mutations" > "$bedfile"

  # run mpileup once for all positions (-l accepts a BED), write stdout to pileupfile
  # use -A -B to include anomalous reads & disable BAQ (same as before)
  echo "[MATCHED CTRL] Running samtools mpileup for all positions..."
  samtools mpileup -f "$Reference" -l "$bedfile" -d 10000000 -A -B "$ctrl" > "$pileupfile" 2>/dev/null || true

  # Now parse the pileup and the mutation file in one awk run:
  #  - first (NR==FNR) load pileup lines into an array keyed by chr:pos
  #  - then, for each mutation line, lookup the pileup entry and compute counts (ref/alt/depth)
  # The counting logic is similar to your previous awk, but runs only once per pileup line.
  awk -v OFS="\t" -v CONTROL="$matched_control" -v OUT="$out" -v MERGED="$merged" '
    function strip_start_end(s,   tmp) {
      # remove start-of-read markers ^X and read-end $ from the mpileup read_bases string
      gsub(/\^./, "", s)
      gsub(/\$/, "", s)
      return s
    }

    # count function: s is read_bases, REF and ALT provided
    function count_ref_alt(s, REF, ALT,   i, L, c, n_ref, n_alt, ref_len, alt_len, ins_seq, del_seq, sign, j, num, indlen, seq, base) {
      n_ref = 0; n_alt = 0
      ref_len = length(REF); alt_len = length(ALT)
      if (alt_len > ref_len) ins_seq = substr(ALT, ref_len+1)
      else ins_seq = ""
      if (ref_len > alt_len) del_seq = substr(REF, alt_len+1)
      else del_seq = ""

      L = length(s); i = 1
      while (i <= L) {
        c = substr(s, i, 1)
        if (c == "." || c == ",") {
          # check following token: could be indel marker +/-
          nextc = substr(s, i+1, 1)
          if (nextc == "+" || nextc == "-") {
            sign = nextc; j = i + 2; num = ""
            while (j <= L && substr(s, j, 1) ~ /[0-9]/) { num = num substr(s, j, 1); j++ }
            indlen = (num == "" ? 0 : num + 0)
            seq = substr(s, j, indlen)
            if (sign == "+") {
              if (ins_seq != "" && seq == ins_seq) n_alt++
            } else {
              if (del_seq != "" && seq == del_seq) n_alt++
            }
            i = j + indlen
          } else {
            n_ref++; i++
          }
        } else if (c == "+" || c == "-") {
          sign = c; j = i + 1; num = ""
          while (j <= L && substr(s, j, 1) ~ /[0-9]/) { num = num substr(s, j, 1); j++ }
          indlen = (num == "" ? 0 : num + 0)
          seq = substr(s, j, indlen)
          if (sign == "+") {
            if (ins_seq != "" && seq == ins_seq) n_alt++
          } else {
            if (del_seq != "" && seq == del_seq) n_alt++
          }
          i = j + indlen
        } else if (c == "*") {
          # placeholder for deletion in mpileup; ignore
          i++
        } else {
          base = toupper(c)
          if (ref_len == 1 && alt_len == 1) {
            if (base == toupper(ALT)) n_alt++
            else if (base == toupper(REF)) n_ref++
            else n_alt++    # conservative: other base = alt
          } else {
            # heuristic for indel/complex: match first base of ALT
            if (alt_len > 0 && base == toupper(substr(ALT,1,1))) n_alt++
            else n_ref++
          }
          i++
        }
      }
      return n_ref " " n_alt
    }

    # --- load pileup into array keyed by chr:pos ---
    FNR==NR {
      # pileup format: chr pos ref depth read_bases qual...
      chr = $1; pos = $2; depth = $4; read_bases = $5
      key = chr ":" pos
      read_bases = strip_start_end(read_bases)
      counts = count_ref_alt(read_bases, "X", "Y")  # placeholder; we will recount per REF/ALT later
      # store raw read_bases and depth for per-mutation counting later
      pile_rb[key] = read_bases
      pile_dp[key] = depth + 0
      next
    }

    # --- process mutations file lines ---
    {
      line = $0
      # preserve original fields reliably even if info contains tabs (we assume tabs separate the standard 6 columns)
      chrom = $1; pos = $2; ref = $3; alt = $4; sample = $5
      # info might be empty or could contain tabs; reconstruct it starting at $6 if present
      info = ""
      if (NF >= 6) {
        info = $6
        for (i = 7; i <= NF; i++) info = info OFS $i
      }

      key = chrom ":" pos
      if (key in pile_rb) {
        rb = pile_rb[key]
        dp = pile_dp[key] + 0
        # compute counts using actual REF/ALT
        cl = count_ref_alt(rb, ref, alt)
        split(cl, a, " ")
        ref_count = a[1] + 0
        alt_count = a[2] + 0
      } else {
        ref_count = 0; alt_count = 0; dp = 0
      }

      # write control reads output (control name only)
      print chrom, pos, ref, alt, sample, CONTROL, ref_count, alt_count, dp > OUT

      # build merged info: append CONTROL and counts to original info
      ctrl_append = "_CONTROL:" CONTROL "_CREF:" ref_count "_CALT:" alt_count "_CDP:" dp
      if (info == "") new_info = ctrl_append
      else new_info = info ctrl_append

      # write merged line: keep six-columns-like shape: chrom pos ref alt sample new_info
      print chrom, pos, ref, alt, sample, new_info > MERGED
    }
  ' "$pileupfile" "$mutations"

  echo "[MATCHED CTRL] Wrote control reads to: $out"
  echo "[MATCHED CTRL] Wrote merged input+control to: $merged"

  # cleanup (trap will remove tmp files)
  return 0
}




################################################### 
# Functions to run each tool, filter on Read Dept, VAF and Mutant allele dept & creating .vcf file

Vardict() { #1: Input .bam file   #2 Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with VarDict..."
  (
    module load R/4.2.2-foss-2022a-bare
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
    module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
    load_modules
    
    gatk --java-options -Xmx15G Mutect2 --QUIET --verbosity ERROR --create-output-variant-index false --native-pair-hmm-threads $THREADS --min-base-quality-score $Qual --callable-depth 10 --minimum-allele-fraction 0.00000000001 -R $Reference -I $1 -O $2${base_name}_Mutect2_RAW.vcf -L $Listregions $3 && \
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
    module load SAMtools/1.19.2-GCCcore-11.3.0
    
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
    module load HTSlib/1.19.1-GCCcore-11.3.0
    sinvict -m 10 -t $2/${base_name}_SiNVICT_inter_bamreadcount/ -o $2
  )
  
  (
    module load SAMtools/1.19.2-GCCcore-11.3.0
    pythonscript ./scripts/sinvict_to_vcf2.py $2/calls_level1.sinvict $2${base_name}_SiNVICT_inter.vcf $Reference
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
################################################### 

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
    pythonscript ./scripts/merge_variants.py $samplename "phased" --infile "sites.txt" --outfile "sites2.txt" $OutputDir
    
    #### MNV Filtering on UNMET score ###
    if [[ "$UNMET" == "true" ]]; then
      unmet_filter "$samplename" "$mnvpath/sites2.txt" "./hg38_genome.UNMET_bin10bp.bed" 0.90 "MNV input -- UNMET filtering"
    fi  
    
    #convert to vcf
    pythonscript ./scripts/isec_to_vcf.py $mnvpath/sites2.txt $mnvpath/sites.vcf 'TRUE'
    bgzip -c $mnvpath/sites.vcf > $mnvpath/sites.vcf.gz && bcftools index $mnvpath/sites.vcf.gz
    echo "MNV input: merging 4 tools to 1 list DONE --> sites.vcf" 
  )
}


mnv_calling() {
  local samplename=$1
  local crosslink_table=$2
  
  mnvista_cmd="./MNVista $bam_file $OutputDir/$samplename/mnvcalling/sites.vcf $OutputDir/$samplename/mnvcalling/results -T 15 -M 6 -R 150 -F 0 -J 0"
  [ -n "$PoN" ] && mnvista_cmd+=" -C './MNV_PoN_No_overlap_new.txt'"
  
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
        module load HTSlib/1.19.1-GCCcore-11.3.0
        ml GCC/13.3.0
        mnvista_cmd+=" -S 1 -N 1 -B 0.00000000001 -P 1"
        eval "$mnvista_cmd"
        echo -e "MNV process for sample $samplename DONE ...\n"
      )
    
    elif [ "$timepoint" == "baseline" ]; then
      # Sample is baseline
      echo "MNV input: Sample $samplename is a baseline sample."
      generate_mnv_input
      (
        module load HTSlib/1.19.1-GCCcore-11.3.0
        ml GCC/13.3.0 
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
              ./MNVista "$ctrl" "$OutputDir/$samplename/mnvcalling/results/results.vcf" "$matched_output" -T 12 -M 6 -R 100 -F 0 -J 0 -S 10 -N 3 -B 0.95 -P $PRIOR -K true
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
      module load HTSlib/1.19.1-GCCcore-11.3.0
      ml GCC/13.3.0
      mnvista_cmd+=" -S 10 -N 3 -B 0.95 -P $PRIOR -K true"
      eval "$mnvista_cmd"
      echo -e "MNV process for sample $samplename DONE ...\n"
    )
  fi
}



################################################### 
# PON INITIATION: Create PoN if argument is given:
if [ -z "$PoN" ]
then
  #No PoN Argument
  echo -e "\nNo PoN argument given\n"
else
  # If PoN Argument is given:
  echo -e "\nPoN argument is given!!\n"
  mkdir -p ./PoN/
  
  
  #Directory PoN with containing samples:
  if [[ -d $PoN ]]; then
    echo 'PoN is directory'
    echo -e "\nPoN included: source for PoN .bam files = $PoN"
    for healthycontrol in "$PoN"/*.bam
    do
      samplename=$(basename "$healthycontrol" .bam)
      bam_dir=$(dirname "$healthycontrol")
      bam_index=(${bam_dir}/${samplename}*.bai)
      
      echo -e "\n\nGenerating VCF for normal control: $samplename"
      mkdir ./PoN/$samplename/
      mkdir -p ./PoN/$samplename/mnvcalling
      
      # Check if BAM index file exists, if not create it
      if [ ! -e "${bam_index[0]}" ]; then
        echo "BAM index file for ${samplename} not found. Creating BAM index..."
        (
          module load SAMtools/1.16.1-GCCcore-11.3.0
          samtools index "$healthycontrol"
        )
      else
        echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
        touch $bam_index
      fi

      Vardict $healthycontrol ./PoN/$samplename/ & Mutect2 $healthycontrol ./PoN/$samplename/ & Lofreq $healthycontrol ./PoN/$samplename/ & Sinvict $healthycontrol ./PoN/$samplename/
      wait
    done
    
    echo -e "\nMerging outputs for PoN..."
    #Merge per tool -> 4 seperate blacklists as .vcf file
    ls ./PoN/*/*Mutect2.vcf.gz | sed 's/^/ /' > ./PoN/Mutect2-normals.dat & ls ./PoN/*/*VarDict.vcf.gz | sed 's/^/ /' > ./PoN/VarDict-normals.dat & \
    ls ./PoN/*/*LoFreq.vcf.gz | sed 's/^/ /' > ./PoN/LoFreq-normals.dat & ls ./PoN/*/*SiNVICT.vcf.gz | sed 's/^/ /' > ./PoN/SiNVICT-normals.dat
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      xargs -a ./PoN/Mutect2-normals.dat bcftools isec -o ./PoN/PoN_Mutect2.txt -O v -n +1 && xargs -a ./PoN/VarDict-normals.dat bcftools isec -o ./PoN/PoN_Vardict.txt -O v -n +1 && \
      xargs -a ./PoN/LoFreq-normals.dat bcftools isec -o ./PoN/PoN_Lofreq.txt -O v -n +1 && xargs -a ./PoN/SiNVICT-normals.dat bcftools isec -o ./PoN/PoN_Sinvict.txt -O v -n +1
    )
    
    pythonscript ./scripts/isec_to_vcf.py ./PoN/PoN_Mutect2.txt ./PoN/PoN_Mutect2.vcf 'FALSE' && pythonscript ./scripts/isec_to_vcf.py ./PoN/PoN_Vardict.txt ./PoN/PoN_Vardict.vcf 'FALSE' && \
    pythonscript ./scripts/isec_to_vcf.py ./PoN/PoN_Lofreq.txt ./PoN/PoN_Lofreq.vcf 'FALSE' && pythonscript ./scripts/isec_to_vcf.py ./PoN/PoN_Sinvict.txt ./PoN/PoN_Sinvict.vcf 'FALSE'
    
    (
      #zipping and indexing Pon-per-tool blacklists for filtering later on;
      load_modules
      bgzip -c ./PoN/PoN_Mutect2.vcf > ./PoN/PoN_Mutect2.vcf.gz && bgzip -c ./PoN/PoN_Vardict.vcf > ./PoN/PoN_Vardict.vcf.gz && \
      bgzip -c ./PoN/PoN_Lofreq.vcf > ./PoN/PoN_Lofreq.vcf.gz && bgzip -c ./PoN/PoN_Sinvict.vcf > ./PoN/PoN_Sinvict.vcf.gz

      bcftools index ./PoN/PoN_Mutect2.vcf.gz && bcftools index ./PoN/PoN_Vardict.vcf.gz && bcftools index ./PoN/PoN_Lofreq.vcf.gz && bcftools index ./PoN/PoN_Sinvict.vcf.gz
      #We now have 4 PoNs, 1 per tool
      
      #Merging into 1 final blacklist...
      bcftools isec -o ./PoN/BLACKLIST_combined.txt -O v -n +2 ./PoN/PoN_Sinvict.vcf.gz ./PoN/PoN_Mutect2.vcf.gz ./PoN/PoN_Lofreq.vcf.gz ./PoN/PoN_Vardict.vcf.gz
    )
      
      #Add mutations from each tool list (N=4) that occur in over 50% of normal samples, to BLACKLIST_combined.txt 
      # Meaning: Inspect PoN_Vardict.txt etc etc, add high-occurance mutations that are not already present in BLACKLIST_combined.txt to BLACKLIST_combined.txt and re-sort.  
      echo -e "Adding mutations to PoN that are found in 20% of normals, but only called with one tool \n$PoN \n"
      pythonscript ./scripts/finish_pon.py 20 ./PoN/

      pythonscript ./scripts/isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      load_modules
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
    echo -e "PoN generation is finished! On to the tumor samples... \n"   
  
  #One-File: a.k.a pre-made blacklist:      
  elif [[ -f $PoN ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp $PoN ./PoN/BLACKLIST_combined_complete.txt
    
    pythonscript ./scripts/isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      load_modules
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
    echo -e "PoN preparation is finished! On to the tumor samples... \n"   

  #Else (error)
  else
    echo "Invalid PoN input. Input folder with .bam files from healthy controls or a pre-generated PoN.txt "
    exit 1
  fi
fi

#PoN generation finished
###################################################



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

###


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
    pythonscript ./scripts/isec_to_vcf.py $OutputDir/$samplename/isec/sites.txt $OutputDir/$samplename/isec/sites.vcf 'FALSE'
    bgzip -c $OutputDir/$samplename/isec/sites.vcf > $OutputDir/$samplename/isec/sites.vcf.gz && bcftools index $OutputDir/$samplename/isec/sites.vcf.gz
    #PoN filtering (or not)
    if [ -z "$PoN" ]; then
      echo -e "No PoN argument given, proceeding with intersect results..."
      echo "Adding mutation metrics back to merged mutation list"
      pythonscript ./scripts/merge_variants.py $samplename "no_pon" $OutputDir
    else
      #Filter combined list with combined blacklist:
      echo -e "PoN argument given, filtering intersect results with PoN..."
      bcftools isec -p $OutputDir/$samplename/isec/PoN-filtered -O v -n~10 $OutputDir/$samplename/isec/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
      
      # Check if variants are left after PoN Filtering ##########################################      
      
      #Add VAF/RD/MRD data to the merged list:
      echo "Adding mutation metrics back to pon-filtered merged mutation list"
      pythonscript ./scripts/merge_variants.py $samplename "pon" $OutputDir
    fi
  )  
  
  #### Filtering on UNMET score ###
  if [[ "$UNMET" == "true" ]]; then
    unmet_filter "$samplename" "$OutputDir/$samplename/merged_vcfs.txt" "./hg38_genome.UNMET_bin10bp.bed" 0.90 "UNMET filtering"
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
  pythonscript ./scripts/post_filtering.py $samplename $OutputDir $UNMET
    
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

# Function for tumor-sample finished
###################################################



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

# Running Tumor samples finished
###################################################

echo 'Finished run!'





