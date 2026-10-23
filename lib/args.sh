# shellcheck shell=bash
# Command-line usage, default settings and argument parsing.
# Sourced by CircuSNV.sh; uses its global settings.

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

parse_args() {
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
}

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
