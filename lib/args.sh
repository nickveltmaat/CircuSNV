# shellcheck shell=bash
# Command-line options: defaults, help text, parsing and validation.
# Sourced by CircuSNV.sh; uses its global settings.

# Default values of the command-line options. A config file (--config) may override them,
# and options on the command line override both.
set_default_options() {
  VAF=0.001            # minimum variant allele frequency
  RDP=100              # minimum read depth
  Calls=1              # minimum number of callers
  PoN=""               # panel-of-normals path (empty = off)
  Qual="25"            # minimum base quality
  MRD="3"              # minimum reads supporting a mutation
  SB="0.9"             # strand-bias threshold
  MNV="false"          # multi-nucleotide variant calling
  PRIOR="0.5"          # Bayesian prior for MNVista
  CrosslinkTable=""    # TSV file for crosslinking (empty = off)
  Matched="false"      # matched-control logic
  UNMET="false"        # UNMET score and repeat filtering
  OutputDir="./output"
  Reference=""
  Listregions=""
  Inputbam=""
  CONFIG_FILE=""
  SKIP_CHECKS="false"
  DRY_RUN="false"
  # Number of threads: from SLURM_CPUS_PER_TASK when running under Slurm
  THREADS=${SLURM_CPUS_PER_TASK:-12}
}

circusnv_version() {
  local version="unknown"
  if [[ -f "$CIRCUSNV_HOME/VERSION" ]]; then
    version=$(head -n 1 "$CIRCUSNV_HOME/VERSION")
  fi
  echo "$version"
}

print_help() {
  cat <<EOF
CircuSNV $(circusnv_version) - somatic SNV/indel calling in cfDNA without matched controls

Usage: CircuSNV.sh -I <bam_or_dir> -R <reference.fasta> -L <regions.bed> [options]
       sbatch CircuSNV.sh -I <bam_or_dir> -R <reference.fasta> -L <regions.bed> [options]

Required:
  -I, --input PATH             BAM file, or directory with BAM files
  -R, --reference FILE         reference FASTA (hg38)
  -L, --regions FILE           BED file with target regions

Filters:
  -V, --min-vaf NUM            minimum variant allele frequency        (default: $VAF)
  -D, --min-depth INT          minimum read depth                      (default: $RDP)
  -C, --min-callers INT        minimum number of callers (1-4)         (default: $Calls)
  -Q, --min-base-quality INT   minimum base quality                    (default: $Qual)
  -M, --min-alt-reads INT      minimum reads supporting a mutation     (default: $MRD)
  -S, --strand-bias NUM        strand-bias threshold; the forward/reverse ratio of the
                               supporting reads must lie between NUM and 1/NUM
                                                                       (default: $SB)
  -P, --pon PATH               panel of normals: directory with normal BAMs, or a
                               pre-made blacklist (chrom pos ref alt, no header)
  -U, --unmet true|false       filter on UNMET score and repeat regions (default: $UNMET)

Optional analyses:
  -H, --mnv true|false         multi-nucleotide variant calling with MNVista (default: $MNV)
  -B, --mnv-prior NUM          Bayesian prior for MNVista              (default: $PRIOR)
  -X, --crosslink-table FILE   TSV: samplename, timepoint (baseline/follow-up),
                               corresponding_baseline, matched_control
  -T, --matched-control true|false
                               add read counts from the matched control (needs -X)
                                                                       (default: $Matched)

Run:
  -O, --output DIR             main output directory                   (default: $OutputDir)
      --threads INT            threads per tool (default: SLURM_CPUS_PER_TASK, else 12)
      --config FILE            settings file, see config/circusnv.example.conf
      --skip-checks            do not check inputs, modules and tools before starting
      --dry-run                check everything and show what would be run, then stop
  -h, --help                  show this help and exit
      --version                show the version and exit

Relative paths are taken from the directory the pipeline is started from.
EOF
}

# Print the help text to stderr and stop with an error (wrong usage).
usage() {
  print_help >&2
  exit 1
}

die_usage() {
  echo "ERROR: $*" >&2
  echo "Run 'CircuSNV.sh --help' for the list of options." >&2
  exit 1
}

# Find --config before the other options, so the config file can set defaults that
# command-line options then override.
find_config_file() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --config) [[ $# -ge 2 ]] || die_usage "--config needs a file"; CONFIG_FILE=$2; shift 2;;
      --config=*) CONFIG_FILE=${1#--config=}; shift;;
      *) shift;;
    esac
  done
}

load_config_file() {
  [[ -n "$CONFIG_FILE" ]] || return 0
  [[ -f "$CONFIG_FILE" ]] || die_usage "config file not found: $CONFIG_FILE"
  # shellcheck source=/dev/null
  source "$CONFIG_FILE"
}

parse_args() {
  local opt val
  while [[ $# -gt 0 ]]; do
    opt=$1
    val=""
    # Accept '--option value', '--option=value' and, as getopts did, '-Ivalue'
    if [[ "$opt" == --*=* ]]; then
      val=${opt#*=}
      opt=${opt%%=*}
      set -- "$opt" "$val" "${@:2}"
    elif [[ "$opt" =~ ^-[IRLVDCPQMSHBXTUO].+ ]]; then
      val=${opt:2}
      opt=${opt:0:2}
      set -- "$opt" "$val" "${@:2}"
    fi
    case "$opt" in
      -h|--help) print_help; exit 0;;
      --version) echo "CircuSNV $(circusnv_version)"; exit 0;;
      --skip-checks) SKIP_CHECKS="true"; shift; continue;;
      --dry-run) DRY_RUN="true"; shift; continue;;
      -I|--input|-R|--reference|-L|--regions|-V|--min-vaf|-D|--min-depth|-C|--min-callers|\
      -P|--pon|-Q|--min-base-quality|-M|--min-alt-reads|-S|--strand-bias|-H|--mnv|\
      -B|--mnv-prior|-X|--crosslink-table|-T|--matched-control|-U|--unmet|-O|--output|\
      --threads|--config)
        [[ $# -ge 2 ]] || die_usage "option $opt needs a value"
        val=$2
        shift 2
        ;;
      -*) die_usage "unknown option: $opt";;
      *) die_usage "unexpected argument: $opt";;
    esac
    case "$opt" in
      -I|--input) Inputbam=$val;;
      -R|--reference) Reference=$val;;
      -L|--regions) Listregions=$val;;
      -V|--min-vaf) VAF=$val;;
      -D|--min-depth) RDP=$val;;
      -C|--min-callers) Calls=$val;;
      -P|--pon) PoN=$val;;
      -Q|--min-base-quality) Qual=$val;;
      -M|--min-alt-reads) MRD=$val;;
      -S|--strand-bias) SB=$val;;
      -H|--mnv) MNV=$val;;
      -B|--mnv-prior) PRIOR=$val;;
      -X|--crosslink-table) CrosslinkTable=$val;;
      -T|--matched-control) Matched=$val;;
      -U|--unmet) UNMET=$val;;
      -O|--output) OutputDir=$val;;
      --threads) THREADS=$val;;
      --config) ;;   # already handled by find_config_file
    esac
  done
}

# Remove trailing slashes from directory arguments, so paths are built the same way everywhere.
normalize_args() {
  [[ "$Inputbam" != "/" ]] && Inputbam=${Inputbam%/}
  [[ "$OutputDir" != "/" ]] && OutputDir=${OutputDir%/}
  [[ -n "$PoN" && "$PoN" != "/" ]] && PoN=${PoN%/}
  return 0
}

is_int() { [[ "$1" =~ ^[0-9]+$ ]]; }
is_number() { [[ "$1" =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; }
# True when $1 lies within [$2, $3] (decimal comparison with awk).
in_range() { awk -v x="$1" -v lo="$2" -v hi="$3" 'BEGIN { exit !(x + 0 >= lo + 0 && x + 0 <= hi + 0) }'; }

# Check all options and report every problem at once.
validate_args() {
  local errors=()

  [[ -n "$Inputbam" ]] || errors+=("-I/--input is required")
  [[ -n "$Reference" ]] || errors+=("-R/--reference is required")
  [[ -n "$Listregions" ]] || errors+=("-L/--regions is required")

  if [[ -n "$Inputbam" && ! -e "$Inputbam" ]]; then
    errors+=("input not found: $Inputbam")
  elif [[ -f "$Inputbam" && "$Inputbam" != *.bam ]]; then
    errors+=("input file is not a .bam file: $Inputbam")
  elif [[ -d "$Inputbam" ]] && ! compgen -G "$Inputbam/*.bam" >/dev/null; then
    errors+=("input directory contains no .bam files: $Inputbam")
  fi
  [[ -z "$Reference" || -f "$Reference" ]] || errors+=("reference FASTA not found: $Reference")
  [[ -z "$Listregions" || -f "$Listregions" ]] || errors+=("regions BED not found: $Listregions")
  [[ -z "$PoN" || -e "$PoN" ]] || errors+=("panel of normals not found: $PoN")
  [[ -z "$CrosslinkTable" || -f "$CrosslinkTable" ]] || errors+=("crosslink table not found: $CrosslinkTable")
  [[ ! -e "$OutputDir" || -d "$OutputDir" ]] || errors+=("output path exists and is not a directory: $OutputDir")

  { is_number "$VAF" && in_range "$VAF" 0 1; } || errors+=("-V/--min-vaf must be a number between 0 and 1 (got '$VAF')")
  is_int "$RDP" || errors+=("-D/--min-depth must be a whole number (got '$RDP')")
  { is_int "$Calls" && (( Calls >= 1 && Calls <= 4 )); } || errors+=("-C/--min-callers must be 1, 2, 3 or 4 (got '$Calls')")
  is_int "$Qual" || errors+=("-Q/--min-base-quality must be a whole number (got '$Qual')")
  is_int "$MRD" || errors+=("-M/--min-alt-reads must be a whole number (got '$MRD')")
  { is_number "$SB" && ! in_range "$SB" 0 0; } || errors+=("-S/--strand-bias must be a number above 0 (got '$SB')")
  { is_number "$PRIOR" && in_range "$PRIOR" 0 1; } || errors+=("-B/--mnv-prior must be a number between 0 and 1 (got '$PRIOR')")
  { is_int "$THREADS" && (( THREADS >= 1 )); } || errors+=("--threads must be a whole number of at least 1 (got '$THREADS')")

  [[ "$MNV" == "true" || "$MNV" == "false" ]] || errors+=("-H/--mnv must be 'true' or 'false' (got '$MNV')")
  [[ "$UNMET" == "true" || "$UNMET" == "false" ]] || errors+=("-U/--unmet must be 'true' or 'false' (got '$UNMET')")
  [[ "${Matched,,}" == "true" || "${Matched,,}" == "false" ]] || errors+=("-T/--matched-control must be 'true' or 'false' (got '$Matched')")
  if [[ "${Matched,,}" == "true" && -z "$CrosslinkTable" ]]; then
    errors+=("-T/--matched-control true needs -X/--crosslink-table with a matched_control column")
  fi

  if (( ${#errors[@]} )); then
    printf 'ERROR: %s\n' "${errors[@]}" >&2
    echo "Run 'CircuSNV.sh --help' for the list of options." >&2
    exit 1
  fi
}
