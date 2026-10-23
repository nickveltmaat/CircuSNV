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

# Initiating settings #
# Locate the installation directory: helper scripts, the Python env, MNVista and resources live
# there. Under sbatch this script runs from a spool copy, so fall back to the submitted path.
locate_home() {
  local dir candidate
  dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  if [[ -f "$dir/scripts/isec_to_vcf.py" ]]; then
    echo "$dir"
    return 0
  fi
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    candidate=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^[[:space:]]*Command=\([^[:space:]]*\).*/\1/p' | head -n 1 || true)
    if [[ -n "$candidate" && "$candidate" != /* && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
      candidate="$SLURM_SUBMIT_DIR/$candidate"
    fi
    if [[ -n "$candidate" && -f "$(dirname "$candidate")/scripts/isec_to_vcf.py" ]]; then
      (cd "$(dirname "$candidate")" && pwd)
      return 0
    fi
  fi
  if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/scripts/isec_to_vcf.py" ]]; then
    echo "$SLURM_SUBMIT_DIR"
    return 0
  fi
  return 1
}

if [[ -z "${CIRCUSNV_HOME:-}" ]] && ! CIRCUSNV_HOME=$(locate_home); then
  echo "ERROR: cannot find the CircuSNV installation directory; set CIRCUSNV_HOME" >&2
  exit 1
fi

for module in config args common callers filters annotate mnv pon samples; do
  source "$CIRCUSNV_HOME/lib/$module.sh"
done

# Settings: built-in defaults, then an optional config file, then the command line.
apply_default_config
set_default_options
find_config_file "$@"
load_config_file
parse_args "$@"
normalize_args
validate_args

# Upper bound of the strand-bias window (two decimals, as before)
SBU=$(echo "scale=2; 1 / $SB" | bc)

echo -e "\nCircuSNV $(circusnv_version)\nInstallation: $CIRCUSNV_HOME\nThreads: $THREADS"
echo -e "\nSettings: \nInput bam(s): $Inputbam\nReference: $Reference\nPanel: $Listregions\nMinimum VAF: $VAF\nMinimum Read Depth: $RDP\nMinimum Base Quality: $Qual\nMinimum Overlapping Calls: $Calls\nPoN: $PoN\nStrand Bias: $SB\nMultiNucleotide Variant Calling: $MNV\nBayesian Prior for MNV calling: $PRIOR\nCrosslink Table: $CrosslinkTable\nMatched Control Logic: $Matched\nUNMET & repeat filtering: $UNMET\nMain output directory: $OutputDir"

prepare_pon
run_samples

echo 'Finished run!'
