# shellcheck shell=bash
# Shared helpers: logging, module loading and running Python helpers in the env.
# Sourced by CircuSNV.sh; uses its global settings.

# Load the modules most steps need (htslib, bcftools).
load_modules() {
  module load "$MODULE_ZLIB"
  module load "$MODULE_HTSLIB"
  module load "$MODULE_BCFTOOLS"
}

# Run a Python helper inside the pipeline's virtual environment.
pythonscript() {
  source "$ENV_DIR/bin/activate"
  python3 "$@"
  deactivate
}

# Create the index of a BAM file if it has none (<name>.bam.bai or <name>.bai).
# Usage: ensure_bam_index <file.bam> <samtools module>
ensure_bam_index() {
  local bam=$1
  local samtools_module=$2
  local name index candidate
  name=$(basename "$bam" .bam)
  index=""
  for candidate in "${bam}.bai" "$(dirname "$bam")/${name}.bai"; do
    if [[ -e "$candidate" ]]; then
      index=$candidate
      break
    fi
  done

  if [[ -z "$index" ]]; then
    echo "BAM index file for ${name} not found. Creating BAM index..."
    (
      module load "$samtools_module"
      samtools index "$bam"
    )
  else
    echo "BAM index file for ${name} already exists: ${index}"
    # Make the index newer than the BAM, so tools do not warn about an outdated index
    touch "$index"
  fi
}

# Print a step heading. The run log puts a timestamp in front of every output line.
log() {
  echo "== $* =="
}

# Directory for run logs and run state inside the output directory. It is hidden so it does
# not look like a sample directory.
run_state_dir() {
  echo "$OutputDir/.circusnv"
}

# Send all further output to the terminal (or Slurm .out/.err) and to a log file per run,
# with a timestamp in front of every line.
start_run_log() {
  local state_dir
  state_dir=$(run_state_dir)
  mkdir -p "$state_dir"
  RUN_LOG="$state_dir/run_$(date '+%Y%m%d_%H%M%S')_$$.log"
  exec 3>&1 4>&2
  exec > >(perl -MPOSIX -pe 'BEGIN { $| = 1 } print strftime("[%Y-%m-%d %H:%M:%S] ", localtime)' | tee -a "$RUN_LOG")
  LOG_PID_OUT=$!
  exec 2> >(perl -MPOSIX -pe 'BEGIN { $| = 1 } print strftime("[%Y-%m-%d %H:%M:%S] ", localtime)' | tee -a "$RUN_LOG" >&4)
  LOG_PID_ERR=$!
  log "Log file: $RUN_LOG"
}

# Restore the original output and let the log writers finish.
finish_run_log() {
  [[ -n "${RUN_LOG:-}" ]] || return 0
  exec 1>&3 2>&4
  wait "${LOG_PID_OUT:-}" "${LOG_PID_ERR:-}" 2>/dev/null || true
  RUN_LOG=""
}
