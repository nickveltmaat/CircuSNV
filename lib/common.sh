# shellcheck shell=bash
# Shared helpers: logging, module loading and running Python helpers in the env.
# Sourced by CircuSNV.sh; uses its global settings.

# Load one environment module. The module system and the venv activate script are not
# written for 'set -u', so unset-variable checks are paused while they run.
load_module() {
  local rc=0
  set +u
  module load "$1" || rc=$?
  set -u
  if (( rc != 0 )); then
    echo "ERROR: could not load module $1" >&2
  fi
  return "$rc"
}

# Load the modules most steps need (htslib, bcftools).
load_modules() {
  load_module "$MODULE_ZLIB"
  load_module "$MODULE_HTSLIB"
  load_module "$MODULE_BCFTOOLS"
}

activate_env() {
  set +u
  # shellcheck source=/dev/null
  source "$ENV_DIR/bin/activate"
  set -u
}

deactivate_env() {
  set +u
  deactivate
  set -u
}

# Run a Python helper inside the pipeline's virtual environment.
pythonscript() {
  activate_env
  python3 "$@"
  deactivate_env
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
      load_module "$samtools_module"
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
  if command -v perl >/dev/null 2>&1; then
    exec > >(perl -MPOSIX -pe 'BEGIN { $| = 1 } print strftime("[%Y-%m-%d %H:%M:%S] ", localtime)' | tee -a "$RUN_LOG")
    LOG_PID_OUT=$!
    exec 2> >(perl -MPOSIX -pe 'BEGIN { $| = 1 } print strftime("[%Y-%m-%d %H:%M:%S] ", localtime)' | tee -a "$RUN_LOG" >&4)
    LOG_PID_ERR=$!
  else
    # Without perl: same log file, no timestamps
    exec > >(tee -a "$RUN_LOG")
    LOG_PID_OUT=$!
    exec 2> >(tee -a "$RUN_LOG" >&4)
    LOG_PID_ERR=$!
  fi
  log "Log file: $RUN_LOG"
}

# Restore the original output and let the log writers finish.
finish_run_log() {
  [[ -n "${RUN_LOG:-}" ]] || return 0
  exec 1>&3 2>&4
  wait "${LOG_PID_OUT:-}" "${LOG_PID_ERR:-}" 2>/dev/null || true
  RUN_LOG=""
}
