# shellcheck shell=bash
# Shared helpers: module loading and running Python helpers in the env.
# Sourced by CircuSNV.sh; uses its global settings.

# Function to load modules
load_modules() {
  module load $MODULE_ZLIB
  module load $MODULE_HTSLIB
  module load $MODULE_BCFTOOLS
}

# Function to run Python script
pythonscript() {
  source "$ENV_DIR/bin/activate"
  python3 "$@"
  deactivate
}
