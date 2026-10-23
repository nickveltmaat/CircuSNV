# shellcheck shell=bash
# Shared helpers: module loading and running Python helpers in the env.
# Sourced by CircuSNV.sh; uses its global settings.

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
