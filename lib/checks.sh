# shellcheck shell=bash
# Pre-run checks: find missing inputs, resources, modules and tools before hours of compute
# are spent. Skipped with --skip-checks.
# Sourced by CircuSNV.sh; uses its global settings.

# Print one line per problem with modules and commands. Every group of modules is loaded in
# its own subshell, so modules from different toolchains do not clash.
tool_problems() {
  # Runs in a subshell (process substitution): the module system is not written for set -u
  set +u
  if ! type module >/dev/null 2>&1; then
    echo "the 'module' command is not available"
    return 0
  fi
  (
    local m t
    for m in "$MODULE_ZLIB" "$MODULE_HTSLIB" "$MODULE_BCFTOOLS" "$MODULE_SAMTOOLS" "$MODULE_R" "$MODULE_GATK"; do
      module load "$m" >/dev/null 2>&1 || echo "cannot load module $m"
    done
    for t in bcftools bgzip samtools gatk vardict-java teststrandbias.R var2vcf_valid.pl \
             lofreq bam-readcount sinvict Rscript perl bc awk; do
      command -v "$t" >/dev/null 2>&1 || echo "command not found: $t"
    done
  )
  if [[ -n "$PoN" && -d "$PoN" ]]; then
    ( module load "$MODULE_SAMTOOLS_PON" >/dev/null 2>&1 || echo "cannot load module $MODULE_SAMTOOLS_PON" )
  fi
  if [[ "$UNMET" == "true" ]]; then
    (
      module load "$MODULE_BEDTOOLS" >/dev/null 2>&1 || echo "cannot load module $MODULE_BEDTOOLS"
      command -v bedtools >/dev/null 2>&1 || echo "command not found: bedtools"
    )
    if [[ ! -f "$UNMET_BED" ]]; then
      (
        command -v wget >/dev/null 2>&1 || echo "command not found: wget (needed to download the UNMET scores)"
        module load "$MODULE_IGV" >/dev/null 2>&1 || echo "cannot load module $MODULE_IGV"
        command -v igvtools >/dev/null 2>&1 || echo "command not found: igvtools (needed to convert the UNMET scores)"
      )
    fi
  fi
  if [[ "$MNV" == "true" ]]; then
    ( module load "$MODULE_GCC_MNVISTA" >/dev/null 2>&1 || echo "cannot load module $MODULE_GCC_MNVISTA" )
  fi
}

# Print one line per problem with the Python environment.
env_problems() {
  # Runs in a subshell (process substitution): the venv activate script is not written for set -u
  set +u
  if [[ ! -f "$ENV_DIR/bin/activate" ]]; then
    echo "Python env not found: $ENV_DIR (set ENV_DIR in a --config file)"
    return 0
  fi
  (
    source "$ENV_DIR/bin/activate"
    command -v python3 >/dev/null 2>&1 || echo "python3 not found in $ENV_DIR"
    command -v oc >/dev/null 2>&1 || echo "OpenCRAVAT ('oc') not found in $ENV_DIR"
    python3 -c 'import pandas, numpy, openpyxl' >/dev/null 2>&1 || echo "Python packages pandas, numpy and openpyxl are not all installed in $ENV_DIR"
  )
}

# Print one line per sample listed in the crosslink table whose BAM is missing.
crosslink_problems() {
  [[ -n "$CrosslinkTable" && -d "$Inputbam" ]] || return 0
  local s
  while IFS= read -r s; do
    [[ -f "$Inputbam/$s.bam" ]] || echo "sample '$s' from the crosslink table has no BAM: $Inputbam/$s.bam"
  done < <(awk -F $'\t' 'NR>1 && ($2=="baseline" || $2=="follow-up") && !seen[$1]++ {print $1}' "$CrosslinkTable")
}

# Print one line per matched-control problem (-T true). The control BAMs are looked up in the
# input directory, as before; a missing control makes that sample fail at the read counts.
matched_control_problems() {
  [[ "${Matched,,}" == "true" ]] || return 0
  local s
  if [[ -f "$Inputbam" ]]; then
    echo "-T true needs -I to be a directory: matched control BAMs are looked up in the input directory (see docs/known-issues.md)"
    return 0
  fi
  while IFS= read -r s; do
    [[ -f "$Inputbam/$s.bam" ]] || echo "matched control BAM not found: $Inputbam/$s.bam"
  done < <(awk -F $'\t' 'NR>1 {gsub(/\r/,"",$4); if ($4 != "" && !seen[$4]++) print $4}' "$CrosslinkTable")
}

# Print one line per warning (problems the pipeline may survive).
check_warnings() {
  [[ -f "$Reference.fai" ]] || echo "reference index not found: $Reference.fai"
  [[ -f "${Reference%.*}.dict" ]] || echo "sequence dictionary for GATK not found: ${Reference%.*}.dict"
  if [[ "$UNMET" == "true" && ! -f "$UNMET_BED" ]]; then
    echo "UNMET scores not found: $UNMET_BED; they will be downloaded and converted"
  fi
  return 0
}

run_checks() {
  local problems=() warnings=() line out_parent

  if [[ "$MNV" == "true" && ! -x "$MNVISTA" ]]; then
    problems+=("MNVista not found or not executable: $MNVISTA (set MNVISTA in a --config file)")
  fi
  # MNVista does not fail on a missing -C file: it silently calls MNVs without the panel
  if [[ "$MNV" == "true" && -n "$PoN" && ! -f "$MNV_PON" ]]; then
    problems+=("MNV panel of normals not found: $MNV_PON (needed for MNVista with -P; set it with --mnv-pon)")
  fi
  if [[ -e "$OutputDir" ]]; then
    [[ -w "$OutputDir" ]] || problems+=("output directory is not writable: $OutputDir")
  else
    out_parent=$(dirname "$OutputDir")
    [[ -w "$out_parent" ]] || problems+=("cannot create the output directory in: $out_parent")
  fi
  if [[ -n "$PoN" ]]; then
    if [[ -e "$PON_DIR" ]]; then
      [[ -w "$PON_DIR" ]] || problems+=("PoN work directory is not writable: $PON_DIR")
    else
      [[ -w "$(dirname "$PON_DIR")" ]] || problems+=("cannot create the PoN work directory: $PON_DIR")
    fi
  fi

  while IFS= read -r line; do [[ -n "$line" ]] && problems+=("$line"); done < <(env_problems)
  while IFS= read -r line; do [[ -n "$line" ]] && problems+=("$line"); done < <(tool_problems)
  while IFS= read -r line; do [[ -n "$line" ]] && problems+=("$line"); done < <(crosslink_problems)
  while IFS= read -r line; do [[ -n "$line" ]] && problems+=("$line"); done < <(matched_control_problems)
  while IFS= read -r line; do [[ -n "$line" ]] && warnings+=("$line"); done < <(check_warnings)

  if (( ${#warnings[@]} )); then
    printf 'WARNING: %s\n' "${warnings[@]}" >&2
  fi
  if (( ${#problems[@]} )); then
    printf 'ERROR: %s\n' "${problems[@]}" >&2
    echo "Pre-run checks failed. Fix the problems above, or run with --skip-checks if a check is wrong." >&2
    exit 1
  fi
  echo "Pre-run checks passed."
}
