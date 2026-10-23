#!/usr/bin/env bash
# preflight.sh - read-only checks before committing and pushing a change.
#
# Usage (from anywhere inside the repository):
#   bash .claude/skills/github-workflow/scripts/preflight.sh
#
# What it looks at: files changed in the working tree (staged, unstaged, untracked
# and not ignored) plus files committed on this branch but not yet on main.
# It changes nothing.
#
# Output lines: PASS / FAIL / WARN / NOTE.
# Exit status: 0 = no FAIL, 1 = at least one FAIL, 2 = could not run.
#
# Environment:
#   PREFLIGHT_ALLOW_VCF=1   do not fail on VCF files (small synthetic test data only)

set -uo pipefail

FAILS=0
WARNS=0
NOTES=0

ok()   { printf 'PASS  %s\n' "$*"; }
fail() { printf 'FAIL  %s\n' "$*"; FAILS=$((FAILS + 1)); }
warn() { printf 'WARN  %s\n' "$*"; WARNS=$((WARNS + 1)); }
note() { printf 'NOTE  %s\n' "$*"; NOTES=$((NOTES + 1)); }

# Run one check function and print PASS only if it produced no FAIL and no WARN.
section() {
  local label="$1" fn="$2" f_before="$FAILS" w_before="$WARNS"
  "$fn"
  if [ "$FAILS" -eq "$f_before" ] && [ "$WARNS" -eq "$w_before" ]; then
    ok "$label"
  fi
}

# ---------------------------------------------------------------- repository
if ! root="$(git rev-parse --show-toplevel 2>/dev/null)"; then
  echo "FAIL  not inside a git repository"
  exit 2
fi
cd "$root" || exit 2

# ------------------------------------------------------------------- branch
branch="$(git branch --show-current 2>/dev/null || true)"
if [ -z "$branch" ]; then
  fail "detached HEAD: switch to a branch first"
elif [ "$branch" = "main" ] || [ "$branch" = "master" ]; then
  fail "currently on '$branch': create a branch first (git switch -c <type>/<description>)"
else
  ok "on branch '$branch'"
fi

# ------------------------------------------------------------ list of files
base=""
for candidate in origin/main main; do
  if git rev-parse --verify -q "$candidate" >/dev/null 2>&1; then
    base="$candidate"
    break
  fi
done

files=()
while IFS= read -r -d '' f; do
  [ -f "$f" ] && files+=("$f")
done < <(
  {
    git diff --name-only -z --diff-filter=ACMR HEAD 2>/dev/null
    git ls-files -z --others --exclude-standard 2>/dev/null
    if [ -n "$base" ]; then
      git diff --name-only -z --diff-filter=ACMR "$base"...HEAD 2>/dev/null
    fi
  } | sort -zu
)

if [ "${#files[@]}" -eq 0 ]; then
  note "no changed files found: nothing to check"
  printf '\nResult: %d FAIL, %d WARN, %d NOTE\n' "$FAILS" "$WARNS" "$NOTES"
  [ "$FAILS" -eq 0 ] && exit 0 || exit 1
fi

printf 'Checking %d changed file(s):\n' "${#files[@]}"
shown=0
for f in "${files[@]}"; do
  printf '      %s\n' "$f"
  shown=$((shown + 1))
  if [ "$shown" -ge 30 ]; then
    printf '      ... and %d more\n' "$((${#files[@]} - shown))"
    break
  fi
done
echo

# --------------------------------------------------------- tools that exist
PY=""
for candidate in python3 python py; do
  if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys' >/dev/null 2>&1; then
    PY="$candidate"
    break
  fi
done
HAVE_SHELLCHECK=0
command -v shellcheck >/dev/null 2>&1 && HAVE_SHELLCHECK=1

# ------------------------------------------------------------------- checks
check_data_files() {
  local f lc
  for f in "${files[@]}"; do
    lc="${f,,}"
    case "$lc" in
      *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.bam|*.bai|*.cram|*.crai|*.sam)
        fail "sequencing data must never be committed: $f" ;;
      *.vcf|*.vcf.gz|*.vcf.gz.tbi|*.bcf|*.tbi)
        if [ "${PREFLIGHT_ALLOW_VCF:-0}" = "1" ]; then
          warn "VCF file allowed by PREFLIGHT_ALLOW_VCF: $f (must be synthetic and small)"
        else
          fail "VCF file in a public repository: $f (only small synthetic test data, agreed with Nick)"
        fi ;;
    esac
    case "$lc" in
      output/*|results/*|tmp/*|*/output/*|*/results/*|*/tmp/*)
        fail "pipeline output directory: $f" ;;
      *.log)
        warn "log file: $f" ;;
    esac
  done
}

check_sensitive_names() {
  local f name
  for f in "${files[@]}"; do
    name="${f,,}"
    name="${name##*/}"
    case "$name" in
      settings.local.json|claude.local.md|.env|.env.*|*.pem|*.pfx|*.p12|*.key|\
      id_rsa|id_rsa.*|id_ed25519|id_ed25519.*|credentials|credentials.*)
        fail "file that should stay private: $f" ;;
    esac
  done
}

SECRET_PATTERNS=(
  '[-]{5}BEGIN [A-Z ]*PRIVATE KEY[-]{5}'
  'ghp_[A-Za-z0-9]{36}'
  'gho_[A-Za-z0-9]{36}'
  'github_pat_[A-Za-z0-9_]{22,}'
  'AKIA[0-9A-Z]{16}'
  'xox[baprs]-[A-Za-z0-9-]{10,}'
  'sk-ant-[A-Za-z0-9_-]{20,}'
)
ASSIGNMENT_PATTERN="(password|passwd|secret|api[_-]?key|token)[[:space:]]*[=:][[:space:]]*['\"][^'\"[:space:]]{6,}['\"]"

check_secrets() {
  local f p lines
  for f in "${files[@]}"; do
    for p in "${SECRET_PATTERNS[@]}"; do
      if grep -Iqe "$p" -E -- "$f" 2>/dev/null; then
        lines="$(grep -InE -e "$p" -- "$f" 2>/dev/null | cut -d: -f1 | head -3 | paste -sd, -)"
        fail "looks like a secret or token in $f (line $lines): remove it and rotate it if it was ever pushed"
      fi
    done
    if grep -Iqie "$ASSIGNMENT_PATTERN" -E -- "$f" 2>/dev/null; then
      lines="$(grep -InEi -e "$ASSIGNMENT_PATTERN" -- "$f" 2>/dev/null | cut -d: -f1 | head -3 | paste -sd, -)"
      warn "possible hard-coded credential in $f (line $lines): check by hand"
    fi
  done
}

check_size() {
  local f size
  for f in "${files[@]}"; do
    size="$(wc -c <"$f" | tr -d ' ')"
    if [ "$size" -gt 10485760 ]; then
      fail "file larger than 10 MB: $f ($((size / 1048576)) MB): data does not belong in the repository"
    elif [ "$size" -gt 1048576 ]; then
      warn "file larger than 1 MB: $f ($((size / 1048576)) MB)"
    fi
  done
}

check_line_endings() {
  local f lc
  for f in "${files[@]}"; do
    lc="${f,,}"
    if grep -Iq $'\r' -- "$f" 2>/dev/null; then
      case "$lc" in
        *.sh|*.py)
          fail "CRLF line endings in $f: a script with CRLF fails on Linux (convert to LF)" ;;
        *.md|*.json|*.txt|*.yml|*.yaml|*.toml|*.cfg|.gitattributes|.gitignore)
          warn "CRLF line endings in $f: the repository uses LF" ;;
      esac
    fi
  done
}

SAW_SH=0

check_bash() {
  local f out first
  for f in "${files[@]}"; do
    case "${f,,}" in
      *.sh)
        SAW_SH=1
        if ! out="$(bash -n "$f" 2>&1)"; then
          fail "bash syntax error in $f: $(printf '%s' "$out" | head -1)"
        fi
        first="$(head -n 1 "$f" 2>/dev/null || true)"
        case "$first" in
          '#!'*) ;;
          *) warn "no shebang on the first line of $f" ;;
        esac
        if [ "$HAVE_SHELLCHECK" -eq 1 ]; then
          out="$(shellcheck -S warning -f gcc "$f" 2>&1 || true)"
          if [ -n "$out" ]; then
            warn "shellcheck found $(printf '%s\n' "$out" | wc -l | tr -d ' ') issue(s) in $f (first: $(printf '%s' "$out" | head -1))"
          fi
        fi ;;
    esac
  done
}

check_python() {
  local f out any=0
  for f in "${files[@]}"; do
    case "${f,,}" in
      *.py)
        any=1
        [ -z "$PY" ] && continue
        if ! out="$("$PY" -I -c 'import ast,sys; ast.parse(open(sys.argv[1],"rb").read(), sys.argv[1])' "$f" 2>&1)"; then
          fail "Python syntax error in $f: $(printf '%s' "$out" | tail -1)"
        fi ;;
    esac
  done
  if [ "$any" -eq 1 ] && [ -z "$PY" ]; then
    note "no working Python found: Python syntax not checked"
  fi
}

check_json() {
  local f out any=0
  for f in "${files[@]}"; do
    case "${f,,}" in
      *.json)
        any=1
        [ -z "$PY" ] && continue
        if ! out="$("$PY" -I -c 'import json,sys; json.load(open(sys.argv[1],encoding="utf-8-sig"))' "$f" 2>&1)"; then
          fail "invalid JSON in $f: $(printf '%s' "$out" | tail -1)"
        fi ;;
    esac
  done
  if [ "$any" -eq 1 ] && [ -z "$PY" ]; then
    note "no working Python found: JSON not checked"
  fi
}

section "no data files or output directories"  check_data_files
section "no private files"                      check_sensitive_names
section "no secrets or tokens"                  check_secrets
section "file sizes"                            check_size
section "line endings"                          check_line_endings
section "bash syntax"                           check_bash
section "python syntax"                         check_python
section "json validity"                         check_json

if [ "$SAW_SH" -eq 1 ] && [ "$HAVE_SHELLCHECK" -eq 0 ]; then
  note "shellcheck is not installed: bash scripts were only checked with 'bash -n'"
fi

# ------------------------------------------------------ pipeline code changed
pipeline_changed=0
for f in "${files[@]}"; do
  case "$f" in
    CircuSNV.sh|*.py)
      case "$f" in */*) ;; *) pipeline_changed=1 ;; esac ;;
  esac
done
if [ "$pipeline_changed" -eq 1 ]; then
  note "pipeline code changed: output may change. A Nibbler test and a baseline comparison are needed before merging. Keep refactors and behaviour changes in separate commits."
fi

printf '\nResult: %d FAIL, %d WARN, %d NOTE\n' "$FAILS" "$WARNS" "$NOTES"
if [ "$FAILS" -gt 0 ]; then
  exit 1
fi
exit 0
