#!/usr/bin/env bash
# Count the records of a 'bcftools isec' sites.txt that lie outside a BED file, per tool
# combination and per distance to the nearest region. Prints counts only (no positions or
# sample names), so the output can be shared.
#
# Usage: off_target_check.sh <sites.txt> <regions.bed>
#
# sites.txt has the columns CHROM POS REF ALT and a tool flag string in the order SiNVICT,
# Mutect2, LoFreq, VarDict, as CircuSNV writes it (<sample>/mnvcalling/sites.txt or
# <sample>/isec/sites.txt). A record is inside when its REF span overlaps a region: the same
# rule the pipeline uses (bcftools --targets-overlap 1).
#
# Exit status: 0 = nothing outside the BED, 1 = records outside, 2 = usage error.
set -euo pipefail

if [[ $# -ne 2 || ! -f "${1:-}" || ! -f "${2:-}" ]]; then
  echo "Usage: $(basename "$0") <sites.txt> <regions.bed>" >&2
  exit 2
fi

awk -F '\t' '
  BEGIN { split("SiNVICT Mutect2 LoFreq VarDict", tool, " ") }
  # First file: the regions (BED: 0-based, end exclusive); header and comment lines skipped
  NR == FNR {
    sub(/\r$/, "")
    if ($0 ~ /^(browser|track|#)/ || NF < 3) next
    n[$1]++
    rs[$1, n[$1]] = $2
    re[$1, n[$1]] = $3
    next
  }
  # Second file: the sites
  {
    total++
    beg = $2 - 1
    end = beg + length($3)
    best = -1
    for (i = 1; i <= n[$1]; i++) {
      if (rs[$1, i] < end && re[$1, i] > beg) { best = 0; break }
      d = (re[$1, i] <= beg) ? beg - re[$1, i] + 1 : rs[$1, i] - end + 1
      if (best < 0 || d < best) best = d
    }
    if (best == 0) next
    outside++

    combo = ""
    for (i = 1; i <= 4; i++) if (substr($5, i, 1) == "1") combo = combo (combo == "" ? "" : "+") tool[i]
    if (combo == "") combo = "(no tool flag)"
    by_combo[combo]++

    if (best < 0) cls = 6
    else if (best <= 5) cls = 1
    else if (best <= 50) cls = 2
    else if (best <= 150) cls = 3
    else if (best <= 300) cls = 4
    else cls = 5
    by_dist[cls]++
  }
  END {
    split("1-5 bp|6-50 bp|51-150 bp|151-300 bp|more than 300 bp|contig without regions", label, "|")
    printf "Records: %d\nOutside the BED: %d\n", total, outside
    if (outside > 0) {
      print "\nBy tool combination:"
      cmd = "sort -k1,1nr -k2"
      for (c in by_combo) printf "  %7d  %s\n", by_combo[c], c | cmd
      close(cmd)
      print "\nBy distance to the nearest region:"
      for (i = 1; i <= 6; i++) if (by_dist[i] > 0) printf "  %7d  %s\n", by_dist[i], label[i]
    }
    exit (outside > 0)
  }
' "$2" "$1"
