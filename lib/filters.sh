# shellcheck shell=bash
# UNMET score annotation and filtering.
# Sourced by CircuSNV.sh; uses its global settings.

# Function to add UNMET score and filter 
unmet_filter() {
  local samp="$1"                             # sample name, e.g. “SAMPLE123”
  local curr_file="$2"                        # merged_vcfs.txt for that sample
  local bed_file="$3"                         # path to the .bed you want
  local threshold="$4"                        # numeric cutoff, e.g. 0.93
  local label="$5"                            # short label for logging
  local tmp_bed_tdf="${bed_file%.bed}.tdf"
  local out_with="${curr_file%.txt}_UNMET.txt"
  

  # 1) Make sure the UNMET.bed exists
  if [[ ! -f "$bed_file" ]]; then
    echo "UNMET is true and $bed_file not found — downloading + converting…"
    wget -q -O "$tmp_bed_tdf" "http://web16.kazusa.or.jp/data/hg38/$(basename "$tmp_bed_tdf")"
    (
      module load IGV
      igvtools tdftobedgraph "$tmp_bed_tdf" "$bed_file" 2>/dev/null
    )
    echo "UNMET: created $bed_file"
  fi

  # 2) Attach UNMET scores
  echo "[$samp] $label: variants before = $(wc -l < "$curr_file")"
  (
    module load BEDTools
    bedtools intersect -a <(
      awk 'BEGIN{FS=OFS="\t"} { print $1, $2, $2+1, $0 }' "$curr_file"
    ) -b "$bed_file" -loj \
    | awk -v OFS="\t" -v th="$threshold" '
        {
          sc = ($12 == "." ? "NA" : $12)
          printf "%s\t%s\t%s\t%s\t%s_UNMET:%s\n",
                 $4, $5, $6, $7, $8, sc
        }' > "$out_with"
  )
  # 3) Filter by threshold
  awk -F'\t' -v th="$threshold" '
    {
      split($5, a, "_UNMET:")
      sc = a[2]
      if (sc == "NA" || sc+0 < th) print
    }
  ' "$out_with" > "$curr_file"

  echo "[$samp] $label: variants after  = $(wc -l < "$curr_file")"
}
