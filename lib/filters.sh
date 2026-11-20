# shellcheck shell=bash
# UNMET score annotation and filtering, and the PoN filter for indels by position.
# Sourced by CircuSNV.sh; uses its global settings.

# Remove the indels of a merged variant list (chrom pos ref alt info) at a position where the
# blacklist has an indel, whatever its alleles. The PoN filter itself (bcftools isec) only
# removes exact allele matches, while polymerase slippage in a short repeat can show up as a
# deletion in the normals and as an insertion in a tumor sample, both left-aligned to the same
# anchor base. SNVs are not touched. The removed records are kept next to the list, in
# <list>_PoN_indel_position_removed.txt, for review.
# Usage: pon_indel_position_filter <samplename> <merged_vcfs.txt> <blacklist.vcf.gz>
pon_indel_position_filter() {
  local samp=$1 curr_file=$2 blacklist=$3
  local removed="${curr_file%.txt}_PoN_indel_position_removed.txt"
  local positions
  positions=$(mktemp "${curr_file%.txt}.pon_indel_positions.XXXXXX")

  # chrom<TAB>pos of every indel in the blacklist
  (
    load_modules
    bcftools view -H "$blacklist" | awk -F '\t' 'length($4) != length($5) { print $1 "\t" $2 }' > "$positions"
  )

  echo "[$samp] PoN indel position filtering: variants before = $(wc -l < "$curr_file") ($(wc -l < "$positions") indel position(s) in the PoN)"
  : > "$removed"
  # FILENAME, not NR==FNR: the positions file may be empty
  awk -F '\t' -v OFS='\t' -v posfile="$positions" -v removed="$removed" '
    FILENAME == posfile { pon[$1 "\t" $2] = 1; next }
    length($3) != length($4) && (($1 "\t" $2) in pon) { print > removed; next }
    { print }
  ' "$positions" "$curr_file" > "$curr_file.tmp"
  mv "$curr_file.tmp" "$curr_file"
  rm -f "$positions"

  echo "[$samp] PoN indel position filtering: variants after  = $(wc -l < "$curr_file"); removed indels listed in $(basename "$removed")"
}

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
    wget -q -O "$tmp_bed_tdf" "$UNMET_URL_BASE/$(basename "$tmp_bed_tdf")"
    (
      load_module "$MODULE_IGV"
      igvtools tdftobedgraph "$tmp_bed_tdf" "$bed_file" 2>/dev/null
    )
    echo "UNMET: created $bed_file"
  fi

  # 2) Attach UNMET scores
  echo "[$samp] $label: variants before = $(wc -l < "$curr_file")"
  (
    load_module "$MODULE_BEDTOOLS"
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
