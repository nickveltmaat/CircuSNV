# shellcheck shell=bash
# Region (BED) files: the cleaned panel BED, the optional MNV flanking, and restricting VCFs to
# a BED. Sourced by CircuSNV.sh; uses its global settings.
#
# prepare_regions sets, once per run and before the panel of normals and the samples:
#   REGIONS_PANEL    the -L BED without header lines; the main path (and the normals of the
#                    PoN) is restricted to it
#   REGIONS_MNV      REGIONS_PANEL with the marked regions flanked (equal to REGIONS_PANEL
#                    without flanking); the MNV input is restricted to it
#   REGIONS_CALLING  the BED the tumor callers run on: the -L file as given, or REGIONS_MNV
#                    with flanking. The normals of the PoN always run on the -L file.

# True when MNV flanking is in effect: a pattern is set and MNV calling is on.
mnv_flank_active() {
  [[ "$MNV" == "true" && -n "$MNV_FLANK_REGIONS" ]]
}

# Print the data lines of a BED file: no browser, track, comment or empty lines, no CR.
bed_data_lines() {
  awk '{ sub(/\r$/, "") } !/^(browser|track|#)/ && NF' "$1"
}

# Print the smallest number of tab-separated columns of the data lines of a BED file.
bed_min_columns() {
  bed_data_lines "$1" | awk -F '\t' 'NR == 1 || NF < m { m = NF } END { print m + 0 }'
}

# Print "<marked> <total>": the number of regions whose last column matches MNV_FLANK_REGIONS
# (an extended regular expression), and the number of regions.
count_flank_regions() {
  bed_data_lines "$Listregions" \
    | FLANK_PATTERN="$MNV_FLANK_REGIONS" awk -F '\t' '$NF ~ ENVIRON["FLANK_PATTERN"] { m++ } END { print m + 0, NR }'
}

# Describe the MNV flanking in one line (run log and --dry-run); no marked region is a warning.
describe_mnv_flank() {
  local marked total
  if ! mnv_flank_active; then
    echo "MNV flanking: off"
    return 0
  fi
  read -r marked total < <(count_flank_regions)
  echo "MNV flanking: $marked of $total region(s) match '$MNV_FLANK_REGIONS' (last BED column) and get $MNV_FLANK_BP bp on each side for calling and the MNV input"
  if (( marked == 0 )); then
    echo "WARNING: no region matches --mnv-flank-regions '$MNV_FLANK_REGIONS'; nothing is flanked" >&2
  fi
}

# Write the region files of this run to <output>/.circusnv/regions/ and set REGIONS_PANEL,
# REGIONS_MNV and REGIONS_CALLING.
prepare_regions() {
  local dir
  dir="$(run_state_dir)/regions"
  mkdir -p "$dir"

  REGIONS_PANEL="$dir/panel.bed"
  bed_data_lines "$Listregions" > "$REGIONS_PANEL"
  echo "Regions: $(wc -l < "$REGIONS_PANEL") panel region(s); per-tool calls outside them are removed"
  describe_mnv_flank

  if ! mnv_flank_active; then
    REGIONS_MNV=$REGIONS_PANEL
    REGIONS_CALLING=$Listregions
    return 0
  fi

  # Extend only the marked regions, then merge them with all other regions. The name column
  # is kept (VarDict reads it as the gene name, -g 4); merged regions get their names joined.
  : > "$dir/marked.bed"
  : > "$dir/other.bed"
  FLANK_PATTERN="$MNV_FLANK_REGIONS" awk -F '\t' -v marked="$dir/marked.bed" -v other="$dir/other.bed" \
    '{ if ($NF ~ ENVIRON["FLANK_PATTERN"]) print > marked; else print > other }' "$REGIONS_PANEL"
  (
    load_module "$MODULE_BEDTOOLS"
    { bedtools slop -i "$dir/marked.bed" -g "$Reference.fai" -b "$MNV_FLANK_BP"; cat "$dir/other.bed"; } \
      | bedtools sort -faidx "$Reference.fai" -i stdin \
      | bedtools merge -i stdin -c 4 -o distinct > "$dir/mnv.bed"
  )
  rm -f "$dir/marked.bed" "$dir/other.bed"
  REGIONS_MNV="$dir/mnv.bed"
  REGIONS_CALLING=$REGIONS_MNV
  echo "Regions: $(wc -l < "$REGIONS_MNV") region(s) after flanking and merging, used for calling and the MNV input"
}

# Keep only the records of a VCF whose REF span overlaps a region of a BED file, and log how
# many were removed. -T streams the VCF (no index needed) and gives no duplicate records for
# overlapping regions, as -R can; bcftools reads the file as 0-based BED because of its .bed
# suffix. The overlap rule is set explicitly: the default for -T is 0 (POS inside a region),
# 1 (record) also keeps a deletion whose padding base lies just before a region.
# Usage: restrict_to_regions <input.vcf[.gz]> <output.vcf[.gz]> <regions.bed> <label>
restrict_to_regions() {
  local input=$1 output=$2 bed=$3 label=$4
  local out_type=v before after
  [[ "$output" == *.gz ]] && out_type=z
  bcftools view -T "$bed" --targets-overlap 1 -O "$out_type" -o "$output" "$input"
  before=$(bcftools view -H "$input" | wc -l)
  after=$(bcftools view -H "$output" | wc -l)
  echo "[$label] outside $(basename "$bed"): $((before - after)) of $before record(s) removed"
}
