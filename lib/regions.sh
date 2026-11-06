# shellcheck shell=bash
# Region (BED) files: the cleaned panel BED and restricting VCFs to a BED.
# Sourced by CircuSNV.sh; uses its global settings.
#
# prepare_regions sets, once per run and before the panel of normals and the samples:
#   REGIONS_PANEL    the -L BED without header lines; the main path (and the normals of the
#                    PoN) is restricted to it
#   REGIONS_MNV      the BED the MNV input is restricted to (REGIONS_PANEL)
#   REGIONS_CALLING  the BED the tumor callers run on (the -L file as given). The normals of
#                    the PoN always run on the -L file.

# Print the data lines of a BED file: no browser, track, comment or empty lines, no CR.
bed_data_lines() {
  awk '{ sub(/\r$/, "") } !/^(browser|track|#)/ && NF' "$1"
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
  REGIONS_MNV=$REGIONS_PANEL
  REGIONS_CALLING=$Listregions
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
