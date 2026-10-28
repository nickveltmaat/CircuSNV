# shellcheck shell=bash
# Settings that used to be hard-coded in the pipeline: environment modules, resource locations
# and fixed thresholds. The defaults below reproduce the original behaviour exactly, except
# MNV_PON, which is off by default (the original always used it together with -P).
# Sourced by CircuSNV.sh; uses its global settings.

apply_default_config() {
  # Environment modules (EasyBuild names as loaded by the original pipeline).
  MODULE_ZLIB="zlib/1.2.12-GCCcore-11.3.0"
  MODULE_HTSLIB="HTSlib/1.19.1-GCCcore-11.3.0"
  MODULE_BCFTOOLS="BCFtools/1.19-GCCcore-11.3.0"
  MODULE_SAMTOOLS="SAMtools/1.19.2-GCCcore-11.3.0"
  MODULE_SAMTOOLS_PON="SAMtools/1.16.1-GCCcore-11.3.0"   # used to index normal BAMs
  MODULE_R="R/4.2.2-foss-2022a-bare"
  MODULE_GATK="GATK/4.3.0.0-GCCcore-11.3.0-Java-11"
  MODULE_GCC_MNVISTA="GCC/13.3.0"
  MODULE_IGV="IGV"
  MODULE_BEDTOOLS="BEDTools"

  # Resources inside the installation directory.
  ENV_DIR="$CIRCUSNV_HOME/env"
  PON_DIR="$CIRCUSNV_HOME/PoN"
  MNVISTA="$CIRCUSNV_HOME/MNVista"
  MNV_PON=""                    # MNV panel of normals for MNVista -C; empty = off (--mnv-pon)
  UNMET_BED="$CIRCUSNV_HOME/hg38_genome.UNMET_bin10bp.bed"
  UNMET_URL_BASE="http://web16.kazusa.or.jp/data/hg38"

  # Fixed values (see docs/known-issues.md for the open decisions about them).
  UNMET_THRESHOLD="0.90"        # keep variants with an UNMET score below this value
  PON_MIN_NORMALS_PCT="20"      # single-tool PoN variants must occur in this % of normals
  MUTECT2_JAVA_MEM="15G"        # Java heap for Mutect2
}
