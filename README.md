# *Circu*SNV
#### A scalable, high-sensitivity variant calling pipeline for cfDNA

CircuSNV detects somatic single-nucleotide variants (SNVs) and small insertions/deletions
(indels) in targeted sequencing data of *circu*lating tumour DNA (ctDNA), without matched
controls. It is built for ultra-low variant allele frequencies (VAF), for example for minimal
residual disease (MRD) monitoring with liquid biopsies.

CircuSNV runs four variant callers in parallel ([VarDict](https://pubmed.ncbi.nlm.nih.gov/27060149/),
[Mutect2](https://www.biorxiv.org/content/10.1101/861054v1), [LoFreq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/)
and [SiNVICT](https://pubmed.ncbi.nlm.nih.gov/27531099/)), filters and merges their calls,
optionally builds and applies a panel of normals (PoN), optionally filters on the
[UNMET](https://doi.org/10.1093/nar/gkad1140) score, optionally calls multi-nucleotide variants
(MNVs) with [MNVista](https://github.com/lasillje/MNVista), and annotates the result with
OpenCRAVAT.

## Workflow

1. **Input and checks.** All options are validated and, unless `--skip-checks` is given, the
   pipeline checks the inputs, the Python environment, MNVista, the environment modules and
   every required command before any work starts. A BAM index is created when it is missing.
2. **Panel of normals (optional, `-P`).**
   - A directory with normal BAMs: the four callers run on every normal; per tool the normals
     are merged; variants found by at least 2 tools form the blacklist, and variants found by
     1 tool are added when they occur in at least 20% of the normals.
   - A file: a pre-made blacklist (`chrom pos ref alt`, tab-separated, no header).
3. **Variant calling.** VarDict, Mutect2, LoFreq and SiNVICT run in parallel per sample, each
   with up to three attempts. Per caller, variants are filtered on strand bias (`-S`), read
   depth (`-D`), VAF (`-V`), supporting reads (`-M`) and base quality (`-Q`), then
   normalised and atomised.
4. **Merging.** Variants called by at least `-C` tools are kept, the PoN blacklist is
   applied, and the depth, VAF and supporting reads of each tool are added back.
5. **UNMET filter (optional, `-U true`).** Variants with an UNMET score of 0.90 or higher are
   removed. The 10-bp binned
   [UNMET scores](http://web16.kazusa.or.jp/data/hg38/hg38_genome.UNMET_bin10bp.tdf) are
   downloaded and converted on first use.
6. **MNV calling (optional, `-H true`).** MNVista runs on the unfiltered calls of all tools.
   With a crosslink table (`-X`), baselines are analysed first and each follow-up is analysed
   only at the SNVs that formed an MNV in its baseline. An MNV panel of normals is used only
   when it is given with `--mnv-pon`; the PoN of `-P` is for the SNV calls.
7. **Annotation and post-filtering.** OpenCRAVAT annotates the variants;
   `scripts/post_filtering.py` removes synonymous, UTR, intronic and 2 kb up/downstream
   variants, variants with a population frequency above 0.01% in 1000 Genomes or gnomAD, and,
   with `-U true`, variants in repeat regions.
8. **Matched-control read counts (optional, `-T true`).** For each remaining variant the
   reference and alternative read counts in the matched control (from the crosslink table)
   are added.
9. **Final annotation.** The remaining variants are annotated again, to an Excel file and a
   VCF. The second annotation step is deliberate: post-filtering loses the layout of the first
   annotation.

See the [pipeline flow diagram](flowchart_pipeline2.drawio.png) for a rough overview.

## Requirements

- A Linux cluster with Slurm and environment modules (Lmod/EasyBuild). The modules the
  pipeline loads are listed, with their default names, in
  [config/circusnv.example.conf](config/circusnv.example.conf).
- Process-ready BAM files, for example pre-processed following the
  [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).
- A reference FASTA (hg38; the annotation and UNMET steps assume hg38) with its `.fai` index
  and GATK `.dict`.
- Tools: GATK 4.3, VarDict Java (`vardict-java`, `teststrandbias.R` and `var2vcf_valid.pl` on
  `PATH`), SAMtools, BCFtools, HTSlib, LoFreq, bam-readcount, SiNVICT, R, Perl, `bc`;
  for `-U true` also BEDTools, `wget` and IGV (`igvtools`); for `-H true` MNVista.
- A Python 3 virtual environment `env/` in the installation directory with the packages from
  [requirements.txt](requirements.txt) and the OpenCRAVAT annotators `cgc cadd_exome clinpred
  cosmic mutation_assessor thousandgenomes gnomad3 repeat`.

## Installation

```
git clone https://github.com/nickveltmaat/CircuSNV
cd CircuSNV
python3 -m venv ./env
source ./env/bin/activate
pip install -r requirements.txt
oc module install-base
oc module install cgc cadd_exome clinpred cosmic mutation_assessor thousandgenomes gnomad3 repeat
deactivate
```

See the [OpenCRAVAT installation instructions](https://open-cravat.readthedocs.io/en/latest/1.-Installation-Instructions.html)
for details.

Pre-built LoFreq and SiNVICT (with bam-readcount) can be
[downloaded](https://drive.google.com/drive/folders/1QBt0NdPqjQU_y-A7omxoyiPfl1DL65Xn?usp=sharing)
as `tools.zip`. Unzip it in the installation directory and add the tools to your `PATH`
(for example in `~/.bashrc`, which Slurm jobs read through `--get-user-env`):

```
export PATH=$PATH:/path/to/CircuSNV/tools/lofreq/src/lofreq/
export PATH=$PATH:/path/to/CircuSNV/tools/lofreq/src/scripts/
export PATH=$PATH:/path/to/CircuSNV/tools/sinvict/bam-readcount/build/bin/
export PATH=$PATH:/path/to/CircuSNV/tools/sinvict/
```

For MNV calling, place the MNVista executable in the installation directory, or point to it
in a config file. An MNV panel of normals is optional: give it with `--mnv-pon`.

## Usage

Submit from the directory where relative paths should be resolved:

```
sbatch /path/to/CircuSNV/CircuSNV.sh -I sample.bam -R hg38.fasta -L targets.bed -O results/
```

or run directly with `bash CircuSNV.sh ...`. The `#SBATCH` header of `CircuSNV.sh` requests
16 CPUs and 16 GB; change it there or with `sbatch` options. Start with `--dry-run` to check
everything and see the plan without running anything.

| Option | Value | Description | Default |
| --- | --- | --- | --- |
| `-I`, `--input` | path | BAM file or directory with BAM files | required |
| `-R`, `--reference` | file | reference FASTA | required |
| `-L`, `--regions` | file | BED file with target regions | required |
| `-V`, `--min-vaf` | number | minimum variant allele frequency | `0.001` |
| `-D`, `--min-depth` | integer | minimum read depth | `100` |
| `-C`, `--min-callers` | 1-4 | minimum number of callers | `1` |
| `-Q`, `--min-base-quality` | integer | minimum base quality | `25` |
| `-M`, `--min-alt-reads` | integer | minimum reads supporting a mutation | `3` |
| `-S`, `--strand-bias` | number | the forward/reverse ratio of the supporting reads must lie between this value and its inverse | `0.9` |
| `-P`, `--pon` | path | directory with normal BAMs, or a pre-made blacklist file | off |
| `-U`, `--unmet` | `true`/`false` | UNMET score and repeat-region filtering | `false` |
| `-H`, `--mnv` | `true`/`false` | MNV calling with MNVista | `false` |
| `-B`, `--mnv-prior` | number | Bayesian prior for MNVista | `0.5` |
| `--mnv-pon` | file | MNV panel of normals for MNVista (its `-C` option), independent of `-P` | off |
| `-X`, `--crosslink-table` | file | crosslink table, see below | off |
| `-T`, `--matched-control` | `true`/`false` | add matched-control read counts (needs `-X` and a directory `-I`) | `false` |
| `-O`, `--output` | directory | main output directory | `./output` |
| `--threads` | integer | threads per tool | `SLURM_CPUS_PER_TASK`, else `12` |
| `--config` | file | settings file, see [config/circusnv.example.conf](config/circusnv.example.conf) | none |
| `--dry-run` | | check everything and show what would run, then stop | |
| `--resume` | | skip samples that finished in an earlier run into the same output directory | |
| `--cleanup` | | remove `isec/` and MNV intermediates after each finished sample | |
| `--skip-checks` | | do not run the pre-run checks | |
| `-h`, `--help` / `--version` | | show help / version | |

Relative paths are taken from the directory the pipeline is started from. Without
`--resume`, a run stops before starting when the output directory already contains results
for one of its samples. If a sample fails, the other samples still run; the run then ends
with a non-zero exit status and lists the failed samples, which `--resume` redoes.

### Examples

Directory of BAMs with a panel of normals:

```
sbatch CircuSNV.sh -I /data/tumors/ -R /refs/hg38.fasta -L /refs/targets.bed \
  -P /data/normals/ -O /results/ -V 0.005 -D 150 -C 2
```

Baselines and follow-ups with MNV calling, matched controls and UNMET filtering:

```
sbatch CircuSNV.sh -I /data/tumors/ -R hg38.fasta -L regions.bed \
  -X crosslink_table.tsv -T true -H true -U true -O analysis_output/
```

### Crosslink table

A tab-separated file with a header line; sample names without file extension. The
`matched_control` column is only needed with `-T true`.

| samplename | timepoint | corresponding_baseline | matched_control |
| --- | --- | --- | --- |
| SAMPLE_A | baseline | | CONTROL_X |
| SAMPLE_B | follow-up | SAMPLE_A | CONTROL_X |
| SAMPLE_C | baseline | | CONTROL_Y |
| SAMPLE_D | follow-up | SAMPLE_C | CONTROL_Y |

With a crosslink table, only the samples marked `baseline` or `follow-up` are processed:
first all baselines, then all follow-ups.

### Pre-made panel of normals

A tab-separated file without header:

| | | | |
| --- | --- | --- | --- |
| chr1 | 2556714 | A | G |
| chr1 | 9710458 | G | GC |

## Output

```
<output>/
├── <sample>/
│   ├── <sample>_<Tool>_RAW.vcf         # per caller, after the strand-bias filter
│   ├── <sample>_<Tool>.vcf.gz(.csi)    # per caller, all filters, normalised
│   ├── isec/                           # merged calls (bcftools isec), PoN-filtered calls
│   ├── mnvcalling/                     # MNV input and MNVista results (with -H true)
│   ├── intermediate-files/             # merged variant lists and the first annotation
│   ├── filtering-info/                 # variants left after each filter, removed variants
│   ├── final-list/                     # annotated final variants (Excel and VCF)
│   └── raw-data/                       # (empty, kept for compatibility)
└── .circusnv/                          # run logs and resume markers
```

A sample without variants after post-filtering gets `NO_MUTS.txt` and an empty
`final-list/`. The panel of normals is built in `PoN/` in the installation directory.

## Development and testing

The pipeline code lives in `CircuSNV.sh` (entry point) and `lib/` (bash modules); the Python
helpers are in `scripts/`. Changes are tested on the cluster by comparing the output with a
reference run: see [docs/nibbler-test.md](docs/nibbler-test.md) and
[tests/nibbler/compare_outputs.py](tests/nibbler/compare_outputs.py). Open questions about
intended behaviour are listed in [docs/known-issues.md](docs/known-issues.md); changes per
version are in [CHANGELOG.md](CHANGELOG.md).
