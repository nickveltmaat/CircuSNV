# *Circu*SNV
#### A Scalable and High-Sensitivity Variant Calling Pipeline
CircuSNV is a robust and scalable bash pipeline designed for the sensitive detection of somatic single-nucleotide variants (SNVs) and small insertions/deletions (indels) from targeted regions. The pipeline is optimized for ultra-low variant allele frequency (VAF) detection and is particularly useful in applications such as minimal residual disease (MRD) monitoring, using *circu*lating tumor DNA (ctDNA) from liquid biopsies. CircuSNV integrates multiple variant calling algorithms, including [VarDict](https://pubmed.ncbi.nlm.nih.gov/27060149/), [LoFreq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/), [Mutect2](https://www.biorxiv.org/content/10.1101/861054v1) & [SiNVICT](https://pubmed.ncbi.nlm.nih.gov/27531099/), applies stringent quality filters, optionally generates a Panel of Normals (PoN), and merges intersected calls. Subsequent processes include matched-control filtering, [UNMET](https://doi.org/10.1093/nar/gkad1140) score exclusion, multi-nucleotide variant (MNV) calling, and functional annotation via OpenCravat. The final outputs are organized per sample into raw, intermediate, and final annotated result directories.
 
#
## Prerequisites:
 * A process-ready BAM file (or a directory with BAM files), e.g. pre-processed following [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
 * [GATK](https://gatk.broadinstitute.org/hc/en-us) > [4.3.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.3.0.0)
 * [Vardict Java](https://bioconda.github.io/recipes/vardict-java/README.html) > 1.8.2 (`teststrandbias.R` and  `var2vcf_valid.pl` need to be in $PATH)
 * [SAMtools](http://www.htslib.org/) > [1.16.1](http://www.htslib.org/download/)
 * BCFtools > [1.11](http://www.htslib.org/download/)
 * HTSlib > [1.11](http://www.htslib.org/download/)
 * Python > [3.9](https://www.python.org/downloads/release/python-390/)
 * R > [4.0.3](https://cran.r-project.org/bin/windows/base/)
 * Perl
 * IGV and BEDTools (OPTIONAL: Only needed if UNMET-based filtering is toggled 'true' using -U)
 * LoFreq* and SiNVICT with their respective prerequisites are pre-built, more info on this in [__'Installation'__](https://github.com/nickveltmaat/SNVcaller/blob/main/README.md#installation)

#
## Workflow Description
#### 1. Data & Command Line Input:
- Read and validate required (`-I`, `-R`, `-L`, `-O`) and optional flags (`-V`, `-D`, `-C`, `-P`, `-Q`, `-M`, `-S`, `-H`, `-X`, `-T`, `-U`)
- `-I` Accepts aligned sequencing data in `.bam` format or a directory containing multiple files. Index `.bai` files will be created using SAMtools if they are not already present in the same folder. 
- Depending on if `-T` and `-X` are provided and `-I` is a directory containing multiple `.bam` files, the order of sample processing will be determined so that matched control files will be processed first. Additionally, if `-H` is true, baseline samples need to be analyzed fist.

#### --- OPTIONAL --- Panel of Normals (PoN) generation:
- If `-P` directs to a directory containing healthy control `.bam` files, variant calling will first be performed on these samples (as described below) to generate tool-specific blacklists, before starting on the tumor samples.
- All blacklists into a unified PoN text file, by keeping variants that are found in at least 2 variant callers, or with 1 variant caller but in at least 20% of the healthy control samples. 
- If `-P` directs to a PoN text file (pre made or user defined), this text file is used for downstream filtering. 

#### 2. Variant Calling with Multi-Algorithm Integration:
- VarDict: Sensitive detection of low-VAF variants.
- Mutect2: Advanced somatic variant caller from GATK.
- LoFreq: A high-specificity caller that adjusts for sequencing error models to refine calls.
- SiNVICT: Optimized for ultra-low frequency detection.

#### 3. Initial Variant Filtering:
- Per caller, variants are filtered based on thresholds for read depth (`-D`), minimum variant-supporting reads (`-M`), Variant Allele Frequency (VAF, `-V`), strand bias (`-S`), and base quality (`-Q`) metrics.
- Variants of all callers are merged, variants called with at least `-C` of these tools will be used for futher analysis.
- PoN filter (Optional). Filters out blacklisted mutations as defined earlier using the `-P` argument if provided.
- Matched control filter (Optional). If `-T` is true, variants found in healthy samples are filtered out using [matched_filtering.py](https://github.com/nickveltmaat/CircuSNV/blob/main/matched_filtering.py). This script filters variants by identifying those that are either unique to the tumor sample or show significantly higher variant allele frequency (VAF) in the tumor compared to the matched normal, using criteria like Fisher’s exact test (p < 0.05), a minimum VAF difference (≥ 0.02), and a relative VAF ratio (≥ 3). Variants with high VAF in the control (germline) are excluded.
- UNMET filter (Optional). If `-U` is true, variants with an unmet score of > 0.93 are filtered out, by automatically downloading the [10-bp-binned UNMET scores](http://web16.kazusa.or.jp/data/hg38/hg38_genome.UNMET_bin10bp.tdf) and converting it to a bed file.

#### --- OPTIONAL --- MNV calling:
- If `-H` is set to true, MNV calling will be performed on unfiltered SNVs using an in-house developed MNV calling algorithm called [MNVista](https://github.com/lasillje/MNVista). See that page for a more detailed explanation. 
- If `-X` is provided, baseline will be analyzed first. MNV calling in follow-up samples will then be restricted to SNVs that were previously identified as part of an MNV in the baseline.

#### 4. Initial annotation & Further filtering:
- OpenCravat is used to annotate variants with functional impact scores and population allele frequencies.
- Integrates databases such as [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [dbSNP](https://www.ncbi.nlm.nih.gov/snp/), [COSMIC](https://cancer.sanger.ac.uk/cosmic), [gnomAD](https://gnomad.broadinstitute.org/) and many more for pathogenicity assessment.
- [post_filtering.py](https://github.com/nickveltmaat/CircuSNV/blob/main/post_filtering.py) is used to filter annotated variant data by removing variants with low functional impact (e.g., synonymous, UTR, intronic, upstream/downstream) and common population variants.
- It outputs a filtered variant list for re-annotation, a summary of filtering steps, and an Excel file of removed variants.
- If `-U`is true, variants found in repeat regions will addtionally be filtered.

#### 6. Output & Reporting:
- The final list is re-annotated, ending up with an annotaed consensus variant set in both `.vcf` and excel, while also keeping per-caller VCF files.
- Both provide annotated variants with key metrics for downstream analyses.
- Per sample, subdirectories are created to store all the data generated by CircuSNV.
```
output_dir(-O)/sample_id/
├── raw-data/                 # per-tool filtered VCFs
├── intermediate-files/       # merged & initial filtered lists
├── filtering-info/           # logs, list of removed variants from post_filtering.py
└── final-list/               # annotated final variants in excel and vcf format. 
```

See the [`Pipeline flow diagram`](https://github.com/nickveltmaat/CircuSNV/blob/main/flowchart_pipeline2.drawio.png) for a rough visual representation of the pipeline.

#
## Installation
**1. Clone the repo**

```
git clone https://github.com/nickveltmaat/CircuSNV
```

**2. Set working directory to the repo**

```
cd /path/to/CircuSNV
```

**3. Create python virtual environment (env)**

```
python3 -m venv ./env
```

**4. Install required packages in env with pip3**

```
source ./env/bin/activate
pip3 install pandas
pip3 install scipy
pip3 install glob
pip3 install xlrd
pip install open-cravat
oc module install-base
oc module ls -a -t annotator  (this generartes a list of available annotators that can be downloaded)
oc module install thousandgenomes gnomad3 repeat cgc cadd_exome clinpred cosmic mutation_assessor ...  (see https://open-cravat.readthedocs.io/en/latest/1.-Installation-Instructions.html for more detailed instructions)
deactivate
```

**5. [Download](https://drive.google.com/drive/folders/1QBt0NdPqjQU_y-A7omxoyiPfl1DL65Xn?usp=sharing) and copy the pre-built tools to `/path/to/CircuSNV/` and unzip**
 
```
unzip ./tools.zip
```

**5. Export paths to pre-built tools**
 
```
export PATH=$PATH:/path/to/CircuSNV/tools/lofreq/src/lofreq/ #Lofreq
export PATH=$PATH:/path/to/CircuSNV/tools/lofreq/src/scripts/ #Lofreq 2?
export PATH=$PATH:/path/to/CircuSNV/tools/sinvict/bam-readcount/build/bin/ #BamReadcount
export PATH=$PATH:/path/to/CircuSNV/tools/sinvict/ #SiNVICT
```
 
**7. Adjust paths in [CircuSNV.sh](https://github.com/nickveltmaat/CircuSNV/blob/main/CircuSNV.sh)**
 
```
cd to CircuSNV installation directory
```
`TIP:  Also check the hardcoded versions of tools loaded via 'module load'. Change them if needed. `
 
#
## Usage
Once all tools and pre-requisites are installed correctly, the pipeline can be called using: 
```
bash ./CircuSNV.sh -I /path/sample.bam -R homo_sapiens.fasta -L targets.bed -O ./results/ OPTIONAL_ARGUMENTS
```

By default, the script requires 12 threads and 12 GB of ram. This can be changed in the CircuSNV.sh code (hard-coded). 

| Flag | Parameter                  | Description                                                                              | Default   |
|------|----------------------------|------------------------------------------------------------------------------------------|-----------|
| `-I` | `<input BAM or directory>` | Path to a single BAM file or directory of BAMs.                                          | **Required** |
| `-R` | `<reference.fasta>`        | Reference genome FASTA (indexed).                                                        | **Required** |
| `-L` | `<regions.bed>`            | BED file containing target regions.                                                      | **Required** |
| `-V` | `<float>`                  | Minimum Variant Allele Frequency (VAF).                                                  | `0.001`   |
| `-D` | `<int>`                    | Minimum read depth.                                                                      | `100`     |
| `-C` | `<int>`                    | Minimum number of tools calling the variant (intersection count).                        | `1`       |
| `-P` | `<PoN directory>`          | Directory of normal BAMs for PoN generation (or pre‑made blacklist.txt).                 | *off*     |
| `-Q` | `<int>`                    | Minimum base quality threshold.                                                          | `25`      |
| `-M` | `<int>`                    | Minimum reads supporting a mutation.                                                     | `3`       |
| `-S` | `<float>`                  | Strand‑bias threshold (min fraction of forward vs. reverse reads).                       | `0.1`     |
| `-H` | `<true/false>`             | Enable multi‑nucleotide variant (MNV) phasing and calling.                               | `false`   |
| `-X` | `<crosslink table.tsv>`    | TSV linking sample names, timepoints, and corresponding baselines for matched/MNV modes. | *off*     |
| `-T` | `<true/false>`             | Matched‑control subtraction logic (requires `-X`).                                       | `false`   |
| `-U` | `<true/false>`             | Use UNMET score filtering to exclude error‑prone and repeat regions.                     | `false`   |
| `-O` | `<output directory>`       | Main directory for all sample outputs.                                                   | `./output`|
| `-h` |                            | Display help message and exit.                                                           |           |

#### Some examples:
- Directory of BAMs with PoN
```
bash CircuSNV.sh \
  -I /data/tumors/ \
  -R /refs/hg38.fasta \
  -L /refs/targets.bed \
  -P /data/normals/ \
  -O /results/ \
  -V 0.005 -D 150 -C 2
```
- Matched‑mode with MNV and UNMET filtering
```
bash CircuSNV.sh \
  -I /data/tumors/ \
  -R hg38.fasta \
  -L regions.bed \
  -X crosslink_table.tsv \
  -T true \
  -H true \
  -U true \
  -O analysis_output/
```

#### Pre-made PoN example
As mentioned, `-P` can be a pre-made `.txt` file rather than a directory. Make sure it is in the following format (no headers): 
| chr1 | 2556714 | A   | G   |
|------|---------|-----|-----|
| chr1 | 2562891 | G   | A   |
| chr1 | 9710455 | G   | A   |
| chr1 | 9710456 | A   | G   |
| chr1 | 9710457 | T   | C   |
| chr1 | 9710458 | G   | GC  |


#### CrosslinkTable example
When providing `-X`, it needs to be a `.tsv` file in the following format: samplenames without file extension

| samplename | timepoint  | corresponding_baseline | matched_control |
|------------|------------|-------------------------|-----------------|
| SAMPLE_A   | baseline   |                         | CONTROL_X       |
| SAMPLE_B   | follow-up  | SAMPLE_A                | CONTROL_X       |
| SAMPLE_C   | follow-up  | SAMPLE_A                | CONTROL_X       |
| SAMPLE_D   | follow-up  | SAMPLE_A                | CONTROL_X       |
| SAMPLE_E   | baseline   |                         | CONTROL_Y       |
| SAMPLE_F   | follow-up  | SAMPLE_E                | CONTROL_Y       |
| SAMPLE_G   | follow-up  | SAMPLE_E                | CONTROL_Y       |
| SAMPLE_H   | follow-up  | SAMPLE_E                | CONTROL_Y       |

`NOTE: If -T is not provided, the matched_control column is not required`
#
