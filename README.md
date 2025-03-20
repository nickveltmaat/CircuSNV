# *Circu*SNV
#### A Scalable and High-Sensitivity Variant Calling Pipeline
CircuSNV is a robust and scalable bioinformatics pipeline designed for the sensitive detection of somatic single-nucleotide variants (SNVs) and small insertions/deletions (indels) from targeted or whole-exome sequencing data. The pipeline is optimized for ultra-low variant allele frequency (VAF) detection and is particularly useful in applications such as minimal residual disease (MRD) monitoring, using *circu*lating tumor DNA (ctDNA) from liquid biopsies. CircuSNV integrates multiple variant calling algorithms, including [VarDict](https://pubmed.ncbi.nlm.nih.gov/27060149/), [LoFreq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/), [Mutect2](https://www.biorxiv.org/content/10.1101/861054v1) & [SiNVICT](https://pubmed.ncbi.nlm.nih.gov/27531099/), to ensure comprehensive variant detection while minimizing false positives, with- or without the use of matched normal samples and/or an optional panel of normal (PoN) blacklist.
 
#
### Prerequisites:
 * A process-ready BAM file (or a directory with BAM files), e.g. pre-processed following [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
 * [GATK](https://gatk.broadinstitute.org/hc/en-us) > [4.2.4.1](https://github.com/broadinstitute/gatk/releases/tag/4.2.4.1)
 * [Vardict Java](https://bioconda.github.io/recipes/vardict-java/README.html) > 1.8.2
 * [SAMtools](http://www.htslib.org/) > [1.16.1](http://www.htslib.org/download/)
 * BCFtools > [1.11](http://www.htslib.org/download/)
 * HTSlib > [1.11](http://www.htslib.org/download/)
 * Python > [3.9](https://www.python.org/downloads/release/python-390/)
 * R > [4.0.3](https://cran.r-project.org/bin/windows/base/)
 * LoFreq* and SiNVICT with its prerequisites are pre-built, more info on this in [__'Installation'__](https://github.com/nickveltmaat/SNVcaller/blob/main/README.md#installation)

## Workflow Description
#### 1. Data Input:
- Accepts aligned sequencing data in `.BAM` format.
- Requires a reference genome (e.g., GRCh38), a targeted regions `.BED` file, and optional panel-of-normals (PoN) files for filtering.

#### 2. Variant Calling with Multi-Algorithm Integration:
- VarDict: Sensitive detection of low-VAF variants.
- Mutect2: Advanced somatic variant caller from GATK.
- LoFreq: A high-specificity caller that adjusts for sequencing error models to refine calls.
- SiNVICT: Optimized for ultra-low frequency detection.
All variants called with `x` or more these tools will be used for futher analysis.

#### 3. Comprehensive Variant Filtering:
- Variants are filtered based on user-defined thresholds for read depth (DP), minimum variant-supporting reads (MRD), VAF cutoff, and strand bias constraints.
- Optional PoN filtering using sequencing data from normal samples removes systematic sequencing artifacts and SNPs. PoN consists of variants found in a group of normal samples.

#### 4. Annotation & Interpretation:
- Uses OpenCravat to annotate variants with functional impact scores and population allele frequencies.
- Integrates databases such as [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [dbSNP](https://www.ncbi.nlm.nih.gov/snp/), [COSMIC](https://cancer.sanger.ac.uk/cosmic), [gnomAD](https://gnomad.broadinstitute.org/) and many more for pathogenicity assessment.

#### 5. Further filtering
- Variants are filtered based on functional effect (e.g. non-synonymous variants are kept)
- Variants are filtered based on population allele frequency

#### 6. Output & Reporting:
- Generates per-caller VCF files and a consensus variant set in both `.vcf` and excel.
- Both provide annotaed variants with key metrics for downstream analysis.

Identification of multinucleotide variants (MNVs) is in development and will be added shortly. 
See [`S1_SNV_pipeline.png`](https://github.com/nickveltmaat/CircuSNV/blob/main/S1_SNV_pipeline.png) for a rough flowchart of the pipeline.


## Installation
**1. Clone the repo**

`git clone https://github.com/nickveltmaat/CircuSNV`

**2. Set working directory to the repo**

`cd /path/to/CircuSNV`

**3. Create python virtual environment (env)**

`python3 -m venv ./env`

**4. Install required packages in env with pip3**

```
source ./env/bin/activate
pip3 install pandas
pip3 install glob
pip3 install xlrd
pip install open-cravat
oc module install-base
oc module ls -a -t annotator  (this generartes a list of available annotators that can be downloaded)
oc module install clinvar cosmic dbsnp ...  (see https://open-cravat.readthedocs.io/en/latest/1.-Installation-Instructions.html for more detailed instructions)
deactivate
```

**5. [Download](https://drive.google.com/drive/folders/1QBt0NdPqjQU_y-A7omxoyiPfl1DL65Xn?usp=sharing) and copy the pre-built tools to `/path/to/CircuSNV/` and unzip**
 
`unzip ./tools.zip`

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
source ../../../env/bin/activate in pythonscript() and annotate() functions
```
 

## Usage
Once all tools and pre-requisites are installed correctly, the pipeline can be called with: 

`bash ./CircuSNV.sh ARGUMENTS`

**Required** arguments:
- `-I` Input:              **String**   --> example: `/path/to/input.bam` Either one-file or directory
- `-R` Reference:          **String**    --> example: `/path/to/reference.fa`
- `-L` Regions List:       **String**    --> example: `/path/to/panel.bed`
- `-D` minimum Read Depth:  **Int**       --> example: `100`
- `-V` minimum VAF:         __float *[0-1]*__  --> example: `0.001`
- `-C` minimum Calls:       __Int *[1-4]*__ --> example: `2`
- `-P` Panel of Normal:     **String** --> example: `/path/to/PoN/directory/` or `path/to/PoN/pre-made-list.txt` Optional
- `-Q` Base Quality:        **Int**  --> example: `18`
- `-M` minimum Mutant Read Depth:   **Int**  --> example: `7`
- `-F` Filter Mode:         **String**  --> `Combined` or `PerTool`
- `-S` Strand bias          __float *[0-1]*__ example `0.2` Lowest fraction of FW reads vs RV reads. (SB = 0.2 --> FW = 1, RV = 5 and vice versa)


Output will be generated in `/path/to/CircuSNV/output/name_of_.bam_file/`
