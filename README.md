# CircuSNV
A pipeline to detect SNVs and other short somatic variants. This tool can work with low variant allele frequencies, ideal for analyzing circulating tumor DNA (ctDNA) without the need of matched controls. 

#### **Pipeline to call SNVs with 4 tools ([`VarDict`](https://pubmed.ncbi.nlm.nih.gov/27060149/), [`LoFreq`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3526318/), [`Mutect2`](https://www.biorxiv.org/content/10.1101/861054v1) & [`SiNVICT`](https://pubmed.ncbi.nlm.nih.gov/27531099/))**


### Prerequisites:
 * A process-ready BAM file (or a directory with BAM files), e.g. pre-processed following [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
 * [GATK](https://gatk.broadinstitute.org/hc/en-us) > [4.2.4.1](https://github.com/broadinstitute/gatk/releases/tag/4.2.4.1)
 * [Vardict Java](https://bioconda.github.io/recipes/vardict-java/README.html) > 1.8.2
 * [SAMtools](http://www.htslib.org/) > [1.16.1](http://www.htslib.org/download/)
 * BCFtools > [1.11](http://www.htslib.org/download/)
 * HTSlib > [1.11](http://www.htslib.org/download/)
 * Python > [3.9](https://www.python.org/downloads/release/python-390/)
 * R > [4.0.3](https://cran.r-project.org/bin/windows/base/)
 * LoFreq* and SiNVICT with its prerequisites are pre-built, more info on this in* [__'Installation'__](https://github.com/nickveltmaat/SNVcaller/blob/main/README.md#installation)

## General Description
This is a pipeline made to reliably generate calls for somatic mutations in Low Variant Allele Frequencies (VAF) samples in specific regions, such as NGS data from cfDNA. This is done by analyzing `.BAM` files with 4 different tools (`VarDict`, `LoFreq`, `Mutect2` & `SiNVICT`). The pipeline will output variants that are called with at least an `x` amount of tools (this can be set from 1-4). Of course, the higher the number, the lower False Positive call rate, the higher the reliability of the call, but also a higher chance of missing relevant somatic variants. 

The general workflow in the pipeline is as follows: 

Firstly, a panel of normals (PoN, a blacklist to filter mutations) can be generated if the argument `-P` is passed, leading to a directory containing normal/healty control samples. This is optional. To generate this, all normal samples are being analayzed with the 4 tools, in the exact same way as the tumor samples, to generate a list of personal variants / SNPs and technical artifacts. All variants found in all of the normal samples, called with at least one of the tools will be included. This blacklist will later be used as a filter.
Ather generating the PoN, the 4 tools will run in parralel, generating raw data for each SNV caller. Every tool has `.vcf` files as output, which are then gunzipped and indexed, so they can be filtered on Read Depth (RD, `-D`), Variant Allele Frequency (VAF, `-V`) and Read Depth of Mutant allele (MRD, `-M`). Subsequently, with all `.vcf` files processed, the variants can be merged on overlapping variants, keeping the RD, VAF and MRD parameters intact. All variants called with `x` or more tools will be kept. Since the tools can call variants in different ways, sometimes a variant is part of a larger variant (MNV, multinucleotide variant). These are therefore duplicate and need to be removed, a proces that a custom python script takes care of. Next, the variants are annotated and directly after filtered based on functional effect (e.g. non-synonymous variants are kept). Finally, the remaining variants will be annotated again using [openCRAVAT](https://opencravat.org/). This is a wrapper around multiple well-known annotating tools, such as [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [dbSNP](https://www.ncbi.nlm.nih.gov/snp/), [COSMIC](https://cancer.sanger.ac.uk/cosmic), [gnomAD](https://gnomad.broadinstitute.org/) and many more. All annotated mutations are saved in an excel file and as `.vcf` files. 


## Installation
**1. Clone the repo**

`git clone https://github.com/nickveltmaat/CircuSNV`

**2. Set working directory to the repo**

`cd /path/to/CircuSNV`

**3. Create python virtual environment (env)**

`python3 -m venv ./env`

**4. Install needed packages in env with pip3**

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
**5. [Download](https://drive.google.com/drive/folders/1QBt0NdPqjQU_y-A7omxoyiPfl1DL65Xn?usp=sharing) and copy the pre-built tools to `/path/to/SNVCaller/` and unzip**
 
 `unzip ./tools.zip`

## Usage
Once all tools and pre-requisites are installed correctly, the pipeline can be called with: 

`bash ./CircuSNV.sh ARGUMENTS`

**Required** arguments:
- `-I` Input:              **String**   --> example: `/path/to/input.bam` Either one-file or directory
- `-R` Reference:          **String**    --> example: `/path/to/reference.fa`
- `-L` Regions List:       **String**    --> example: `/path/to/panel.bed`
- `-D` minimum Read Depth:  **Int**       --> example: `100`
- `-V` minimum VAF:         __float *[0-1]*__  --> example: `0.004`
- `-C` minimum Calls:       __Int *[1-4]*__ --> example: `2`
- `-P` Panel of Normal:     **String** --> example: `/path/to/PoN/directory/` Optional
- `-Q` Base Quality:        **Int**  --> example: `18`
- `-M` minimum Mutant Read Depth:   **Int**  --> example: `7`

Output will be generated in `/path/to/CircuSNV/output/name_of_.bam_file/`
