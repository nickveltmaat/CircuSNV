#!/bin/bash
#SBATCH --job-name=CircuSNV_PTLD2
#SBATCH --output="CircuSNV_PTLD2.out"
#SBATCH --error="CircuSNV_PTLD2.err"
#SBATCH --time=159:59:00
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --get-user-env=L60
#SBATCH --export=NONE


###֖ Helper: print usage and exit ֖###
usage() {
  cat <<EOF

Usage: bash $0 -I <input_bam_or_dir> -R <reference.fasta> -L <regions.bed> [options]

Required arguments:
  -I   Input BAM file or directory
  -R   Reference FASTA
  -L   BED of target regions

Optional arguments:
  -V   Minimum VAF (default: 0.001)
  -D   Minimum read depth (default: [100])
  -C   Minimum number of variant calls (default: [1])
  -P   Panel-of-Normals directory
  -Q   Minimum base quality (default: [25])
  -M   Minimum reads supporting a mutation (default: [3])
  -S   Strand-bias threshold (default: [0.1])
  -H   Enable MNV calling (true/[false])
  -X   Crosslink table TSV
  -T   Use matched control logic (true/[false])
  -U   Use UNMET score to filter error-prone variants (true/[false])
  -O   Main output directory (default: [./output])
  -h   Show this help message and exit

EOF
  exit 1
}

###֖ Default (optional) parameter values ֖###
VAF=0.001            # minimum variant allele frequency
RDP=100              # minimum read depth
Calls=1              # minimum number of variant-calls
PoN=""               # panel-of-normals path (empty = off)
Qual="25"            # minimum base quality
MRD="3"              # minimum reads supporting a mutation
SB="0.1"             # strand-bias threshold
MNV="false"          # enable multi-nucleotide-variant logic?
CrosslinkTable=""    # TSV file for crosslinking (empty = off)
Matched="false"      # use matched control logic?
UNMET="false"        # use UNMET score filtering
OutputDir="./output"

###֖ Parsing parameters ֖###
while getopts "R:L:I:V:D:C:P:Q:B:M:S:H:X:T:U:O:h:" arg; do
  case $arg in
    R) Reference=$OPTARG;;      # "/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly38.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation_hg38.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
    Q) Qual=$OPTARG;;           # BaseQuality (25)
    M) MRD=$OPTARG;;            # Minimum reads supporting a mutation
    S) SB=$OPTARG;;             # Strand bias. Lowest fraction of FW reads vs RV reads. (SB = 0.2 --> FW = 1, RV = 5 and vice versa)
    H) MNV=$OPTARG;;            # if "true", apply MultiNucleotide Variant calling (detect variants in pHase)
    X) CrosslinkTable=$OPTARG;; # tab separated table (.tsv file) with 'samplename', 'timepoint' and 'corresponding_baseline' columns
    T) Matched=$OPTARG;;        # if "true", apply matched_control logic
    U) UNMET=$OPTARG;;          # if "true", apply UNMET score to filter hard to sequence regions
    O) OutputDir=$OPTARG;;      # Sets the path to the main directory in which the data will be outputted. 
    h) usage;;
    *) usage;;
  esac
done

#֖ Bail out if any required flags are missing
check_required_args() {
  local miss=()
  for var in Reference Listregions Inputbam; do
    [[ -z "${!var:-}" ]] && miss+=("$var")
  done
  if (( ${#miss[@]} )); then
    echo "ERROR: Missing required arguments: ${miss[*]}"
    usage
  fi
}

check_required_args

#֖ Initiating settings
cd /groups/umcg-pmb/tmp02/projects/hematopathology/Nick/CircuSNV/ #CHANGE THIS
envpath="./env"
SBU=$(echo "scale=2; 1 / $SB" | bc)

echo -e "\nSettings: \nInput bam(s): $Inputbam\nReference: $Reference\nPanel: $Listregions\nMinimum VAF: $VAF\nMinimum Read Depth: $RDP\nMinimum Base Quality: $Qual\nMinimum Overlapping Calls: $Calls\nPoN: $PoN\nStrand Bias: $SB\nMultiNucleotide Variant Calling: $MNV\nCrosslink Table: $CrosslinkTable\nMatched Control Logic: $Matched\nUNMET & repeat filtering: $UNMET\nMain output directory: $OutputDir"


# Function to load modules
load_modules() {
  module load zlib/1.2.12-GCCcore-11.3.0
  module load HTSlib/1.19.1-GCCcore-11.3.0
  module load BCFtools/1.19-GCCcore-11.3.0
}

# Function to run Python script
pythonscript() {
  source "$envpath/bin/activate"
  python3 "$@"
  deactivate
}

# Function to annotate with OpenCravat
annotate() { #1 input.txt (tab delim snv data)   #2 Output filenames    #3 & 4 output format
  echo -e "annotating $1 with OpenCravat..."
  source "$envpath/bin/activate"
  oc run $1 \
    -l hg38 \
    -n $2 --silent --mp 12 \
    -a cgc cadd_exome clinpred cosmic mutation_assessor thousandgenomes gnomad3 repeat \
    -t $3 $4
  deactivate 
}

# Function to normalize, atomize, compress, and index filtered VCF files
process_vcf() {
  local input_vcf=$1
  local output_vcf=$2

  bcftools norm -O v -m -any -f "$Reference" -a --atom-overlaps . -c w "$input_vcf" 2>/dev/null | bcftools norm -d all > "$output_vcf" 2>/dev/null
  bgzip "$output_vcf"
  bcftools index "${output_vcf}.gz"
}


################################################### 
# Functions to run each tool, filter on Read Dept, VAF and Mutant allele dept & creating .vcf file

Vardict() { #1: Input .bam file   #2 Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with VarDict..."
  (
    module load R/4.2.2-foss-2022a-bare
    load_modules
    
    vardict-java -th 12 -r 1 -G $Reference -f 0.00000000001 -N $base_name -I 20 -L 20 -Y 20 --nosv --deldupvar -q $Qual -Q $Qual -b $1 -c 1 -S 2 -E 3 -g 4 $Listregions | \
    teststrandbias.R | var2vcf_valid.pl -S -N $base_name -E -f 0.00000000001 > $2${base_name}_VarDict_RAW_inter.vcf
    
    bcftools filter -O v -i "(FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) >= ${SB} & (FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) <= ${SBU} & (INFO/ODDRATIO[0]) > 0 & INFO/PMEAN < 67" \
    $2${base_name}_VarDict_RAW_inter.vcf 2>/dev/null > $2${base_name}_VarDict_RAW.vcf
    
    bgzip -c $2${base_name}_VarDict_RAW.vcf > $2${base_name}_VarDict_inter.vcf.gz && bcftools index $2${base_name}_VarDict_inter.vcf.gz && cp $2${base_name}_VarDict_inter.vc* $2/mnvcalling
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/VD > ${MRD}" $2${base_name}_VarDict_inter.vcf.gz > $2${base_name}_VarDict_inter2.vcf.gz && \
    process_vcf "$2${base_name}_VarDict_inter2.vcf.gz" "$2${base_name}_VarDict.vcf"
    rm -rf $2${base_name}*inter*.*
    echo -e "VarDict DONE for $base_name\n"
  )
}

######
Mutect2() { #1: Input .bam file    #2: Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with Mutect2..."
  (
    module load GATK/4.3.0.0-GCCcore-11.3.0-Java-11
    load_modules
    
    gatk --java-options -Xmx12G Mutect2 --QUIET --verbosity ERROR --create-output-variant-index false --native-pair-hmm-threads 12 --min-base-quality-score $Qual --callable-depth 10 --minimum-allele-fraction 0.00000000001 -R $Reference -I $1 -O $2${base_name}_Mutect2_RAW.vcf -L $Listregions $3 && \
    gatk FilterMutectCalls --QUIET --verbosity ERROR --create-output-variant-index false --min-allele-fraction 0.00000000001 -R $Reference -V $2${base_name}_Mutect2_RAW.vcf -O $2${base_name}_M2_Fstats1.vcf
    bcftools filter -O v -i "(FORMAT/SB[0:2]/FORMAT/SB[0:3]) >= ${SB} & (FORMAT/SB[0:2]/FORMAT/SB[0:3]) <= ${SBU}" $2${base_name}_M2_Fstats1.vcf > $2${base_name}_M2_Fstats.vcf
    
    grep -v '^##contig=<ID=chrUn\|^##contig=<ID=HLA' $2${base_name}_M2_Fstats.vcf | grep -v '_alt' | grep -v '_random' > $2${base_name}_M2_F2stats.vcf && cp $2${base_name}_M2_F2stats.vc* $2/mnvcalling
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/AD[0:1] > ${MRD}" $2${base_name}_M2_F2stats.vcf > $2${base_name}M2_F3stats.vcf.gz && \
    process_vcf "$2${base_name}M2_F3stats.vcf.gz" "$2${base_name}_Mutect2.vcf"
    
    rm -rf $2${base_name}*tats.*
    echo -e "Mutect2 DONE for $base_name\n"
  )
}


######
Lofreq() { #1: Input .bam file    #2: Output folder
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with LoFreq..."
  (
    load_modules
    
    lofreq call-parallel --pp-threads 12 --no-default-filter --call-indels -C 10 -f $Reference -q $Qual -Q $Qual -o $2${base_name}_LoFreq_inter.vcf -l $Listregions $1 && \
    awk '!seen[$0]++' $2${base_name}_LoFreq_inter.vcf > $2${base_name}_LoFreq_inter_RAW.vcf && \
    
    bcftools filter -O v -i "(INFO/DP4[2]/INFO/DP4[3]) >= ${SB} & (INFO/DP4[2]/INFO/DP4[3]) <= ${SBU}" $2${base_name}_LoFreq_inter_RAW.vcf > $2${base_name}_LoFreq_RAW.vcf

    bgzip -c $2${base_name}_LoFreq_RAW.vcf > $2${base_name}_LoFreq_inter.vcf.gz && bcftools index $2${base_name}_LoFreq_inter.vcf.gz && cp $2${base_name}_LoFreq_inter.vcf.* $2/mnvcalling

    bcftools filter -O z -i "INFO/DP > ${RDP} & INFO/AF > ${VAF} & INFO/DP4[2]+INFO/DP4[3] > ${MRD}" $2${base_name}_LoFreq_inter.vcf.gz > $2${base_name}_LoFreq_inter2.vcf.gz  && \
    process_vcf "$2${base_name}_LoFreq_inter2.vcf.gz" "$2${base_name}_LoFreq.vcf"   
    
    rm -rf $2${base_name}_LoFreq_inter*.*
    echo -e "LoFreq DONE for $base_name\n"
  )
}

######
Sinvict() { #1: Input .bam file    #2: Output folder
  InputSinvict=$1
  local base_name=$(basename $1 .bam)
  echo "Analyzing ${base_name} with bam-readcount && SiNVICT..."
  
  mkdir -p $2/${base_name}_SiNVICT_inter_bamreadcount/
  
  # Define the number of chunks for parallel processing
  N=12
  
  # Split the BED file into N chunks (excluding the header lines) and add the header lines back to each chunk
  header_lines=$(grep -P "^browser|^track" $Listregions | wc -l)
  tail -n +$((header_lines + 1)) $Listregions | split -l $((($(wc -l < $Listregions) - header_lines) / N + 1)) -d -a 2 --additional-suffix=.bed - $2/bed_chunk_

  for chunk_file in $2/bed_chunk_*.bed; do
    (grep -P "^browser|^track" $Listregions; cat $chunk_file) > ${chunk_file}.tmp && mv ${chunk_file}.tmp $chunk_file
  done
  
  # Function to run bam-readcount for a chunk
  run_bam_readcount() {
    local chunk_file=$1
    local output_file=$2
    load_modules
    bam-readcount -l $chunk_file -w 1 -b $Qual -f $Reference $InputSinvict 2>/dev/null > $output_file
  }
  
  # Run bam-readcount in parallel using background jobs
  (
    for chunk_file in $2/bed_chunk_*.bed; do
      run_bam_readcount $chunk_file $2/${base_name}_SiNVICT_inter_bamreadcount/output_$(basename $chunk_file .bed).bamreadcount &
    done
    wait
    echo 'bamreadcount finished'
  )
  
  # Merge bam-readcount output files
  cat $2/${base_name}_SiNVICT_inter_bamreadcount/output_*.bamreadcount > $2/${base_name}_SiNVICT_inter_bamreadcount/output.bamreadcount
  
  # Run SiNVICT with the merged bam-readcount output
  (
    module load HTSlib/1.19.1-GCCcore-11.3.0
    sinvict -m 10 -t $2/${base_name}_SiNVICT_inter_bamreadcount/ -o $2
  )
  
  (
    module load SAMtools/1.19.2-GCCcore-11.3.0
    pythonscript ./sinvict_to_vcf2.py $2/calls_level1.sinvict $2${base_name}_SiNVICT_inter.vcf $Reference
  )

  (
    load_modules

    bcftools sort -o $2${base_name}_SiNVICT_inter_RAW1.vcf -O v $2${base_name}_SiNVICT_inter.vcf
    bcftools filter -O v -i "(FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) >= ${SB} & (FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) <= ${SBU}" $2${base_name}_SiNVICT_inter_RAW1.vcf > $2${base_name}_SiNVICT_RAW.vcf
    bgzip -c $2${base_name}_SiNVICT_RAW.vcf > $2${base_name}_SiNVICT_inter_RAW.vcf.gz

    bcftools index $2${base_name}_SiNVICT_inter_RAW.vcf.gz && cp $2${base_name}_SiNVICT_inter_RAW.vc* $2/mnvcalling  
    
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/VAF > ${VAF} & FORMAT/MDP > ${MRD}" $2${base_name}_SiNVICT_inter_RAW.vcf.gz > $2${base_name}_SiNVICT_inter3.vcf.gz && \
    process_vcf "$2${base_name}_SiNVICT_inter3.vcf.gz" "$2${base_name}_SiNVICT.vcf"  

    rm -rf $2*SiNVICT_inter* && rm -rf $2*alls_level* rm -rf $2*bed_chunk_*
    echo -e "SiNVICT DONE for $base_name\n"
  )
}
################################################### 

generate_mnv_input() {
  (
    load_modules
    
    local mnvpath="$OutputDir/$samplename/mnvcalling"
    
    echo 'MNV input: normalizing, atomizing, zipping and indexing 4 inputs for MNV analysis --> merging into one list'
    bgzip -c $mnvpath/*_M2_F2stats.vcf > $mnvpath/${samplename}_M2_F2stats.vcf.gz && bcftools index $mnvpath/${samplename}_M2_F2stats.vcf.gz 
    
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*VarDict_inter.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/VarDict_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*M2_F2stats.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/Mutect2_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*LoFreq_inter.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/LoFreq_MNV.vcf 2>/dev/null
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $mnvpath/*SiNVICT_inter_RAW.vcf.gz 2>/dev/null | bcftools norm -d all > $mnvpath/SiNVICT_MNV.vcf 2>/dev/null
    
    bgzip -c $mnvpath/VarDict_MNV.vcf > $mnvpath/VarDict_MNV.vcf.gz && bcftools index $mnvpath/VarDict_MNV.vcf.gz
    bgzip -c $mnvpath/Mutect2_MNV.vcf > $mnvpath/Mutect2_MNV.vcf.gz && bcftools index $mnvpath/Mutect2_MNV.vcf.gz
    bgzip -c $mnvpath/LoFreq_MNV.vcf > $mnvpath/LoFreq_MNV.vcf.gz && bcftools index $mnvpath/LoFreq_MNV.vcf.gz
    bgzip -c $mnvpath/SiNVICT_MNV.vcf > $mnvpath/SiNVICT_MNV.vcf.gz && bcftools index $mnvpath/SiNVICT_MNV.vcf.gz
    
    
    #isec to merge into one list --> now with 1 tools. 
    bcftools isec -p $mnvpath -O z -n +1 $mnvpath/SiNVICT_MNV.vcf.gz $mnvpath/Mutect2_MNV.vcf.gz $mnvpath/LoFreq_MNV.vcf.gz $mnvpath/VarDict_MNV.vcf.gz
    #Merge sites with data: --> sites2.txt
    pythonscript ./merge_variants.py $samplename "phased" --infile "sites.txt" --outfile "sites2.txt" $OutputDir
    
    
    ##### MATCHED MODE MNV FILTERING #####
    if [[ "$Matched" == "true" ]]; then
      matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{print $4}' "$CrosslinkTable")
      ctrl_file="$OutputDir/${matched_control}/mnvcalling/sites2.txt"
      curr_file="$mnvpath/sites2.txt"
      #echo -e "MNV input: [$samplename] Matched filtering: ctrl_file = $ctrl_file \t curr_file = $curr_file"
      
      if [[ -f "$ctrl_file" && -f "$curr_file" ]]; then
        total_before=$(wc -l < "$curr_file")
        echo "MNV input: [$samplename] Matched filtering: variants before = $total_before"
        
        pythonscript ./matched_filtering.py $curr_file $ctrl_file $curr_file #ZOIETS?? zou ctrl_file moeten overwriten
        
        total_after=$(wc -l < "$curr_file")
        echo "MNV input: [$samplename] Matched filtering: variants after = $total_after"
      else
        echo "MNV input: [$samplename] Matched filtering skipped -- missing input file(s) or sample is a control"
      fi
    fi
    #### END OF MATCHED MODE MNV FILTERING ####
    
    
    #convert to vcf
    pythonscript ./isec_to_vcf.py $mnvpath/sites2.txt $mnvpath/sites.vcf 'TRUE'
    bgzip -c $mnvpath/sites.vcf > $mnvpath/sites.vcf.gz && bcftools index $mnvpath/sites.vcf.gz
    echo "MNV input: merging 4 tools to 1 list DONE --> sites.vcf"
   
    
#    ############ OPTIONAL MNV FILTERING ############
#    #Filtering in GNOMAD & 1000G HERE: look into opencravat env for .vcf files with snps
#    echo 'MNV input: Filtering: filtering mnv input list on SNP databases'
#    bash ./filter_snps_mnv_input2.sh -a 0.001 -g "$envpath/lib/python3.9/site-packages/cravat/modules/annotators/gnomad3/data/" -b $Listregions -d "$envpath/lib/python3.9/site-packages/cravat/modules/annotators/thousandgenomes/data/"
#    #Filter $OutputDir/$samplename/mnvcalling/sites.txt on 1000G.txt and Gnomad3.txt outputs
#    grep -vF -f <(cut -f1-4 "$envpath/lib/python3.9/site-packages/cravat/modules/annotators/gnomad3/data/gnomad3_$(basename $Listregions)_0.001.txt" "$envpath/lib/python3.9/site-packages/cravat/modules/annotators/thousandgenomes/data/1000g_$(basename $Listregions)_0.001.txt") $mnvpath/sites2.txt > $mnvpath/sites_snp_filtered.txt
#    pythonscript ./isec_to_vcf.py $mnvpath/sites_snp_filtered.txt $mnvpath/sites_snp_filtered.vcf 'TRUE'
#    echo "MNV input: Filtering: filtering on SNP dbs DONE --> sites_snp_filtered.vcf"
#    
#    
#    if [ -z "$PoN" ]; then
#      cp $mnvpath/sites_snp_filtered.txt $mnvpath/merged_vcfs_phased_mnv.txt
#      echo "MNV input: Filtering: filtering on PoN was NOT PERFORMED --> merged_vcfs_phased_mnv.txt (copy from input sites_snp_filtered.txt)"
#    else
#      echo "MNV input: Filtering: filtering mnv input list on PoN..."
#      
#      bgzip -c $mnvpath/sites_snp_filtered.vcf > $mnvpath/sites_snp_filtered.vcf.gz && bcftools index $mnvpath/sites_snp_filtered.vcf.gz #only bgzip and index if pon filtering is performed
#      mkdir $mnvpath/ponfiltered/
#      bcftools isec -p $mnvpath/ponfiltered/ -O v -n~10 $mnvpath/sites_snp_filtered.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
#      bcftools isec -p $mnvpath/ponfiltered2/ -O v -n~10 $mnvpath/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
#      pythonscript ./merge_variants.py $samplename "phased_pon" $OutputDir
#      pythonscript ./isec_to_vcf.py $mnvpath/sites_pon_filtered.txt $mnvpath/sites_pon_filtered.vcf 'TRUE'
#      echo "MNV input: Filtering: filtering on PoN DONE --> merged_vcfs_phased_mnv.txt"
#    fi
#    ############ OPTIONAL MNV FILTERING END ###########
#    
#    # Convert merged_vcfs_phased_mnv.txt into merged_vcfs_phased_mnv.vcf with VAF/RD/MRD metrics
#    echo -e "MNV input: converting mnv input list to vcf...\nMNV input: Generating list for MNV analysis done!"
#    pythonscript ./isec_to_vcf.py $mnvpath/merged_vcfs_phased_mnv.txt $mnvpath/merged_vcfs_phased_mnv.vcf 'TRUE'
#    #bgzip -c $mnvpath/merged_vcfs_phased_mnv.vcf > $mnvpath/merged_vcfs_phased_mnv.vcf.gz && bcftools index $mnvpath/merged_vcfs_phased_mnv.vcf.gz
    
    echo -e "sites.vcf  ==  merged 4 tools MNV input\nsites_snp_filtered.vcf  ==  Filtered on 1000G & Gnomad3 dbs\nmerged_vcfs_phased_mnv.vcf  ==  Subsequently filtered on PoN if available\nsites_pon_filtered.vcf  ==  sites.vcf filtered on ONLY PoN (so without SNP dbs)" > $mnvpath/README2.txt
       
  )
}


mnv_calling() {
  local samplename=$1
  local crosslink_table=$2

  if [ -n "$crosslink_table" ]; then
    echo "MNV input: Crosslink table provided: $crosslink_table"
    
    # Extract the timepoint and corresponding baseline for the sample
    sample_info=$(awk -v sample="$samplename" '$1 == sample {print $2, $3}' "$crosslink_table")
    timepoint=$(echo "$sample_info" | awk '{print $1}')
    corresponding_baseline=$(echo "$sample_info" | awk '{print $2}' | tr -d '\r')
    
    if [ "$timepoint" == "follow-up" ]; then
      # Sample is follow-up
      echo "MNV input: Sample $samplename is a follow-up sample."
      echo -e "MNV input: corresponding baseline sample = $corresponding_baseline"
      echo -e "MNV input: Copying results.vcf from $corresponding_baseline and analyzing that with MNV algorithm..."      

      #fetch list of SNVs (new vcf) that make up the phased MNVs in the baseline sample (so it is pre-filtered already and dont need to look at all variants)
      cp $OutputDir/$corresponding_baseline/mnvcalling/results/results.vcf "$OutputDir/$samplename/mnvcalling/sites.vcf"
      (
        module load HTSlib/1.19.1-GCCcore-11.3.0
        ml GCC/13.3.0
        ./mnv $bam_file "$OutputDir/$samplename/mnvcalling/sites.vcf" $OutputDir/$samplename/mnvcalling/results -T 12 -W 50 -M 6 -S 2 -L 2 -B 0.5 -N 1 -C "./MNV_PoN_2sample_overlap.txt" 
        #As sensitive as possible!! (FU sample)
        echo -e "MNV process for sample $samplename DONE ...\n"
      ) &
    
    elif [ "$timepoint" == "baseline" ]; then
      # Sample is baseline
      echo "MNV input: Sample $samplename is a baseline sample."
      generate_mnv_input
      (
        module load HTSlib/1.19.1-GCCcore-11.3.0
        ml GCC/13.3.0
        ./mnv $bam_file "$OutputDir/$samplename/mnvcalling/sites.vcf" $OutputDir/$samplename/mnvcalling/results -T 12 -W 50 -M 6 -S 15 -L 4 -N 3 -B 0.8 -C "./MNV_PoN_No_overlap_new.txt"
        echo -e "MNV process for sample $samplename DONE ...\n"
      ) &
      
    else
      echo -e "MNV input: Crosslink table is provided but sample $samplename is not annotated as either a 'baseline' or 'follow-up' sample. \nMNV input will be generated for $samplename but not be analyzed"
      generate_mnv_input
    fi
  
  else
    #Treat as baseline sample
    echo "MNV input: Crosslink table not provided. Treating sample as baseline."
    generate_mnv_input
    (
      module load HTSlib/1.19.1-GCCcore-11.3.0
      ml GCC/13.3.0
      #echo -e "\n\nbam_file: $bam_file\nvcf: $vcf\nOutputDir: $OutputDir\nsamplename: $samplename\n\n"
      ./mnv $bam_file "$OutputDir/$samplename/mnvcalling/sites.vcf" $OutputDir/$samplename/mnvcalling/results -T 12 -W 50 -M 6 -S 15 -L 4 -N 3 -B 0.8 -C "./MNV_PoN_No_overlap_new.txt"
      echo -e "MNV process for sample $samplename DONE ...\n"
    ) &
  fi
}



################################################### 
# PON INITIATION: Create PoN if argument is given:
if [ -z "$PoN" ]
then
  #No PoN Argument
  echo -e "\nNo PoN argument given\n"
else
  # If PoN Argument is given:
  echo -e "\nPoN argument is given!!\n"
  mkdir -p ./PoN/
  
  #Directory PoN with containing samples:
  if [[ -d $PoN ]]; then
    echo 'PoN is directory'
    echo -e "\nPoN included: source for PoN .bam files = $PoN"
    for healthycontrol in "$PoN"/*.bam
    do
      samplename=$(basename "$healthycontrol" .bam)
      bam_dir=$(dirname "$healthycontrol")
      bam_index=(${bam_dir}/${samplename}*.bai)
      
      echo -e "\n\nGenerating VCF for normal control: $samplename"
      mkdir ./PoN/$samplename/
      
      # Check if BAM index file exists, if not create it
      if [ ! -e "${bam_index[0]}" ]; then
        echo "BAM index file for ${samplename} not found. Creating BAM index..."
        (
          module load SAMtools/1.16.1-GCCcore-11.3.0
          samtools index "$healthycontrol"
        )
      else
        echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
        touch $bam_index
      fi

      Vardict $healthycontrol ./PoN/$samplename/ & Mutect2 $healthycontrol ./PoN/$samplename/ & Lofreq $healthycontrol ./PoN/$samplename/ & Sinvict $healthycontrol ./PoN/$samplename/
      wait
    done
    
    echo -e "\nMerging outputs for PoN..."
    #Merge per tool -> 4 seperate blacklists as .vcf file
    ls ./PoN/*/*Mutect2.vcf.gz | sed 's/^/ /' > ./PoN/Mutect2-normals.dat & ls ./PoN/*/*VarDict.vcf.gz | sed 's/^/ /' > ./PoN/VarDict-normals.dat & \
    ls ./PoN/*/*LoFreq.vcf.gz | sed 's/^/ /' > ./PoN/LoFreq-normals.dat & ls ./PoN/*/*SiNVICT.vcf.gz | sed 's/^/ /' > ./PoN/SiNVICT-normals.dat
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      xargs -a ./PoN/Mutect2-normals.dat bcftools isec -o ./PoN/PoN_Mutect2.txt -O v -n +1 && xargs -a ./PoN/VarDict-normals.dat bcftools isec -o ./PoN/PoN_Vardict.txt -O v -n +1 && \
      xargs -a ./PoN/LoFreq-normals.dat bcftools isec -o ./PoN/PoN_Lofreq.txt -O v -n +1 && xargs -a ./PoN/SiNVICT-normals.dat bcftools isec -o ./PoN/PoN_Sinvict.txt -O v -n +1
    )
    pythonscript ./isec_to_vcf.py ./PoN/PoN_Mutect2.txt ./PoN/PoN_Mutect2.vcf 'FALSE' && pythonscript ./isec_to_vcf.py ./PoN/PoN_Vardict.txt ./PoN/PoN_Vardict.vcf 'FALSE' && \
    pythonscript ./isec_to_vcf.py ./PoN/PoN_Lofreq.txt ./PoN/PoN_Lofreq.vcf 'FALSE' && pythonscript ./isec_to_vcf.py ./PoN/PoN_Sinvict.txt ./PoN/PoN_Sinvict.vcf 'FALSE'
    
    (
      #zipping and indexing Pon-per-tool blacklists for filtering later on;
      load_modules
      bgzip -c ./PoN/PoN_Mutect2.vcf > ./PoN/PoN_Mutect2.vcf.gz && bgzip -c ./PoN/PoN_Vardict.vcf > ./PoN/PoN_Vardict.vcf.gz && \
      bgzip -c ./PoN/PoN_Lofreq.vcf > ./PoN/PoN_Lofreq.vcf.gz && bgzip -c ./PoN/PoN_Sinvict.vcf > ./PoN/PoN_Sinvict.vcf.gz

      bcftools index ./PoN/PoN_Mutect2.vcf.gz && bcftools index ./PoN/PoN_Vardict.vcf.gz && bcftools index ./PoN/PoN_Lofreq.vcf.gz && bcftools index ./PoN/PoN_Sinvict.vcf.gz
      #We now have 4 PoNs, 1 per tool
      
      #Merging into 1 final blacklist...
      bcftools isec -o ./PoN/BLACKLIST_combined.txt -O v -n +2 ./PoN/PoN_Sinvict.vcf.gz ./PoN/PoN_Mutect2.vcf.gz ./PoN/PoN_Lofreq.vcf.gz ./PoN/PoN_Vardict.vcf.gz
    )
      
      #Add mutations from each tool list (N=4) that occur in over 50% of normal samples, to BLACKLIST_combined.txt 
      # Meaning: Inspect PoN_Vardict.txt etc etc, add high-occurance mutations that are not already present in BLACKLIST_combined.txt to BLACKLIST_combined.txt and re-sort.  
      echo -e "Adding mutations to PoN that are found in 20% of normals, but only called with one tool \n$PoN \n"
      pythonscript ./finish_pon.py 20 ./PoN/

      pythonscript ./isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      load_modules
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
    echo -e "PoN generation is finished! On to the tumor samples... \n"   
  
  #One-File: a.k.a pre-made blacklist:      
  elif [[ -f $PoN ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp $PoN ./PoN/BLACKLIST_combined_complete.txt
    
    pythonscript ./isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      load_modules
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
    echo -e "PoN preparation is finished! On to the tumor samples... \n"   

  #Else (error)
  else
    echo "Invalid PoN input. Input folder with .bam files from healthy controls or a pre-generated PoN.txt "
    exit 1
  fi
fi

#PoN generation finished
###################################################





###################################################
# Function for running a tumor-sample
run_tumor_sample() {
  local bam_file=$1
  local samplename=$(basename "$bam_file" .bam)
  local bam_dir=$(dirname "$bam_file")
  local bam_index=(${bam_dir}/${samplename}*.bai)

  echo -e "\nProcessing sample $samplename\n" 
  mkdir -p $OutputDir/$samplename
  mkdir -p $OutputDir/$samplename/mnvcalling
  
  # Check if BAM index file exists, if not create it
  if [ ! -e "${bam_index[0]}" ]; then
    echo "BAM index file for ${samplename} not found. Creating BAM index..."
    (
      module load SAMtools/1.19.2-GCCcore-11.3.0
      samtools index "$bam_file"
    )
  else
    echo "BAM index file for ${samplename} already exists: ${bam_index[0]}"
    touch $bam_index
  fi
    
  echo "Variant calling in tumor sample..."
  Vardict $bam_file "$OutputDir/$samplename/" &
  Mutect2 $bam_file "$OutputDir/$samplename/" &
  Lofreq $bam_file "$OutputDir/$samplename/" &
  Sinvict $bam_file "$OutputDir/$samplename/" ##################
  wait
  
  ### Merge outputs from all tools, zip and index:
  (
    load_modules
    
    #Combine 4 tools' output into one textfile using intersect
    bcftools isec -p $OutputDir/$samplename/isec -O z -n +$Calls $OutputDir/$samplename/*_SiNVICT.vcf.gz $OutputDir/$samplename/*_Mutect2.vcf.gz $OutputDir/$samplename/*_LoFreq.vcf.gz $OutputDir/$samplename/*_VarDict.vcf.gz
    #Convert to vcf, zip and index
    pythonscript ./isec_to_vcf.py $OutputDir/$samplename/isec/sites.txt $OutputDir/$samplename/isec/sites.vcf 'FALSE'
    bgzip -c $OutputDir/$samplename/isec/sites.vcf > $OutputDir/$samplename/isec/sites.vcf.gz && bcftools index $OutputDir/$samplename/isec/sites.vcf.gz
    #PoN filtering (or not)
    if [ -z "$PoN" ]; then
      echo -e "No PoN argument given, proceeding with intersect results..."
      echo "Adding mutation metrics back to merged mutation list"
      pythonscript ./merge_variants.py $samplename "no_pon" $OutputDir
    else
      #Filter combined list with combined blacklist:
      echo -e "PoN argument given, filtering intersect results with PoN..."
      bcftools isec -p $OutputDir/$samplename/isec/PoN-filtered -O v -n~10 $OutputDir/$samplename/isec/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
      
      # Check if variants are left after PoN Filtering ##########################################      
      
      #Add VAF/RD/MRD data to the merged list:
      echo "Adding mutation metrics back to pon-filtered merged mutation list"
      pythonscript ./merge_variants.py $samplename "pon" $OutputDir
    fi
  )  
  
  #### FILTERING BASED ON MATCHED NORMAL SAMPLE ####
  # --- ADDED: Filter annotated results by matched control in Matched mode ---
  if [[ "$Matched" == "true" ]]; then
    matched_control=$(awk -F $'\t' -v s="$samplename" '$1==s{print $4}' "$CrosslinkTable")
    ctrl_file="$OutputDir/$matched_control/intermediate-files/merged_vcfs.txt"
    curr_file="$OutputDir/$samplename/merged_vcfs.txt"
    echo -e "ctrl_file = $ctrl_file \ncurr_file = $curr_file"
    if [[ -f "$ctrl_file" && -f "$curr_file" ]]; then
      echo "[$samplename] Matched filtering: variants before = $(wc -l < ${curr_file})"
      pythonscript ./matched_filtering.py $curr_file $ctrl_file $curr_file
      echo "[$samplename] Matched filtering: variants after = $(wc -l < ${curr_file})"
    else
      echo "[$samplename] Matched filtering skipped -- missing input file(s) or sample is a control"
    fi
  fi
  # -----------------------------------------------------------------------
  #### END OF FILTERING BASED ON MATCHED NORMAL SAMPLE ####
  
  #### Filtering on UNMET score ###
  if [[ "$UNMET" == "true" ]]; then
  BED_UNMET="./hg38_genome.UNMET_bin10bp.bed"
  curr_file="$OutputDir/$samplename/merged_vcfs.txt"
    #Create UNMET .bed file if needed
    if [[ ! -f "$BED_UNMET" ]]; then
      echo "UNMET is true and $BED_UNMET not yet created. Creating $BED_UNMET..."
      wget http://web16.kazusa.or.jp/data/hg38/hg38_genome.UNMET_bin10bp.tdf
      (
        module load IGV
        igvtools tdftobedgraph ./hg38_genome.UNMET_bin10bp.tdf ./hg38_genome.UNMET_bin10bp.bed 2>/dev/null
      )
      echo "UNMET: $BED_UNMET created"
    else
      echo "UNMET is true and $BED_UNMET already exists. "
    fi
    # Attaching UNMET scores to merged_vcfs.txt
    echo "[$samplename] UNMET filtering: variants before = $(wc -l < ${curr_file})"
    (
      module load BEDTools
      bedtools intersect -a <(
        awk 'BEGIN{FS=OFS="\t"} { print $1, $2, $2+1, $0 }' "$curr_file"
      ) \
        -b "$BED_UNMET" -loj \
      | awk 'BEGIN{FS=OFS="\t"} {
          sc = ($12 == "." ? "NA" : $12)
          printf "%s\t%s\t%s\t%s\t%s_UNMET:%s\n",
            $4, $5, $6, $7, $8, sc
        }' \
      > "${curr_file%.txt}_UNMET.txt"
    )
    # Filtering on UNMET > 0.93
    awk -F'\t' '
      {
        split($5, a, "_UNMET:")
        sc = a[2]
        if (sc == "NA" || sc+0 < 0.93) print
      }
    ' "${curr_file%.txt}_UNMET.txt" > "$curr_file"
    echo "[$samplename] UNMET filtering: variants after = $(wc -l < ${curr_file})"
  fi
    
  ### End of filtering on UNMET score ###
    
  #MNVCalling
  if [ "$MNV" == "true" ]; then
    #If phased argument is true: 
    echo -e "\nMNV process start: -H argument is set to 'true', proceeding with creating input list for MNV calling"
    mnv_calling $samplename $CrosslinkTable
  elif [ "$MNV" == "false" ]; then
    echo "-H argument is set to 'false'; not proceeding with MNV calling"
    #rm -rf $OutputDir/$samplename/mnvcalling/
  else
    echo "-H argument not supplied or is not 'true' nor 'false'; not proceeding with MNV calling"
    #rm -rf $OutputDir/$samplename/mnvcalling/
  fi
  
  #1st annotation step
  echo "First step of annotating mutations..."
  awk -v samplename="$samplename" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" samplename "\t" $5}' "$OutputDir/$samplename/merged_vcfs.txt" > $OutputDir/$samplename/merged_vcfs2.txt
  annotate $OutputDir/$samplename/merged_vcfs2.txt "merged_vcfs_annotated" excel
    
  #Filtering annotated output
  echo "Post filtering based on functional effect and public SNP databases..." 
  pythonscript ./post_filtering.py $samplename $OutputDir $UNMET
    
  #Annotating again if there are mutations left:
  MUTS=$(awk '{ FS="\t" } END { print $2 }' $OutputDir/$samplename/filtering_info.txt)
  echo "$MUTS mutations left after filtering..."
  if [[ $MUTS -eq 0 ]]; then #If 0 mutations left:
    echo 'NO mutations left after filtering...' 
    echo "After post-filtering, no mutations were left to annotate or plot, therefore, these steps are skipped for this sample" > $OutputDir/$samplename/NO_MUTS.txt
  else #If mutations are found:
    echo 'At least one mutation remains, on to annotating final list...' 
    annotate $OutputDir/$samplename/merged_vcfs_annotated-filtered.txt "Final_list_$samplename" excel vcf
  fi
  
  #Wait for MNV algorithm to finish
  wait
  
  #Cleaning:
  mkdir $OutputDir/$samplename/raw-data/ && mkdir $OutputDir/$samplename/filtering-info/ && mkdir $OutputDir/$samplename/intermediate-files/ && mkdir $OutputDir/$samplename/final-list/
  mv $OutputDir/$samplename/filtering_info.txt $OutputDir/$samplename/filtering-info/ && mv $OutputDir/$samplename/removed_variants.xlsx $OutputDir/$samplename/filtering-info/
  mv $OutputDir/$samplename/merged* $OutputDir/$samplename/intermediate-files/ && mv $OutputDir/$samplename/Final_list* $OutputDir/$samplename/final-list/
  
  echo -e "DONE analyzing $samplename!!\n\n\n" 
}

# Function for tumor-sample finished
###################################################



###################################################
# Running Tumor samples from directory or as single file
mkdir -p $OutputDir
if [[ -d $Inputbam ]]; then
  if [[ -n "$CrosslinkTable" && -f "$CrosslinkTable" ]]; then
    if [[ "$Matched" == "true" ]]; then
      echo "Matched mode: controls --> baselines --> follow-ups"

      # 1) Controls
      awk -F $'\t' 'NR>1 && $2=="baseline"{print $4}' "$CrosslinkTable" | sort -u | \
      while read -r ctrl; do
        run_tumor_sample "$Inputbam/${ctrl}.bam"
      done

      # 2) Baselines
      awk -F $'\t' 'NR>1 && $2=="baseline"{print $1}' "$CrosslinkTable" | \
      while read -r samp; do
        run_tumor_sample "$Inputbam/${samp}.bam"
      done

      # 3) Follow-ups
      awk -F $'\t' 'NR>1 && $2=="follow-up"{print $1}' "$CrosslinkTable" | \
      while read -r samp; do
        run_tumor_sample "$Inputbam/${samp}.bam"
      done
    else
      # Non-matched: baselines then follow-ups
      echo "Standard order: baselines --> follow-ups"
      samples=( $(awk -F $'\t' 'NR>1{if($2=="baseline")print$1}NR>1{if($2=="follow-up")print$1}' "$CrosslinkTable") )
      echo -e "\n\n\nSamples order = \n${samples[@]}\n\n\n"
      for s in "${samples[@]}"; do
        run_tumor_sample "$Inputbam/${s}.bam"
      done
    fi
  else
    # No table: all BAMs
    echo "No CrosslinkTable: processing all BAMs"
    for bam in $Inputbam/*.bam; do
      run_tumor_sample "$bam"
    done
  fi
elif [[ -f $Inputbam ]]; then
  run_tumor_sample "$Inputbam"
else
  echo "Invalid Inputbam: $Inputbam"
  exit 1
fi

# Running Tumor samples finished
###################################################

echo 'Finished run!  '











