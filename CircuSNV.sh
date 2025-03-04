#!/bin/bash
#SBATCH --job-name=CircuSNV
#SBATCH --output="CircuSNV_run.out"
#SBATCH --error="CircuSNV_run.err"
#SBATCH --time=100:00:00
#SBATCH --mem=12gb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --get-user-env=L60
#SBATCH --export=NONE

while getopts "R:L:I:O:V:D:C:P:Q:B:M:F:S:N:" arg; do 
  case $arg in
    R) Reference=$OPTARG;;      # "/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly38.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation_hg38.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
    Q) Qual=$OPTARG;;           # BaseQuality (18)
    M) MRD=$OPTARG;;            # Minimum reads supporting a mutation
    F) FilterMode=$OPTARG;;     # Pon Filtering mode: "PerTool" or "Combined"
    S) SB=$OPTARG;;             # Strand bias. Lowest fraction of FW reads vs RV reads. (SB = 0.2 --> FW = 1, RV = 5 and vice versa)
    N) MNVlength=$OPTARG;;      # Maximum MNV length. Determines how many bp's SNVs can be apart before merging into MNVs. If set to 1, it will only output SNvs, so no complex variants (apart from InDels)
  esac
done

ref2bit="/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly38.2bit" #Change if needed
SBU=$(echo "scale=2; 1 / $SB" | bc)

cd /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/CircuSNV/

echo -e "\nSettings: \n\nInput bam(s): $Inputbam \nReference: $Reference \nPanel: $Listregions\nminimum VAF: $VAF\nminimum Read Depth: $RDP\nminimum BaseQuality: $Qual\nminimum overlapping calls: $Calls\nPoN:$PoN"

##### Python and annotate function: ---Change env directory eventually...---
pythonscript() {
  source /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/activate
  python3 $1 $2 $3 $4 $5 $6 $7
  deactivate
}

annotate() { #1 input.txt (tab delim snv data)   #2 Output filenames    #3&#4 output 
  echo -e "annotating $1 with OpenCravat..."
  source /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/activate
  oc run $1 \
    -l hg38 \
    -n $2 --silent --mp 8 \
    -a clinvar cgc cadd chasmplus clinpred cosmic dbsnp mutation_assessor\
        thousandgenomes vest gnomad3 \
    -t $3 $4
  deactivate 
}

################################################### 
# Functions to run each tool, filter on Read Dept, VAF and Mutant allele dept & creating .vcf file

Vardict() { #1: Input .bam file   #2 Output folder
  echo "Analyzing $(basename $1 .bam) with VarDict..."
    (
      module load R/4.2.2-foss-2022a-bare   
    
      vardict-java -th 8 -r $MRD -G $Reference -f $VAF -N $(basename $1 .bam) -X $MNVlength -I 20 -L 20 -Y 20 --nosv --deldupvar -q $Qual -Q $Qual -b $1 -c 1 -S 2 -E 3 -g 4 $Listregions | \
      teststrandbias.R | var2vcf_valid.pl -S -N $(basename $1 .bam) -E -f $VAF > $2$(basename $1 .bam)_VarDict_RAW.vcf
    )
    
    (
      module load zlib/1.2.12-GCCcore-11.3.0
      module load HTSlib/1.19.1-GCCcore-11.3.0
      module load BCFtools/1.19-GCCcore-11.3.0
      
      bgzip -c $2$(basename $1 .bam)_VarDict_RAW.vcf > $2$(basename $1 .bam)_VarDict_inter.vcf.gz && bcftools index $2$(basename $1 .bam)_VarDict_inter.vcf.gz && \
      bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/VD > ${MRD} & INFO/PMEAN < 67 & (FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) >= ${SB} & (FORMAT/ALD[0:0]/FORMAT/ALD[0:1]) <= ${SBU}" $2$(basename $1 .bam)_VarDict_inter.vcf.gz > $2$(basename $1 .bam)_VarDict_inter2.vcf.gz && \
      bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $2$(basename $1 .bam)_VarDict_inter2.vcf.gz | bcftools norm -d all > $2$(basename $1 .bam)_VarDict.vcf && \
      bgzip $2$(basename $1 .bam)_VarDict.vcf && bcftools index $2$(basename $1 .bam)_VarDict.vcf.gz && \
      rm -rf $2$(basename $1 .bam)_VarDict_inter*.*
      echo -e 'VarDict DONE\n'
    )
}

######
Mutect2() { #1: Input .bam file    #2: Output folder
  echo "Analyzing $(basename $1 .bam) with Mutect2..."
  (
    module load GATK/4.2.4.1-Java-8-LTS
    module load BCFtools/1.19-GCCcore-11.3.0
    module load zlib/1.2.12-GCCcore-11.3.0
    module load HTSlib/1.19.1-GCCcore-11.3.0
    
    gatk --java-options -Xmx8g Mutect2 --QUIET --verbosity ERROR --create-output-variant-index false --native-pair-hmm-threads 8 --max-mnp-distance $MNVlength \
      --min-base-quality-score $Qual --callable-depth $RDP --minimum-allele-fraction $VAF -R $Reference -I $1 -O $2$(basename $1 .bam)_Mutect2_RAW.vcf -L $Listregions $3 && \
    gatk FilterMutectCalls --QUIET --verbosity ERROR --create-output-variant-index false \
      --min-allele-fraction $VAF -R $Reference -V $2$(basename $1 .bam)_Mutect2_RAW.vcf -O $2$(basename $1 .bam)_M2_Fstats.vcf && \
    grep -v '^##contig=<ID=chrUn\|^##contig=<ID=HLA' $2$(basename $1 .bam)_M2_Fstats.vcf | grep -v '_alt' | grep -v '_random' > $2$(basename $1 .bam)_M2_F2stats.vcf && \
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/AD[0:1] > ${MRD} & (FORMAT/SB[0:2]/FORMAT/SB[0:3]) >= ${SB} & (FORMAT/SB[0:2]/FORMAT/SB[0:3]) <= ${SBU}" $2$(basename $1 .bam)_M2_F2stats.vcf > $2$(basename $1 .bam)M2_F3stats.vcf.gz && \
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $2$(basename $1 .bam)M2_F3stats.vcf.gz | bcftools norm -d all > $2$(basename $1 .bam)_Mutect2.vcf && \
    bgzip $2$(basename $1 .bam)_Mutect2.vcf && bcftools index $2$(basename $1 .bam)_Mutect2.vcf.gz && \
    rm -rf $2$(basename $1 .bam)*tats*
    echo -e 'Mutect2 DONE\n'
  )
}


######
Lofreq() { #1: Input .bam file    #2: Output folder
  echo "Analyzing $(basename $1 .bam) with LoFreq..."
  (
    module load BCFtools/1.19-GCCcore-11.3.0
    module load zlib/1.2.12-GCCcore-11.3.0
    module load HTSlib/1.19.1-GCCcore-11.3.0
    
    lofreq call-parallel --pp-threads 8 --no-default-filter --call-indels -C $RDP -f $Reference -q $Qual -Q $Qual -o $2$(basename $1 .bam)_LoFreq_inter.vcf -l $Listregions $1 && \
    awk '!seen[$0]++' $2$(basename $1 .bam)_LoFreq_inter.vcf > $2$(basename $1 .bam)_LoFreq_RAW.vcf && \

    bgzip -c $2$(basename $1 .bam)_LoFreq_RAW.vcf > $2$(basename $1 .bam)_LoFreq_inter.vcf.gz && bcftools index $2$(basename $1 .bam)_LoFreq_inter.vcf.gz
    bcftools filter -O z -i "INFO/DP > ${RDP} & INFO/AF > ${VAF} & INFO/DP4[2]+INFO/DP4[3] > ${MRD} & (INFO/DP4[2]/INFO/DP4[3]) >= ${SB} & (INFO/DP4[2]/INFO/DP4[3]) <= ${SBU}" $2$(basename $1 .bam)_LoFreq_inter.vcf.gz > $2$(basename $1 .bam)_LoFreq_inter2.vcf.gz && \
    bcftools norm -O v -m -any -f $Reference -a --atom-overlaps . -c w $2$(basename $1 .bam)_LoFreq_inter2.vcf.gz | bcftools norm -d all > $2$(basename $1 .bam)_LoFreq.vcf && \
    bgzip $2$(basename $1 .bam)_LoFreq.vcf && bcftools index $2$(basename $1 .bam)_LoFreq.vcf.gz && \
    rm -rf $2$(basename $1 .bam)_LoFreq_inter*.*
    echo -e 'LoFreq DONE\n'
  )
}

######
Sinvict() {
  InputSinvict=$1
  #1: Input .bam file    #2: Output folder
  echo "Analyzing $(basename $1 .bam) with bam-readcount && SiNVICT..."
  
  # Create output directories
  mkdir -p $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/
  
  # Define the number of chunks for parallel processing
  N=8  # Adjust this value based on your system's capabilities
  
  # Split the BED file into N chunks (excluding the header lines)
  header_lines=$(grep -P "^browser|^track" $Listregions | wc -l)
  tail -n +$((header_lines + 1)) $Listregions | split -l $((($(wc -l < $Listregions) - header_lines) / N + 1)) -d -a 2 --additional-suffix=.bed - $2/bed_chunk_
  
  # Add the header lines back to each chunk
  for chunk_file in $2/bed_chunk_*.bed; do
    (grep -P "^browser|^track" $Listregions; cat $chunk_file) > ${chunk_file}.tmp && mv ${chunk_file}.tmp $chunk_file
  done
  
  # Function to run bam-readcount for a chunk
  run_bam_readcount() {
    local chunk_file=$1
    local output_file=$2
    module load zlib/1.2.12-GCCcore-11.3.0
    module load HTSlib/1.19.1-GCCcore-11.3.0
    bam-readcount -l $chunk_file -w 1 -b $Qual -f $Reference $InputSinvict > $output_file
  }
  
  # Run bam-readcount in parallel using background jobs
  (
    for chunk_file in $2/bed_chunk_*.bed; do
      run_bam_readcount $chunk_file $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/output_$(basename $chunk_file .bed).bamreadcount &
    done
    
    # Wait for all background jobs to finish
    wait
  )
  
  # Merge bam-readcount output files
  cat $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/output_*.bamreadcount > $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/output.bamreadcount
  
  # Run SiNVICT with the merged bam-readcount output
  (
    module load HTSlib/1.19.1-GCCcore-11.3.0
    sinvict -m $RDP -t $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/ -o $2
  )
  
  (
    module load SAMtools/1.16.1-GCCcore-11.3.0
    pythonscript ./sinvict_to_vcf2.py $2/calls_level1.sinvict $2$(basename $1 .bam)_SiNVICT_inter.vcf $Reference
  )

  (
    module load BCFtools/1.19-GCCcore-11.3.0
    module load zlib/1.2.12-GCCcore-11.3.0
    module load HTSlib/1.19.1-GCCcore-11.3.0

    bcftools sort -o $2$(basename $1 .bam)_SiNVICT_RAW.vcf -O v $2$(basename $1 .bam)_SiNVICT_inter.vcf
    bgzip -c $2$(basename $1 .bam)_SiNVICT_RAW.vcf > $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz
    bcftools index $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz && \
    bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/VAF > ${VAF} & FORMAT/MDP > ${MRD} & (FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) >= ${SB} & (FORMAT/MSB[0:0]/FORMAT/MSB[0:1]) <= ${SBU}" $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz > $2$(basename $1 .bam)_SiNVICT_inter3.vcf.gz && \
    bcftools norm -O v -m-any -f $Reference -a --atom-overlaps . -c w $2$(basename $1 .bam)_SiNVICT_inter3.vcf.gz | bcftools norm -d all > $2$(basename $1 .bam)_SiNVICT.vcf && \
    bgzip $2$(basename $1 .bam)_SiNVICT.vcf && bcftools index $2$(basename $1 .bam)_SiNVICT.vcf.gz && \
    rm -rf $2*SiNVICT_inter* && rm -rf $2*alls_level*
    echo -e 'SiNVICT DONE\n'
  )
}
################################################### 



################################################### 
# PON INITIATION: Create PoN if argument is given:
if [ -z "$PoN" ]
then
  #No PoN Argument
  echo -e "\nNo PoN argument given\n"
else
  # If PoN Argument is given:
  echo -e "\nPoN argument is given!!\n"
  mkdir ./PoN/
  
  #Directory PoN with containing samples:
  if [[ -d $PoN ]]; then
    echo 'PoN is directory'
    echo -e "\nPoN included: source for PoN .bam files = $PoN"
    for healthycontrol in "$PoN"/*.bam
    do
      echo -e '\n\nGenerating VCF for normal control: '
      echo -e normalcontrol= $(basename $healthycontrol .bam) '\n'
      (
        module load SAMtools/1.16.1-GCCcore-11.3.0
        samtools index $healthycontrol
      )
      Vardict $healthycontrol ./PoN/ & Mutect2 $healthycontrol ./PoN/ & Lofreq $healthycontrol ./PoN/ & Sinvict $healthycontrol ./PoN/
    done
    
    echo -e "\nMerging outputs for PoN..."
    #Merge per tool -> 4 seperate blacklists as .vcf file
    ls ./PoN/*Mutect2.vcf.gz | sed 's/^/ /' > ./PoN/Mutect2-normals.dat & ls ./PoN/*VarDict.vcf.gz | sed 's/^/ /' > ./PoN/VarDict-normals.dat & \
    ls ./PoN/*LoFreq.vcf.gz | sed 's/^/ /' > ./PoN/LoFreq-normals.dat & ls ./PoN/*SiNVICT.vcf.gz | sed 's/^/ /' > ./PoN/SiNVICT-normals.dat
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      xargs -a ./PoN/Mutect2-normals.dat bcftools isec -o ./PoN/PoN_Mutect2.txt -O v -n +1 && xargs -a ./PoN/VarDict-normals.dat bcftools isec -o ./PoN/PoN_Vardict.txt -O v -n +1 && \
      xargs -a ./PoN/LoFreq-normals.dat bcftools isec -o ./PoN/PoN_Lofreq.txt -O v -n +1 && xargs -a ./PoN/SiNVICT-normals.dat bcftools isec -o ./PoN/PoN_Sinvict.txt -O v -n +1
    )
    pythonscript ./isec_to_vcf.py ./PoN/PoN_Mutect2.txt ./PoN/PoN_Mutect2.vcf 'FALSE' && pythonscript ./isec_to_vcf.py ./PoN/PoN_Vardict.txt ./PoN/PoN_Vardict.vcf 'FALSE' && \
    pythonscript ./isec_to_vcf.py ./PoN/PoN_Lofreq.txt ./PoN/PoN_Lofreq.vcf 'FALSE' && pythonscript ./isec_to_vcf.py ./PoN/PoN_Sinvict.txt ./PoN/PoN_Sinvict.vcf 'FALSE'
    
    (
      #zipping and indexing Pon-per-tool blacklists for filtering later on;
      module load BCFtools/1.19-GCCcore-11.3.0
      module load HTSlib/1.19.1-GCCcore-11.3.0
      bgzip -c ./PoN/PoN_Mutect2.vcf > ./PoN/PoN_Mutect2.vcf.gz && bgzip -c ./PoN/PoN_Vardict.vcf > ./PoN/PoN_Vardict.vcf.gz && \
      bgzip -c ./PoN/PoN_Lofreq.vcf > ./PoN/PoN_Lofreq.vcf.gz && bgzip -c ./PoN/PoN_Sinvict.vcf > ./PoN/PoN_Sinvict.vcf.gz

      bcftools index ./PoN/PoN_Mutect2.vcf.gz && bcftools index ./PoN/PoN_Vardict.vcf.gz && bcftools index ./PoN/PoN_Lofreq.vcf.gz && bcftools index ./PoN/PoN_Sinvict.vcf.gz
      #We now have 4 PoNs, 1 per tool
      
      #Merging into 1 final blacklist...
      bcftools isec -o ./PoN/BLACKLIST_combined.txt -O v -n +2 ./PoN/PoN_Sinvict.vcf.gz ./PoN/PoN_Mutect2.vcf.gz ./PoN/PoN_Lofreq.vcf.gz ./PoN/PoN_Vardict.vcf.gz
    )
      
      #Add mutations from each tool list (N=4) that occur in over 50% of normal samples, to BLACKLIST_combined.txt 
      # Meaning: Inspect PoN_Vardict.txt etc etc, add high-occurance mutations that are not already present in BLACKLIST_combined.txt to BLACKLIST_combined.txt and re-sort.  
      echo -e "Adding mutations to PoN that are found in 50% of normals, but only called with one tool \n$PoN \n"
      pythonscript ./finish_pon.py 10 ./PoN/


      pythonscript ./isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      module load HTSlib/1.19.1-GCCcore-11.3.0
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
      
  
  #One-File: a.k.a pre-made blacklist:      
  elif [[ -f $PoN ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp $PoN ./PoN/BLACKLIST_combined_complete.txt
    
    pythonscript ./isec_to_vcf.py ./PoN/BLACKLIST_combined_complete.txt ./PoN/BLACKLIST_combined.vcf 'FALSE' && \
      
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      module load HTSlib/1.19.1-GCCcore-11.3.0
      bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz
      bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
    ) 
      
    
  
  #Else (error)
  else
    echo "Invalid PoN input. Input folder with .bam files from healthy controls or a pre-generated PoN.txt "
    exit 1
  fi
fi
echo -e "\nPoN generation is finished! On to the tumor samples... \n" 
#PoN generation finished
###################################################



###################################################
# Function for running a tumor-sample
run_tumor_sample() {
  samplename=$(basename $1 .bam)
  echo -e "\n\n\nProcessing sample $samplename\n" 
  mkdir ./output/$samplename
  echo "Variant calling in tumor sample..."
  (
    module load SAMtools/1.16.1-GCCcore-11.3.0
    samtools index $1
  )
  Vardict $1 "./output/$samplename/" &
  Mutect2 $1 "./output/$samplename/" &
  Lofreq $1 "./output/$samplename/" &
  Sinvict $1 "./output/$samplename/" ##################
  wait
  
  # If -P is empty:
    # Don't filter on PoN, prepare list for MNV handling directly (adjust merge_variants):
    #- isec variants (similar to Combined mode, but without the PoN filtering)
    #merge_variants.py
    
  # elif pon argument exists: do below 
  
  #Merging, Annotating & Post-filtering:
  if [ "$FilterMode" = "PerTool" ]; then
    ###  PerTool PoN Mode ###
    echo -e "\n\nPoN Filtering mode is set to 'PerTool'"
    # Deprecated
    
  elif [ "$FilterMode" = "Combined" ]; then
    ### Combined PoN mode: ### 
    echo -e "\n\nPoN Filtering mode is set to 'Combined'"
    #Merge outputs from all tools, zip and index:
    (
      module load BCFtools/1.19-GCCcore-11.3.0
      module load zlib/1.2.12-GCCcore-11.3.0
      module load HTSlib/1.19.1-GCCcore-11.3.0
      
      bcftools isec -p ./output/$samplename/Combined_PoN-filtered -O z -n +$Calls ./output/$samplename/*_SiNVICT.vcf.gz ./output/$samplename/*_Mutect2.vcf.gz ./output/$samplename/*_LoFreq.vcf.gz ./output/$samplename/*_VarDict.vcf.gz
      pythonscript ./isec_to_vcf.py ./output/$samplename/Combined_PoN-filtered/sites.txt ./output/$samplename/Combined_PoN-filtered/sites.vcf 'FALSE'
      bgzip -c ./output/$samplename/Combined_PoN-filtered/sites.vcf > ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz
      bcftools index ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz
      #Filter combined list with combined blacklist:
      bcftools isec -p ./output/$samplename/Combined_PoN-filtered/merged_PoN-filtered -O v -n~10 ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz

    )  
    
    ############################################ Check if variants are left after PoN Filtering ########################################## 
    #Add VAF/RD/MRD data to the merged list:
    pythonscript ./merge_variants.py $samplename "Combined"
    # The same, but for the non-PoN filtered mutations for phased variant calling. 
    pythonscript ./merge_variants.py $samplename "phased"
    
    
    #
  else
    echo "$FilterMode is not a valid option for FilterMode... please provide either 'PerTool' or 'Combined'... "
    exit 1
  fi
  
  
  # Merging SNVs into MNVs if provided
  if [ -z "$MNVlength" ] || [ "$MNVlength" -eq 1 ]; then
    echo "MNVlength is not provided or equal to 1. Skipping joining SNVs."
    awk -v samplename="$samplename" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" samplename "\t" $5}' "./output/$samplename/merged_vcfs.txt" > ./output/$samplename/merged_vcfs2.txt
    annotate ./output/$samplename/merged_vcfs2.txt "merged_vcfs_annotated" excel
    
  else
    echo "Combining SNVs that should be MNVs in a range of $MNVlength bp..."
    pythonscript ./isec_to_vcf.py ./output/$samplename/merged_vcfs.txt ./output/$samplename/merged_vcfs.vcf 'TRUE'
    (
        module load SAMtools/1.16.1-GCCcore-11.3.0
        pythonscript ./joinAdjacentSNPs_test.py $MNVlength $1 $ref2bit -v "./output/$samplename/merged_vcfs.vcf" -o="./output/$samplename/merged_vcfs_MNV.vcf"
    )
    grep -v '^#' ./output/$samplename/merged_vcfs_MNV.vcf | awk -v samplename="$samplename" '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" samplename "\t" $8}' > ./output/$samplename/merged_vcfs_MNV.txt
    annotate ./output/$samplename/merged_vcfs_MNV.txt "merged_vcfs_annotated" excel
  fi
  
  #Filtering annotated output
  echo "Post filtering..." 
  pythonscript ./post_filtering.py $samplename
  
  #Annotating again if there are mutations left:
  MUTS=$(awk '{ FS="\t" } END { print $2 }' ./output/$samplename/filtering_info.txt)
  echo "$MUTS mutations left after filtering..."
  if [[ $MUTS -eq 0 ]]; then #If 0 mutations left:
    echo 'NO mutations left after filtering...' 
    echo "After post-filtering, no mutations were left to annotate or plot, therefore, these steps are skipped for this sample" > ./output/$samplename/NO_MUTS.txt
  else #If mutations are found:
    echo 'At least one mutation remains, on to annotating final list...' 
    annotate ./output/$samplename/merged_vcfs_annotated-filtered.txt "Final_list_$samplename" excel vcf
  fi
  
  #Cleaning:
  mkdir ./output/$samplename/raw-data/ && mkdir ./output/$samplename/filtering-info/ && mkdir ./output/$samplename/intermediate-files/ && mkdir ./output/$samplename/final-list/
  mv ./output/$samplename/filtering_info.txt ./output/$samplename/filtering-info/ && mv ./output/$samplename/removed_variants.xlsx ./output/$samplename/filtering-info/
  mv ./output/$samplename/merged* ./output/$samplename/intermediate-files/ && mv ./output/$samplename/Final_list* ./output/$samplename/final-list/
  mv ./output/$samplename/*-filtered/ ./output/$samplename/raw-data/
  mv ./output/$samplename/*.log ./output/$samplename/filtering-info
  
  echo -e "DONE analyzing $samplename!!\n\n\n" 
}
# Function for tumor-sample finished
###################################################



###################################################
# Running Tumor samples from directory or as single file
mkdir ./output
#Directory:
if [[ -d $Inputbam ]]; then
  for tumorsample in "$Inputbam"/*.bam
  do
    run_tumor_sample $tumorsample
  done
  #Batch analytics here
#One-File:     
elif [[ -f $Inputbam ]]; then
  run_tumor_sample $Inputbam
#Else (error)
else
  echo "$Inputbam is not valid"
  exit 1
fi
# Running Tumor samples finished
###################################################

echo 'Finished run!  '
