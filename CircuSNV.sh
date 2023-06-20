#!/bin/bash

while getopts "R:L:I:O:V:D:C:P:Q:B:M:F:" arg; do 
  case $arg in
    R) Reference=$OPTARG;;      # "/groups/umcg-pmb/tmp01/apps/data/reference_sequences/Homo_sapiens_assembly38.fasta"
    L) Listregions=$OPTARG;;    # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/2020_covered_probe_notranslocation_hg38.bed"
    I) Inputbam=$OPTARG;;       # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/103584/"
    V) VAF=$OPTARG;;            # 0.004
    D) RDP=$OPTARG;;            # 100
    C) Calls=$OPTARG;;          # 1
    P) PoN=$OPTARG;;            # "/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/GenomeScan_SequenceData/104103_normals/"
    Q) Qual=$OPTARG;;           # BaseQuality (18)
    M) MRD=$OPTARG;;
    F) FilterMode=$OPTARG;;     # Pon Filtering mode: "PerTool" or "Combined"
  esac
done

cd /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/CircuSNV/ ### Change this directory ###

echo -e "\nSettings: \n\nInput bam(s): $Inputbam \nReference: $Reference \nPanel: $Listregions\nminimum VAF: $VAF\nminimum Read Depth: $RDP\nminimum BaseQuality: $Qual\nminimum overlapping calls: $Calls\nPoN:$PoN"
echo -e '\nLoading modules: \n'

module load R/4.0.3-foss-2018b-bare
module load GATK/4.2.4.1-Java-8-LTS
module load SAMtools/1.16.1-GCCcore-11.3.0
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 #??
module load HTSlib/1.11-GCCcore-7.3.0


##### Python and annotate function: ---Change env directory eventually...---
pythonscript() {
  source /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/activate ### Change this directory ###
  python3 $1 $2 $3 $4 $5 $6 $7
  deactivate
}

annotate() {
  echo -e "annotating $1 with OpenCravat..."
  source /groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/activate  ### Change this directory ###
  oc run $1 \
    -l hg38 \
    -n $2 --silent \
    -a clinvar cgc cadd chasmplus clinpred cosmic dbsnp mutation_assessor\
        thousandgenomes thousandgenomes_european vest cadd_exome gnomad3 \
    -t $3 $4
  deactivate 
}

################################################### 
# Functions to run each tool, filter on Read Dept, VAF and Mutant allele dept & creating .vcf file
######
Vardict() { #1: Input .bam file   #2 Output folder
echo "Analyzing $(basename $1 .bam) with VarDict..."
vardict-java -th 14 -r $MRD -G $Reference -f $VAF -N $(basename $1 .bam) -X 20 -q $Qual -b $1 -c 1 -S 2 -E 3 -g 4 $Listregions | \
teststrandbias.R | var2vcf_valid.pl -S -N $(basename $1 .bam) -E -f $VAF > $2$(basename $1 .bam)_VarDict_RAW.vcf && \
bgzip -c $2$(basename $1 .bam)_VarDict_RAW.vcf > $2$(basename $1 .bam)_VarDict_inter.vcf.gz && bcftools index $2$(basename $1 .bam)_VarDict_inter.vcf.gz && \
bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/VD > ${MRD} & INFO/PMEAN < 67" $2$(basename $1 .bam)_VarDict_inter.vcf.gz > $2$(basename $1 .bam)_VarDict.vcf.gz && \
bcftools index $2$(basename $1 .bam)_VarDict.vcf.gz && \
rm -rf $2$(basename $1 .bam)_VarDict_inter.*
}

######
Mutect2() { #1: Input .bam file    #2: Output folder
echo "Analyzing $(basename $1 .bam) with Mutect2..."
gatk Mutect2 --QUIET --verbosity ERROR --create-output-variant-index false --native-pair-hmm-threads 14 --max-mnp-distance 20 \
  --min-base-quality-score $Qual --callable-depth $RDP --minimum-allele-fraction $VAF -R $Reference -I $1 -O $2$(basename $1 .bam)_Mutect2_RAW.vcf -L $Listregions $3 && \
gatk FilterMutectCalls --QUIET --verbosity ERROR --create-output-variant-index false \
  --min-allele-fraction $VAF -R $Reference -V $2$(basename $1 .bam)_Mutect2_RAW.vcf -O $2$(basename $1 .bam)_M2_Fstats.vcf && \
grep -v '^##contig=<ID=chrUn\|^##contig=<ID=HLA' $2$(basename $1 .bam)_M2_Fstats.vcf | grep -v '_alt' | grep -v '_random' > $2$(basename $1 .bam)_M2_F2stats.vcf && \
bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/AF > ${VAF} & FORMAT/AD[0:1] > ${MRD}" $2$(basename $1 .bam)_M2_F2stats.vcf > $2$(basename $1 .bam)M2_F3stats.vcf.gz && \

bcftools norm -m-both -O z $2$(basename $1 .bam)M2_F3stats.vcf.gz > $2$(basename $1 .bam)_Mutect2.vcf.gz && \

bcftools index $2$(basename $1 .bam)_Mutect2.vcf.gz && \
rm -rf $2$(basename $1 .bam)*tats*
}

######
Lofreq() { #1: Input .bam file    #2: Output folder
echo "Analyzing $(basename $1 .bam) with LoFreq..."
lofreq call-parallel --pp-threads 14 --no-default-filter --call-indels -C $RDP -f $Reference -q $Qual -Q $Qual -o $2$(basename $1 .bam)_LoFreq_RAW.vcf -l $Listregions $1 && \
bgzip -c $2$(basename $1 .bam)_LoFreq_RAW.vcf > $2$(basename $1 .bam)_LoFreq_inter.vcf.gz && bcftools index $2$(basename $1 .bam)_LoFreq_inter.vcf.gz && \
bcftools filter -O z -i "INFO/DP > ${RDP} & INFO/AF > ${VAF} & INFO/DP4[2]+INFO/DP4[3] > ${MRD}" $2$(basename $1 .bam)_LoFreq_inter.vcf.gz > $2$(basename $1 .bam)_LoFreq.vcf.gz && \
bcftools index $2$(basename $1 .bam)_LoFreq.vcf.gz && \
rm -rf $2$(basename $1 .bam)_LoFreq_inter.*
}

######
Sinvict() { #1: Input .bam file    #2: Output folder
echo "Analyzing $(basename $1 .bam) with bam-readcount && SiNVICT..."
mkdir $2/$(basename $1 .bam)_SiNVICT_inter_bamreadcount/
bam-readcount -l $Listregions -w 1 -b $Qual -f $Reference $1 > $2$(basename $1 .bam)_SiNVICT_inter_bamreadcount/output.bamreadcount && \

sinvict -m $RDP -t $2$(basename $1 .bam)_SiNVICT_inter_bamreadcount/ -o $2 && ./sinvict_to_vcf2.py $2/calls_level1.sinvict $2$(basename $1 .bam)_SiNVICT_inter.vcf $Reference && \
cat $2$(basename $1 .bam)_SiNVICT_inter.vcf | vcf-sort -c > $2$(basename $1 .bam)_SiNVICT_RAW.vcf && \
bgzip -c $2$(basename $1 .bam)_SiNVICT_RAW.vcf > $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz && \
bcftools index $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz && \
bcftools filter -O z -i "FORMAT/DP > ${RDP} & FORMAT/VAF > ${VAF} & FORMAT/MDP > ${MRD}" $2$(basename $1 .bam)_SiNVICT_inter_RAW.vcf.gz > $2$(basename $1 .bam)_SiNVICT.vcf.gz && \
bcftools index $2$(basename $1 .bam)_SiNVICT.vcf.gz && \
rm -rf $2*SiNVICT_inter* && rm -rf $2*alls_level*
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
      #samtools index $healthycontrol
      #Vardict $healthycontrol ./PoN/ & Mutect2 $healthycontrol ./PoN/ & Lofreq $healthycontrol ./PoN/ & Sinvict $healthycontrol ./PoN/ ####################
    done
    echo -e "\nMerging outputs for PoN..."
    #Merge per tool -> 4 seperate blacklists as .vcf file
    ls ./PoN/*Mutect2.vcf.gz | sed 's/^/ /' > ./PoN/Mutect2-normals.dat & ls ./PoN/*VarDict.vcf.gz | sed 's/^/ /' > ./PoN/VarDict-normals.dat & \
    ls ./PoN/*LoFreq.vcf.gz | sed 's/^/ /' > ./PoN/LoFreq-normals.dat & ls ./PoN/*SiNVICT.vcf.gz | sed 's/^/ /' > ./PoN/SiNVICT-normals.dat
    xargs -a ./PoN/Mutect2-normals.dat bcftools isec -o ./PoN/PoN_Mutect2.txt -O v -n +1 && xargs -a ./PoN/VarDict-normals.dat bcftools isec -o ./PoN/PoN_Vardict.txt -O v -n +1 && \
    xargs -a ./PoN/LoFreq-normals.dat bcftools isec -o ./PoN/PoN_Lofreq.txt -O v -n +1 && xargs -a ./PoN/SiNVICT-normals.dat bcftools isec -o ./PoN/PoN_Sinvict.txt -O v -n +1 && \
    ./isec_to_vcf.py ./PoN/PoN_Mutect2.txt ./PoN/PoN_Mutect2.vcf && ./isec_to_vcf.py ./PoN/PoN_Vardict.txt ./PoN/PoN_Vardict.vcf && \
    ./isec_to_vcf.py ./PoN/PoN_Lofreq.txt ./PoN/PoN_Lofreq.vcf && ./isec_to_vcf.py ./PoN/PoN_Sinvict.txt ./PoN/PoN_Sinvict.vcf
    #zipping and indexing Pon-per-tool blacklists for filtering later on;
    bgzip -c ./PoN/PoN_Mutect2.vcf > ./PoN/PoN_Mutect2.vcf.gz && bgzip -c ./PoN/PoN_Vardict.vcf > ./PoN/PoN_Vardict.vcf.gz && \
    bgzip -c ./PoN/PoN_Lofreq.vcf > ./PoN/PoN_Lofreq.vcf.gz && bgzip -c ./PoN/PoN_Sinvict.vcf > ./PoN/PoN_Sinvict.vcf.gz
    bcftools index ./PoN/PoN_Mutect2.vcf.gz && bcftools index ./PoN/PoN_Vardict.vcf.gz && bcftools index ./PoN/PoN_Lofreq.vcf.gz && bcftools index ./PoN/PoN_Sinvict.vcf.gz
    #We now have 4 PoNs, 1 per tool
    
    #if needed: Merging into 1 final blacklist occurs here...
    bcftools isec -o ./PoN/BLACKLIST_combined.txt -O v -n +1 ./PoN/PoN_Sinvict.vcf.gz ./PoN/PoN_Mutect2.vcf.gz ./PoN/PoN_Lofreq.vcf.gz ./PoN/PoN_Vardict.vcf.gz
    ./isec_to_vcf.py ./PoN/BLACKLIST_combined.txt ./PoN/BLACKLIST_combined.vcf && bgzip -c ./PoN/BLACKLIST_combined.vcf > ./PoN/BLACKLIST_combined.vcf.gz && bcftools index ./PoN/BLACKLIST_combined.vcf.gz 
      
      
  #One-File: a.k.a pre-made blacklist:      
  elif [[ -f $PoN ]]; then
    echo -e "PoN is given as pre-made blacklist: \n$PoN \n"
    cp $PoN ./PoN/BLACKLIST.txt
  
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
#samtools index $1
Vardict $1 "./output/$samplename/" & Mutect2 $1 "./output/$samplename/" & Lofreq $1 "./output/$samplename/" & Sinvict $1 "./output/$samplename/" ##################

#Merging, Annotating & Post-filtering:
if [ "$FilterMode" = "PerTool" ]; then
  ###  PerTool PoN Mode ###
  echo -e "\n\nPoN Filtering mode is set to 'PerTool'"
  #Filter first per tool on Pon,
  bcftools isec -p ./output/$samplename/Mutect2_PoN-filtered -O z -n~10 ./output/$samplename/*_Mutect2.vcf.gz ./PoN/PoN_Mutect2.vcf.gz
  bcftools isec -p ./output/$samplename/LoFreq_PoN-filtered -O z -n~10 ./output/$samplename/*_LoFreq.vcf.gz ./PoN/PoN_Lofreq.vcf.gz
  bcftools isec -p ./output/$samplename/SiNVICT_PoN-filtered -O z -n~10 ./output/$samplename/*_SiNVICT.vcf.gz ./PoN/PoN_Sinvict.vcf.gz
  bcftools isec -p ./output/$samplename/VarDict_PoN-filtered -O z -n~10 ./output/$samplename/*_VarDict.vcf.gz ./PoN/PoN_Vardict.vcf.gz
  #Merge pon-filtered-per-tool vcfs into one:
  bcftools isec -p ./output/$samplename/PerTool_PoN-filtered -O v -n +$Calls ./output/$samplename/SiNVICT_PoN-filtered/0000.vcf.gz ./output/$samplename/Mutect2_PoN-filtered/0000.vcf.gz \
  ./output/$samplename/LoFreq_PoN-filtered/0000.vcf.gz ./output/$samplename/VarDict_PoN-filtered/0000.vcf.gz
  #Add VAF/RD/MRD data to the merged list:
  echo "Merging sites.txt with vcfs to retreive data..."
  pythonscript ./merge_variants.py $samplename "PerTool"

elif [ "$FilterMode" = "Combined" ]; then
  ### Combined PoN mode: ### 
  echo -e "\n\nPoN Filtering mode is set to 'Combined'"
  #Merge outputs from all tools, zip and index:
  bcftools isec -p ./output/$samplename/Combined_PoN-filtered -O z -n +$Calls ./output/$samplename/*_SiNVICT.vcf.gz ./output/$samplename/*_Mutect2.vcf.gz ./output/$samplename/*_LoFreq.vcf.gz ./output/$samplename/*_VarDict.vcf.gz
  ./isec_to_vcf.py ./output/$samplename/Combined_PoN-filtered/sites.txt ./output/$samplename/Combined_PoN-filtered/sites.vcf
  bgzip -c ./output/$samplename/Combined_PoN-filtered/sites.vcf > ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz
  bcftools index ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz
  #Filter combined list with combined blacklist:
  bcftools isec -p ./output/$samplename/Combined_PoN-filtered/merged_PoN-filtered -O v -n~10 ./output/$samplename/Combined_PoN-filtered/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
  #Add VAF/RD/MRD data to the merged list:
  pythonscript ./merge_variants.py $samplename "Combined"

elif [ "$FilterMode" = "Both" ]; then
  ### Both PoN mode: ###
  echo -e "\n\nPoN Filtering mode is set to 'Both': PerTool and Combined lists are PoN filtered\nFiltering per tool first..."
  #Filter first per tool on Pon,
  bcftools isec -p ./output/$samplename/Mutect2_PoN-filtered -O z -n~10 ./output/$samplename/*_Mutect2.vcf.gz ./PoN/PoN_Mutect2.vcf.gz
  bcftools isec -p ./output/$samplename/LoFreq_PoN-filtered -O z -n~10 ./output/$samplename/*_LoFreq.vcf.gz ./PoN/PoN_Lofreq.vcf.gz
  bcftools isec -p ./output/$samplename/SiNVICT_PoN-filtered -O z -n~10 ./output/$samplename/*_SiNVICT.vcf.gz ./PoN/PoN_Sinvict.vcf.gz
  bcftools isec -p ./output/$samplename/VarDict_PoN-filtered -O z -n~10 ./output/$samplename/*_VarDict.vcf.gz ./PoN/PoN_Vardict.vcf.gz
  #Merge pon-filtered-per-tool vcfs into one, zip and index:
  echo -e "\nMerging PerTool filtered lists..."
  bcftools isec -p ./output/$samplename/PerTool_PoN-filtered -O z -n +$Calls ./output/$samplename/SiNVICT_PoN-filtered/0000.vcf.gz ./output/$samplename/Mutect2_PoN-filtered/0000.vcf.gz \
  ./output/$samplename/LoFreq_PoN-filtered/0000.vcf.gz ./output/$samplename/VarDict_PoN-filtered/0000.vcf.gz
  #Create vcf from merged PoN filtered sites.txt
  ./isec_to_vcf.py ./output/$samplename/PerTool_PoN-filtered/sites.txt ./output/$samplename/PerTool_PoN-filtered/sites.vcf
  bgzip -c ./output/$samplename/PerTool_PoN-filtered/sites.vcf > ./output/$samplename/PerTool_PoN-filtered/sites.vcf.gz
  bcftools index ./output/$samplename/PerTool_PoN-filtered/sites.vcf.gz
  #Filter merged list with combined blacklist:
  echo -e "\nFiltering merged list with combined PoN..."
  bcftools isec -p ./output/$samplename/PerTool_PoN-filtered/Both_PoN-filtered -O v -n~10 ./output/$samplename/PerTool_PoN-filtered/sites.vcf.gz ./PoN/BLACKLIST_combined.vcf.gz
  #Add VAF/RD/MRD data to the merged list:
  pythonscript ./merge_variants.py $samplename "Both"
  
else
  echo "$FilterMode is not a valid option for FilterMode... please provide either 'PerTool' or 'Combined'... "
  exit 1
fi

#MERGE SNVs on merged_vcfs.txt
echo "Removing duplicate SNVs that are part of MNVs and merging similar SNVs into MNVs..."
pythonscript ./SNV-MNV_handling.py $samplename $Reference


#Annotate merged_vcfs.txt
annotate ./output/$samplename/merged_vcfs_MNV.txt "merged_vcfs_annotated" excel
#Filtering annotated output
echo "Post filtering..." 
pythonscript ./post_filtering.py $samplename

#Annotating again if there are mutations left:
MUTS=$(awk '{ FS="\t" } { print $2 }' ./output/$samplename/filtering_info.txt | sed '10q;d')
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

