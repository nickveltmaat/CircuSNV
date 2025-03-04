#! /usr/bin/python

"""
This module reads a annotated (using OpenCravat) merged_variants sites.xlsx.
 
Variants are filtered on functional effect (non-synonymous, UTR variants, intronic variants etc...)
Output is a filtered sites.txt that needs to be annotated again, as well as filtering info (text file showing remaining variants after each filtering step) and an excel file containing filtered-out variants

Author: Nick Veltmaat
Date: 15-2-2023
"""

import sys
import pandas as pd
import os

if len(sys.argv) != 2:
  print("Usage:\t" + sys.argv[0] + "\t<samplename>")
  exit(0)

samplename = str(sys.argv[1])
path = './output/'



#Load annotated variants
annotated = pd.read_excel(path + samplename + '/merged_vcfs_annotated.xlsx', sheet_name='Variant', skiprows=1, engine='openpyxl')
annotated['VAF'] = annotated['Tags'].str.split('AF:', 1).str[1].str.split('_').str[0].astype(float)
annotated['DP'] = annotated['Tags'].str.split('DP:', 1).str[1].str.split('_').str[0].astype(int)
annotated['MRD'] = annotated['Tags'].str.split('MRD:', 1).str[1].str.split(';').str[0].astype(int)

#Initiate filter counts df
filtereddf = pd.DataFrame(columns=['filter', 'variants_remaining'])
filtereddf.loc[len(filtereddf)] = ['start', len(pd.read_csv(path + samplename + '/merged_vcfs.txt', sep='\t', header=None))]
if os.path.exists(path + samplename + '/merged_vcfs_MNV.txt'): #check if MNV handling was performed
    filtereddf.loc[len(filtereddf)] = ['SNV/MNV', len(pd.read_csv(path + samplename + '/merged_vcfs_MNV.txt', sep='\t', header=None))]
else:
    pass

###START FILTERING & recording filtered variants
#Synonymous
filteredvariants1 = annotated[annotated['Sequence Ontology'] == 'synonymous_variant'].copy() #Save synonymous variants in other df
filteredvariants1['FilterReason'] = 'Synonymous'
annotated = annotated[annotated['Sequence Ontology'] != 'synonymous_variant'] #Filter synonymous
filtereddf.loc[len(filtereddf)] = ['Synonymous', len(annotated)]

#3' variants
filteredvariants2 = annotated[annotated['Sequence Ontology'] == '3_prime_UTR_variant'].copy()
filteredvariants2['FilterReason'] = "3' UTR"
annotated = annotated[annotated['Sequence Ontology'] != '3_prime_UTR_variant'] #Filter 3' UTR variants
#5' variants
filteredvariants3 = annotated[annotated['Sequence Ontology'] == '5_prime_UTR_variant'].copy()
filteredvariants3['FilterReason'] = "5' UTR"
annotated = annotated[annotated['Sequence Ontology'] != '5_prime_UTR_variant'] #Filter 5' UTR variants
filtereddf.loc[len(filtereddf)] = ["5'/3' UTR", len(annotated)]

#Intronic
filteredvariants4 = annotated[annotated['Sequence Ontology'] == 'intron_variant'].copy()
filteredvariants4['FilterReason'] = "Intronic"
annotated = annotated[annotated['Sequence Ontology'] != 'intron_variant'] #Filter Intronic variants
filtereddf.loc[len(filtereddf)] = ['Intronic', len(annotated)]

#2kb upstream
filteredvariants5 = annotated[annotated['Sequence Ontology'] == '2kb_upstream_variant'].copy()
filteredvariants5['FilterReason'] = "2kb upstream"
annotated = annotated[annotated['Sequence Ontology'] != '2kb_upstream_variant'] #Filter 2kb Upstream variants
#2kb Downstream
filteredvariants6 = annotated[annotated['Sequence Ontology'] == '2kb_downstream_variant'].copy()
filteredvariants6['FilterReason'] = "2kb downstream"
annotated = annotated[annotated['Sequence Ontology'] != '2kb_downstream_variant'] #Filter 2kb Downstream variants
filtereddf.loc[len(filtereddf)] = ['2kb up-/downstream', len(annotated)]

#in SNP db > 0.1% AF
filteredvariants7 = annotated[((annotated['AF'] > 0.0001) | #1000G global
                (annotated['EUR AF'] > 0.0001) | #1000G European
                (annotated['Non-Fin Eur AF'] > 0.0001) | #GnomAD European
                (annotated['Global AF'] > 0.0001))].copy() #GnomAD global
                
filteredvariants7['FilterReason'] = "> 0.01% in GnoMAD3 / 1000G"
#annotated = annotated.loc[(annotated['Global AF'].fillna(0) <= 0.0001) | (annotated['Non-Fin Eur AF'].fillna(0) <= 0.0001) | (annotated['AF'] <= 0.0001) | (annotated['EUR AF'] <= 0.0001)] 
annotated = annotated[~((annotated['AF'] > 0.0001) | #1000G global
                (annotated['EUR AF'] > 0.0001) | #1000G European
                (annotated['Non-Fin Eur AF'] > 0.0001) | #GnomAD European
                (annotated['Global AF'] > 0.0001))] #GnomAD global
#Filter out SNPS with db VAF of > 1%


filtereddf.loc[len(filtereddf)] = ['> 0.01% in GnoMAD3 / 1000G', len(annotated)]  

#VAF < 0.5%
#filteredvariants8 = annotated[annotated['VAF'] < 0.005].copy()
#filteredvariants8['FilterReason'] = "VAF < 0.5%"
#annotated = annotated[annotated['VAF'] >= 0.005]
#filtereddf.loc[len(filtereddf)] = ['VAF < 0.5%', len(annotated)]    
    
#VAF < 1.0%
#filteredvariants9 = annotated[annotated['VAF'] < 0.01].copy()
#filteredvariants9['FilterReason'] = "VAF < 1%"
#annotated = annotated[annotated['VAF'] >= 0.01]
#filtereddf.loc[len(filtereddf)] = ['VAF < 1%', len(annotated)] 



#Merging all filtered-out variants in one dataframe
#filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, filteredvariants4, 
#                        filteredvariants5, filteredvariants6, filteredvariants7, filteredvariants8, filteredvariants9]    
filteredvariantslist = [filteredvariants1, filteredvariants2, filteredvariants3, filteredvariants4, 
                        filteredvariants5, filteredvariants6, filteredvariants7]                            
                        
filteredvariants = pd.concat(filteredvariantslist)

#Saving files: 
filtereddf.to_csv(path + samplename + '/filtering_info.txt', sep='\t', index = False, header= True) #Filtering info
annotated['sample'] = samplename
annotated[['Chrom.1', 'Pos', 'Reference allele', 'Alternate allele', 'sample','Tags']].to_csv(path + samplename + '/merged_vcfs_annotated-filtered.txt', sep= '\t', index = False, header= False) #Filtered merged variants
filteredvariants.to_excel(path + samplename + '/removed_variants.xlsx', index=False) #Filterd-out variants
