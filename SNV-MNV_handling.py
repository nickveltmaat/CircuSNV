#! /usr/bin/python

"""
This module reads a merged bcftools isec output from 'merge_variants.py' and does two things to reduce the amount of duplicate calls for the same mutational event
  - Remove duplicate SNVs that are part of MNV (most often called by another tool)
  - Combine close-proximity SNVs into MNVs (Must have similar Mutant Read Depth, as a control that these variants are actually the same mutational event)

Author: Nick Veltmaat
Date: 17-2-2023
"""

import pandas as pd
import sys
import subprocess

if len(sys.argv) != 3:
  print("Usage:\t" + sys.argv[0] + "\t<samplename>\t<reference_fasta>")
  exit(0)
  
samplename = sys.argv[1]
ref_genome_file = sys.argv[2]

path = './output/'

mergedvcf = pd.read_csv(path + samplename + '/merged_vcfs.txt', sep='\t', header=None)

mergedvcf['tools'] = mergedvcf[4].str.split('VD:', 2).str[1].str.split('_', 2).str[0]
mergedvcf['AF'] = mergedvcf[4].str.split('AF:', 2).str[1].str.split('_', 2).str[0].astype(float)
mergedvcf['RD'] = mergedvcf[4].str.split('DP:', 2).str[1].str.split('_', 2).str[0].astype(int)
mergedvcf['MRD'] = mergedvcf[4].str.split('MRD:', 2).str[1].astype(int)
mergedvcf.columns = ['CHR', 'POS', 'REF', 'ALT', 'remove','tools', 'AF', 'RD', 'MRD']
mergedvcf = mergedvcf.drop('remove', axis=1)
mergedvcf['loc'] = mergedvcf['CHR'] + '_' + mergedvcf['POS'].astype(str) + '_' + mergedvcf['REF'] + '>' + mergedvcf['ALT']

snvdf = mergedvcf[(mergedvcf['REF'].apply(len) == 1) & (mergedvcf['ALT'].apply(len) == 1)].copy()

mnvdf = mergedvcf[(mergedvcf['REF'].apply(len) > 1) & (mergedvcf['ALT'].apply(len) > 1)].copy()
mnvdf['MNV_length'] = mnvdf['REF'].apply(len)
mnvdf['MNV_altlength'] = mnvdf['ALT'].apply(len)
mnvdf = mnvdf[mnvdf['MNV_altlength'] == mnvdf['MNV_length']].drop(['MNV_altlength', 'MNV_length'], axis=1)


#Generate dataframe of SNVs that are part of MNVs (mnvdf -> snvs that can be removed if present in original list)
if len(mnvdf) > 0:
    #Generate dataframe of SNVs that are part of MNVs (mnvdf -> snvs that can be removed if present in original list)
    df = mnvdf.copy()

    refs = df['REF'].apply(list).apply(pd.Series, 1).stack()
    alts = df['ALT'].apply(list).apply(pd.Series, 1).stack()

    refs.index = refs.index.droplevel(-1) # to line up with df's index
    alts.index = alts.index.droplevel(-1)
    refs.name = 'REF' # needs a name to join
    alts.name = 'ALT'
    del df['REF']
    del df['ALT']
    df2 = df.join(refs)
    df2['ALT'] = alts

    df2 = df2[df2['REF'].astype(bool)]    

    POSs = set(df2['POS']) #Create set of positions

    dfs = [] #Add each variant with same position to its own df and store in list
    for i in POSs:
        dfs.append(df2[df2['POS'] == i]) 

    dfs2 = [] #Increase the number of POS for every generated SNV
    for i in dfs:
        df = i.copy()
        ln = len(df)
        df['ln'] = list(range(ln))
        df['POS'] = df['POS'] + df['ln']
        df = df.drop(columns='ln')
        dfs2.append(df)

    mnvstosnv = pd.concat(dfs2) #Combine into one df
    mnvstosnv = mnvstosnv[mnvstosnv['REF'] != mnvstosnv['ALT']] #Remove non-variants (A>A)
    mnvstosnv['loc'] = mnvstosnv['CHR'] + '_' + mnvstosnv['POS'].astype(str) + '_' + mnvstosnv['REF'] + '>' + mnvstosnv['ALT']

    #Compare loc with SNVs from merged_vcfs and remove SNVs from merged_vcfs which are duplicate
    tokeep = mergedvcf[~mergedvcf['loc'].isin(mnvstosnv['loc'])]

elif len(mnvdf) == 0:
    tokeep = mergedvcf
else:
  pass

### PART 2: Combine SNVs into MNVs ###
# Read the input dataframe from a file
merged_df = tokeep.copy()

#Sort on position (chrom, pos)
merged_df['chrnum'] = merged_df['CHR'].str.extract(r'(\d+|[XY])').replace({'X': 23, 'Y': 24}).astype(int)
merged_df = merged_df.sort_values(by=['chrnum', 'POS']).drop(columns =['chrnum'])

#Reshuffle columns:
merged_df['info'] = 'SV-M2-LF-VD:' + merged_df['tools'] + '_AF:' + merged_df['AF'].astype(str) + '_DP:' + merged_df['RD'].astype(str)+ '_MRD:' + merged_df['MRD'].astype(str)
merged_df = merged_df[['CHR', 'POS','REF','ALT','info']] 

# Save the merged variants to a new file
merged_df.to_csv(path + samplename + '/merged_vcfs_MNV.txt', sep='\t', index = False, header= False)
