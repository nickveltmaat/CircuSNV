
"""
This module reads the merged_vcf.txt files from the tumor- and matched normal samples, compares VAF and exports tumor-only and higher-VAF-in-tumor variants.

Author: Nick Veltmaat
Date: 7-5-2025
"""

import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import argparse

parser = argparse.ArgumentParser(description="Filter merged_vcf.txt based on matched control")

parser.add_argument("tumor", type=str, help="Path to tumor merged_vcf.txt")
parser.add_argument("control", type=str, help="Path to matched normal merged_vcf.txt")
parser.add_argument("outfile", type=str, help="Output file (should be .txt")

args = parser.parse_args()

tumorpath = args.tumor
controlpath = args.control
outfile = args.outfile


#General Functions:

def read_sites(path_to_sites):
    '''
    This functions reads a sites.txt (output of bcftools isec) and returns a pandas dataframe
    '''
    df = pd.read_csv(path_to_sites, sep='\t', header=None)
    df[4] = df[4].astype(str)
    df[4] = df[4].apply(lambda x: x.zfill(4))
    df['loc'] = df[0].astype(str) + '_' + df[1].astype(str) + '_' + df[2] + '>' + df[3]
    
    df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS",'loc']
    df['VAF'] = df['TOOLS'].str.split('_AF:').str[1].str.split('_').str[0].astype(float)
    df['RDP'] = df['TOOLS'].str.split('_DP:').str[1].str.split('_').str[0].astype(int)
    df['MRD'] = df['TOOLS'].str.split('_MRD:').str[1].str.split('_').str[0].astype(int)
    return df


def chrom_sort_key(chrom):
    if pd.isna(chrom):
        return 25  # Missing chrom goes last
    if isinstance(chrom, (pd.Series, list)):
        chrom = chrom.iloc[0] if not chrom.empty else ''
    chrom = str(chrom).replace('chr', '')  # remove 'chr' if present
    if chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    else:
        try:
            return int(chrom)
        except ValueError:
            return 25  # unknown chromosomes go last



tumordf = read_sites(tumorpath)
controldf = read_sites(controlpath)

# Filter tumordf to include only rows where 'loc' is not in controldf
tumoronly = tumordf[~tumordf['loc'].isin(controldf['loc'])].copy()
tumoronly['somatic'] = True
# And vice versa
overlap = tumordf[tumordf['loc'].isin(controldf['loc'])].copy()
overlap = overlap.merge(controldf[['loc', 'VAF', 'RDP', 'MRD']], on='loc', suffixes=('_tumor', '_control'))

if overlap.empty:
  print("matched_filtering.py: No overlap between variants in tumor and control file, exiting...")
  sys.exit(0)




###############

# Function to compute and return all metrics
def somatic_metrics(row, alpha=0.05, min_diff=0.02, min_rel_ratio=3, max_germline_vaf=0.05):
    # Tumor and control counts
    tumor = [row['MRD_tumor'], row['RDP_tumor'] - row['MRD_tumor']]
    control = [row['MRD_control'], row['RDP_control'] - row['MRD_control']]
    
    # Fisher's exact test
    _, pvalue = fisher_exact([tumor, control], alternative='greater')
    
    # VAFs
    vaf_tumor = row['VAF_tumor']
    vaf_control = row['VAF_control']
    
    # Differences
    diff = vaf_tumor - vaf_control
    rel_ratio = vaf_tumor / vaf_control if vaf_control > 0 else float('inf')
    
    # 1) germline exclusion mask
    germline_ok = (vaf_control <= max_germline_vaf)
    # 2) somatic-enrichment mask (either abs OR rel)
    enriched    = (diff >= min_diff) or (rel_ratio >= min_rel_ratio)
    enriched    = (diff >= min_diff) or ((rel_ratio >= min_rel_ratio) and (diff >= 0.002))
    # 3) Fisher-significance mask
    sig         = (pvalue < alpha)
    
    somatic = bool(germline_ok and enriched and sig)
    
    return pd.Series({
        # leave these numeric
        'pvalue':    round(pvalue,    4),
        'diff':  round(diff,  4),
        'rel_ratio': round(rel_ratio, 2),
        # return true booleans here
        'germline_ok': germline_ok,
        'enriched':    enriched,
        'sig':         sig,
        'somatic':     somatic
    })


# Apply to dataframe
df_metrics = overlap.apply(somatic_metrics, axis=1)
# Combine with original df
overlap = pd.concat([overlap, df_metrics], axis=1)

###############


# Add missing columns to tumoronly, filling with NaN
missing_cols = [col for col in overlap.columns if col not in tumoronly.columns]
for col in missing_cols:
    tumoronly[col] = np.nan

tumoronly = tumoronly[overlap.columns]

df_final = pd.concat([overlap, tumoronly], ignore_index=True)


for c in ['RDP_tumor', 'MRD_tumor', 'RDP_control', 'MRD_control']:  # cast integer columns
    if c in df_final.columns:
        df_final[c] = df_final[c].astype('Int64')

for c in ['germline_ok', 'enriched', 'sig', 'somatic']:  # cast boolean columns
    if c in df_final.columns:
        df_final[c] = df_final[c].astype('boolean')


#Re-sort and save the final df
if not df_final.empty and 'CHROM' in df_final.columns:
    df_final['chrom_sort'] = df_final['CHROM'].apply(chrom_sort_key)
    df_final = df_final.sort_values(by=['chrom_sort', 'POS']).drop(columns=['chrom_sort']).reset_index(drop=True)
else:
    df_final = df_final.reset_index(drop=True)  # Just reset index if nothing to sort

df_final.to_csv(outfile.replace(".txt", "_matched_fltr.txt"), sep='\t', index=False)

# Filter for somatic mutations only
if 'somatic' in df_final.columns:
    df_somatic_only = df_final[df_final['somatic'] == True]
    df_somatic_only = df_somatic_only[['CHROM', 'POS', 'REF', 'ALT', 'TOOLS']]
else:
    df_somatic_only = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'TOOLS'])

df_somatic_only.to_csv(outfile, sep='\t', index=False, header=False)

