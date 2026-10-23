#!/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/python3

"""
This module loads the generated PoN and adds mutations that are solely found with one tool, but in at least x % of the normal samples. 
Returns a sorted PoN. 

Author: Nick Veltmaat
Date: 26-7-2024
"""

import sys
import pandas as pd
import numpy as np


if len(sys.argv) != 3: #Changed
  print("Usage:\t" + sys.argv[0] + "\t minimal percentage of normal samples that should contain a mutation to be added to PoN, that were only found using 1 tool" + "\t<path/to/PoN/directory>") #Changed
  exit(0)
  
percentage = int(sys.argv[1])
path = str(sys.argv[2])


pon = pd.read_csv(path + '/BLACKLIST_combined.txt', sep='\t', header=None)
pon = pd.read_csv(path + '/BLACKLIST_combined.txt', sep='\t', header=None, dtype={len(pon.columns) - 1: str})
pon_original = pon.copy()
pon['mut'] = pon[0] +"_" + pon[1].astype(str) +"_" + pon[2] +">" + pon[3]


def g(toolpon):
    df = pd.read_csv(toolpon, sep='\t', header=None)
    df = pd.read_csv(toolpon, sep='\t', header=None, dtype={len(df.columns)-1: str})
    df['mut'] = df[0] +"_" + df[1].astype(str) +"_" + df[2] +">" + df[3]
    
    split_df = df[4].apply(lambda x: pd.Series(list(x)))
    split_df.columns = [f'{6+i}' for i in range(len(split_df.columns))]
    
    # Concatenate the DataFrame
    newdf = pd.concat([df, split_df], axis=1)
    newdf.iloc[:, 6:] = newdf.iloc[:, 6:].astype(int)
    Nnormals = len(newdf.iloc[:, 6:].columns)
    newdf['sum'] = newdf.iloc[:, 6:].sum(axis=1)
    newdf['N'] = Nnormals
    newdf['pct'] = round(newdf['sum'] / newdf['N'] * 100, 10)
    
    #Keep only mutations present in X % or more normal samples
    newdf = newdf[newdf['pct'] >= percentage]
    
    #Keep only mutations that are not already in the original pon
    newdf = newdf[~newdf['mut'].isin(set(pon['mut']))]
    
    return newdf


sinvictpon = g(path + '/PoN_Sinvict.txt')
sinvictpon = sinvictpon[[0,1,2,3]]
sinvictpon[4] = '1000'

mutectpon = g(path + '/PoN_Mutect2.txt')
mutectpon = mutectpon[[0,1,2,3]]
mutectpon[4] = '0100'

lofreqpon = g(path + '/PoN_Lofreq.txt')
lofreqpon = lofreqpon[[0,1,2,3]]
lofreqpon[4] = '0010'

vardictpon = g(path + '/PoN_Vardict.txt')
vardictpon = vardictpon[[0,1,2,3]]
vardictpon[4] = '0001'

complete_pon = pd.concat([pon_original, lofreqpon, sinvictpon, mutectpon, vardictpon], ignore_index=True)

### Sorting final list
def chromosome_sort_key(chr_str):
    if chr_str.startswith('chr'):
        try:
            return int(chr_str[3:])
        except ValueError:
          # Handle special cases like 'chrX', 'chrY'
            if chr_str == 'chrX':
                return 23
            elif chr_str == 'chrY':
                return 24
            else:
                return 99  #for unexpected cases

    return 99  # For non-chr prefixes

def sort_df(df):

  df['chr_sort_key'] = df[0].apply(chromosome_sort_key)
  df = df.sort_values(['chr_sort_key', 1])
  df = df.drop('chr_sort_key', axis=1)
  return(df)

complete_pon_sorted = sort_df(complete_pon)

complete_pon_sorted.to_csv(path + '/BLACKLIST_combined_complete.txt', sep='\t', index = False, header= False) 
