"""
This module reads a 'bcftools isec' output file of merged variants (sites.txt) and adds the data from the raw VCF files to it (RD, VAF & MRD)

Author: Nick Veltmaat
Date: 14-2-2023
"""

import pandas as pd
import io
import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description="Merge VCF variant data with bcftools isec output.")

parser.add_argument("samplename", type=str, help="Name of the sample.")
parser.add_argument("filtermode", type=str, help="Filtering mode: 'phased', 'phased_pon', etc.")
parser.add_argument("--infile", type=str, default=None, help="Input file (optional, depending on mode).")
parser.add_argument("--outfile", type=str, default=None, help="Output file (optional, depending on mode).")
parser.add_argument("outdir", type=str, help="Output directory.")

args = parser.parse_args()

samplename = args.samplename
filtermode = args.filtermode
infile = args.infile
outfile = args.outfile
outdir = args.outdir
if not outdir.endswith('/'):
    outdir += '/'
    
#General Functions:
#####
def read_vcf(path):
    '''
    This functions reads a .vcf file and returns a pandas dataframe
    '''
    with gzip.open(path, 'rt') as f:
      lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)), 
                       dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 
                              'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str}, 
                       sep='\t').rename(columns={'#CHROM': 'CHROM'})

def read_sites(path_to_sites):
    '''
    This functions reads a sites.txt (output of bcftools isec) and returns a pandas dataframe
    '''
    df = pd.read_csv(path_to_sites, sep='\t', header=None)
    df[4] = df[4].astype(str)
    df[4] = df[4].apply(lambda x: x.zfill(4))
    df['loc'] = df[0].astype(str) + '_' + df[1].astype(str) + '_' + df[2] + '>' + df[3]
    df.columns = ["CHROM", "POS", 'REF', 'ALT', "TOOLS",'loc']
    return df

#Generating RD, VAF & MDP columns from .vcf
def vcf_info_vardict(df):
    '''
    This functions loads a Vardict .vcf df and extracs 'RD', 'VAF' & 'MRD' data into new columns
    '''
    df['MRD'] = df[df.columns[9]].str.split(':').str[2]
    df['AF'] = df[df.columns[9]].str.split(':').str[4]
    df['RD'] = df[df.columns[9]].str.split(':').str[1]
    return df

def vcf_info_lofreq(df):
    '''
    This functions loads a Lofreq .vcf df and extracs 'RD', 'VAF' & 'MRD' data into new columns
    '''
    df['AF'] = df['INFO'].str.split('AF=').str[1].str.split(';SB').str[0].astype(float)
    df['RD'] = df['INFO'].str.split('DP=').str[1].str.split(';AF').str[0].astype(int)
    df['MRD'] = df['INFO'].str.split('DP4=').str[1].str.split(';').str[0].str.split(',').str[2:].apply(lambda x: sum([int(i) for i in x])).astype(int)
    return df

def vcf_info_mutect(df):
    '''
    This functions loads a Mutect2 .vcf df and extracs 'RD', 'VAF' & 'MRD' data into new columns
    '''
    df['AF'] = df[df.columns[9]].str.split(':').str[2].astype(float)
    df['RD'] = df[df.columns[9]].str.split(':').str[3].astype(int)
    df['MRD'] = df[df.columns[9]].str.split(':').str[1].str.split(',').str[1]
    return df

def vcf_info_sinvict(df):
    '''
    This functions loads a Sinvict .vcf df and extracs 'RD', 'VAF' & 'MRD' data into new columns
    '''
    #df['MRD'] = df['calls_level1.sinvict'].str.split(':',4).str[2].astype(int)
    df['MRD'] = df['calls_level1.sinvict'].str.split(':').str[2].astype(int)
    df['RD'] = df['calls_level1.sinvict'].str.split(':').str[0].astype(int)
    df['AF'] = df['calls_level1.sinvict'].str.split(':').str[1].astype(float)
    return df
########


# Reading vcf files:
if filtermode == 'no_pon':
    sites = read_sites(outdir + samplename + '/isec/sites.txt') #Reads PerTool-PoN-filtered sites
    dfsv = read_vcf(outdir + samplename + '/isec/0000.vcf.gz') #reads .vcf files of merged vcfs to retreive RD, VAF & MRD info
    dfm2 = read_vcf(outdir + samplename + '/isec/0001.vcf.gz')
    dflf = read_vcf(outdir + samplename + '/isec/0002.vcf.gz')
    dfvd = read_vcf(outdir + samplename + '/isec/0003.vcf.gz')
elif filtermode == 'pon':
    sites = read_sites(outdir + samplename + '/isec/PoN-filtered/sites.txt') #Reads Combined-PoN-filtered sites
    sites2 = read_sites(outdir + samplename + '/isec/sites.txt') #Reads unfiltered sites, to see which tools called the mutation. 
    sites['TOOLS'] = pd.merge(sites, sites2, on='loc')['TOOLS_y']  #replaces 'TOOLS column'
    dfsv = read_vcf(outdir + samplename + '/isec/0000.vcf.gz') #reads .vcf files of merged vcfs to retreive RD, VAF & MRD info
    dfm2 = read_vcf(outdir + samplename + '/isec/0001.vcf.gz')
    dflf = read_vcf(outdir + samplename + '/isec/0002.vcf.gz')
    dfvd = read_vcf(outdir + samplename + '/isec/0003.vcf.gz')
elif filtermode == 'phased_pon':
    sites = read_sites(outdir + samplename + '/mnvcalling/ponfiltered/sites.txt') #Reads unfiltered sites.
    sites2 = read_sites(outdir + samplename + '/mnvcalling/ponfiltered2/sites.txt') #Reads unfiltered sites. 
    dfsv = read_vcf(outdir + samplename + '/mnvcalling/SiNVICT_MNV.vcf.gz') #reads .vcf files of merged vcfs to retreive RD, VAF & MRD info
    dfm2 = read_vcf(outdir + samplename + '/mnvcalling/Mutect2_MNV.vcf.gz')
    dflf = read_vcf(outdir + samplename + '/mnvcalling/LoFreq_MNV.vcf.gz')
    dfvd = read_vcf(outdir + samplename + '/mnvcalling/VarDict_MNV.vcf.gz')
elif filtermode == 'phased':
    sites = read_sites(outdir + samplename + '/mnvcalling/' + infile) #Reads unfiltered sites. 
    dfsv = read_vcf(outdir + samplename + '/mnvcalling/SiNVICT_MNV.vcf.gz') #reads .vcf files of merged vcfs to retreive RD, VAF & MRD info
    dfm2 = read_vcf(outdir + samplename + '/mnvcalling/Mutect2_MNV.vcf.gz')
    dflf = read_vcf(outdir + samplename + '/mnvcalling/LoFreq_MNV.vcf.gz')
    dfvd = read_vcf(outdir + samplename + '/mnvcalling/VarDict_MNV.vcf.gz')
    
#Elif filtermode = non existant...
    
#Extracting RD, VAF & MRD data per tool.vcf
dfsv = vcf_info_sinvict(dfsv).rename(columns = {'AF':'AFsv', 'RD':'RDsv', 'MRD':'MRDsv'})
dfsv['loc'] = dfsv['CHROM'] + '_' + dfsv['POS'].astype(str) + '_' + dfsv['REF'] + '>' + dfsv['ALT']
dfsv = dfsv.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO', 'FORMAT', 'calls_level1.sinvict'])

dfm2 = vcf_info_mutect(dfm2).rename(columns = {'AF':'AFm2', 'RD':'RDm2', 'MRD':'MRDm2'}).drop(dfm2.columns[9], axis=1)
dfm2['loc'] = dfm2['CHROM'] + '_' + dfm2['POS'].astype(str) + '_' + dfm2['REF'] + '>' + dfm2['ALT']
dfm2 = dfm2.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL',  'INFO', 'FORMAT', 'FILTER'])

dflf = vcf_info_lofreq(dflf).rename(columns = {'AF':'AFlf', 'RD':'RDlf', 'MRD':'MRDlf'})
dflf['loc'] = dflf['CHROM'] + '_' + dflf['POS'].astype(str) + '_' + dflf['REF'] + '>' + dflf['ALT']
dflf = dflf.drop(columns =['CHROM', 'POS','REF', 'ALT', 'ID', 'FILTER', 'INFO', 'QUAL'])

dfvd = vcf_info_vardict(dfvd).rename(columns = {'AF':'AFvd', 'RD':'RDvd', 'MRD':'MRDvd'}).drop(dfvd.columns[9], axis=1)
dfvd['loc'] = dfvd['CHROM'] + '_' + dfvd['POS'].astype(str) + '_' + dfvd['REF'] + '>' + dfvd['ALT']
dfvd = dfvd.drop(columns =['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL', 'INFO', 'FORMAT'])


#Merging dfs per tool with sites.txt
df_merged = pd.merge(pd.merge(pd.merge(pd.merge(dfm2, dflf, on='loc', how='outer'), 
                                       dfsv, on='loc', how='outer'), 
                              dfvd, on='loc', how='outer'), 
                     sites, on='loc', how='inner')







def process_dataframe(df_merged):
  df_merged["AF_mean"] = df_merged[["AFm2", "AFlf", "AFsv", "AFvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(4)
  df_merged["DP_mean"] = df_merged[["RDm2", "RDlf", "RDsv", "RDvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(0).astype(int)
  df_merged["MRD_mean"] = df_merged[["MRDm2", "MRDlf", "MRDsv", "MRDvd"]].apply(pd.to_numeric, args=['coerce']).mean(axis=1, skipna = True).round(0).astype(int)
  df_merged['info'] = 'SV-M2-LF-VD:' + df_merged['TOOLS'] + '_AF:' + df_merged['AF_mean'].astype(str) + '_DP:' + df_merged['DP_mean'].astype(str)+ '_MRD:' + df_merged['MRD_mean'].astype(str)
  df_merged = df_merged[['CHROM', 'POS','REF','ALT','info']].copy()
  
  #Sort on position (chrom, pos)
  df_merged['chrnum'] = df_merged['CHROM'].str.extract(r'(\d+|[XY])').replace({'X': 23, 'Y': 24}).astype(int)
  df_merged = df_merged.sort_values(by=['chrnum', 'POS']).drop(columns =['chrnum'])
    
  return df_merged

    
    
# Example usage:
df_merged = process_dataframe(df_merged).copy()




if filtermode == 'phased_pon':
  df_merged2 = pd.merge(pd.merge(pd.merge(pd.merge(dfm2, dflf, on='loc', how='outer'), 
                                         dfsv, on='loc', how='outer'), 
                                dfvd, on='loc', how='outer'), 
                       sites2, on='loc', how='inner')
  df_merged2 = process_dataframe(df_merged2).copy()
else:
  pass


if filtermode == 'no_pon':
  df_merged.to_csv(outdir + samplename + '/merged_vcfs.txt', sep='\t', index = False, header= False)
elif filtermode == 'pon':
  df_merged.to_csv(outdir + samplename + '/merged_vcfs.txt', sep='\t', index = False, header= False)
elif filtermode == 'phased':
  df_merged.to_csv(outdir + samplename + '/mnvcalling/' + outfile, sep='\t', index = False, header= False)
elif filtermode == 'phased_pon':
  df_merged.to_csv(outdir + samplename + '/mnvcalling/merged_vcfs_phased_mnv.txt', sep='\t', index = False, header= False) #SNP filtered --> pon filtered
  df_merged2.to_csv(outdir + samplename + '/mnvcalling/sites_pon_filtered.txt', sep='\t', index = False, header= False) #Is extra: Only PoN filtered

 
