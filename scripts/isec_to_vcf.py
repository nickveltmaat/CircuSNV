#! /usr/bin/python


"""
This module reads a 'bcftools isec' output file of merged variants, and then converts to standard .vcf format.

Author: Nick Veltmaat
Date: 10-2-2023
"""

import sys
import re

if len(sys.argv) != 4:
  print("Usage:\t" + sys.argv[0] + "\t<path/to/input_file>\t<path/to/output_filename>" + "\t<path/to/input_file>\t<path/to/output_filename>" + 'Parse Info (TRUE/FALSE)') #Changed
  sys.exit(1)

outfile = open(sys.argv[2], "w")
outfile.write("##fileformat=VCFv4.3\n")

if sys.argv[3].lower() == 'true':
  outfile.write('''##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Read Depth">\n##INFO=<ID=MRD,Number=1,Type=Integer,Description="Mutation supporting Read Depth">\n''')
else:
  pass

#Printing chr contigs in order
def custom_sort_key(chrom):
    match = re.match(r'chr(\d+|X|Y)', chrom)
    if match:
        value = match.group(1)
        if value.isdigit():
            return int(value)
        elif value == 'X':
            return 23
        elif value == 'Y':
            return 24
    return chrom

chroms = []
with open(sys.argv[1]) as infile:
    for line in infile:
        line = line.rstrip()
        tokens = line.split()
        chromosome = tokens[0]
        chroms.append(chromosome)
    
chroms2 = sorted(set(chroms), key=custom_sort_key)

for i in chroms2:
  outfile.write("##contig=<ID=" + str(i) + ">\n")

#Printing vcf lines
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    position = tokens[1]
    ref = tokens[2]
    alt = tokens[3]
    if sys.argv[3].lower() == 'true':
      info = tokens[4]
      splitinfo = info.split('_', 4)
      vaf = splitinfo[1].split(':')[1]
      dp = splitinfo[2].split(':')[1]
      mdp = splitinfo[3].split(':')[1]
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + "AF="+vaf+";DP="+dp+";MRD="+mdp + "\n")
    elif sys.argv[3].lower() == 'false':
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + ".\n")
    else:
      print('Parse TRUE or FALSE to determine whether you want the INFO column kept or not')
