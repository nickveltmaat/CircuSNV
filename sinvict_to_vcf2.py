#! /usr/bin/python

"""
This module reads a raw SiNVICT output file and converts to standard .vcf format.

Author: Nick Veltmaat
Date: 6-2-2023
"""

import os
import sys
import subprocess

if len(sys.argv) != 4:
  print("Usage:\t" + sys.argv[0] + "\t<input_file>\t<output_filename>\t<reference_fasta>")
  exit(0)
  
outfile = open(sys.argv[2], "w")
outfile.write("##fileformat=VCFv4.3\n")

chroms = []
with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    chroms.append(chromosome)

outfile.write('''##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n''')
outfile.write('''##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">\n''')
outfile.write('''##FORMAT=<ID=MDP,Number=1,Type=Integer,Description="Number of reads supporting the mutation">\n''')

chroms2 = set(chroms)
    
for i in chroms2:
    outfile.write("##contig=<ID="+str(i)+">\n")

outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+os.path.basename(sys.argv[1])+"\n")

with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    position = tokens[1]
    ref = tokens[3]
    alt = tokens[5]
    vaf = round(float(tokens[7])/100, 3)
    depth = int(tokens[4])
    mutdp = int(tokens[6])
    if alt[0] == "+" :
      #insertion
      alt = ref + alt[1:]
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp)  + "\n")
    elif alt[0] == "-" :
      #deletion
      position = int(position)
      position -= 1
      fa = subprocess.check_output(["samtools", "faidx", sys.argv[3], chromosome+":"+str(position)+"-"+str(position)])
      base = fa.split("\n".encode())[1]
      base.rstrip()
      ref = base + str(alt[1:]).encode()
      alt = base
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.decode("utf-8").upper() + "\t" + alt.decode("utf-8").upper() + "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp) + "\n")
    else: 
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS"  + "\t" + "." + "\t" + "DP:VAF:MDP" + "\t" + str(depth)+':'+str(vaf)+':'+str(mutdp) + "\n")

# original copied and modified at 6 feb 2023 from: https://github.com/sfu-compbio/sinvict/blob/master/sinvict_to_vcf.py
