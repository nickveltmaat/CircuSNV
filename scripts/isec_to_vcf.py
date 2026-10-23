#! /usr/bin/python

#!/groups/umcg-pmb/tmp01/projects/hematopathology/Lymphoma/Nick/SNVcallingPipeline/env/bin/python3

"""
This module reads a 'bcftools isec' output file of merged variants, and then converts to standard .vcf format.

Author: Nick Veltmaat
Date: 10-2-2023
"""

import sys

if len(sys.argv) != 4: #Changed
  print("Usage:\t" + sys.argv[0] + "\t<path/to/input_file>\t<path/to/output_filename>" + "\t<path/to/input_file>\t<path/to/output_filename>" + 'Parse Info (TRUE/FALSE)') #Changed
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
    
chroms2 = set(chroms)
for i in chroms2:
    outfile.write("##contig=<ID="+str(i)+">\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

with open(sys.argv[1]) as infile:
  for line in infile:
    line = line.rstrip()
    tokens = line.split()
    chromosome = tokens[0]
    position = tokens[1]
    ref = tokens[2]
    alt = tokens[3]
    #outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS"  + "\t" + ".\n")
    if sys.argv[3].lower() == 'true':
      info = tokens[4] #changed
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + info + "\n")
    elif sys.argv[3].lower() == 'false':
      outfile.write(chromosome + "\t" + str(position) + "\t" + "." + "\t" + ref.upper() + "\t" + alt.upper() + "\t" + "." + "\t" + "PASS" + "\t" + ".\n")
    else:
      print('Parse TRUE or FALSE to determine whether you want the INFO column kept or not')
