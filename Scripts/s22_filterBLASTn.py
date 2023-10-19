#!/usr/bin/env python3

"""Filter the BLASTn output files by size


Usage: python3 s22_filterBLASTn.py <BLASTn file> <output file>

Args:
	BLASTn file: BLASTn output file with hits/matches
  output file: file where the filtered outputs should be storred in
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 3:
  print(__doc__)
  exit()

input_file = sys.argv[1]
output_path = sys.argv[2]

output_file = open(output_path, "w")

for line in open(input_file, "r"):
    line = line.strip()
    length = line.split("\t")[3]
    if int(length) >= 500:
      print(line, file = output_file)

output_file.close()