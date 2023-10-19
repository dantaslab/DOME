#!/usr/bin/env python3

"""Filter out the duplicates (A B vs B A) from the joint BLASTn output file


Usage: python3 s22_nodupBLASTn.py <BLASTn file> <filtered output>

Args:
	BLASTn file: joint BLASTn output file; .txt file
  filtered output: outputs with duplicates filtered out; .txt file
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 3:
  print(__doc__)
  exit()

input_file = sys.argv[1]
output_path = sys.argv[2]

pairs = []

BLASTn_nodup_file = open(output_path, "w")

for line in open(input_file, "r"):
  line = line.strip()
  contig1 = line.split("\t")[0]
  contig2 = line.split("\t")[1]
  pair = [contig1, contig2]
  pair.sort()
  combo = pair[0] + "-" + pair[1]
  if combo not in pairs:
    pairs.append(combo)
    print(line, file = BLASTn_nodup_file)
  else:
    pairs.remove(combo)

BLASTn_nodup_file.close()