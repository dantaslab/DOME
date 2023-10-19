#!/usr/bin/env python3

"""Merge extracted metaphlan output


Usage: python3 s23_metaphlanMerge.py <samples list> <output dir>

Args:
	samples list: A txt file (mapping file) with all samples
  output dir: directory where the extracted profiles are stored and where the merged file will be saved (ends in "/")
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 3:
  print(__doc__)
  exit()

mapping = sys.argv[1]
outdir = sys.argv[2]

samples = []

for sample in open(mapping, "r"):
  samples.append(sample.strip())

out_path = outdir + "Merged_profile_family.txt"
out_file = open(out_path, "w")

print("Sample\tFamily\tAbundance", file = out_file)

for sample in samples:
  profile = outdir + sample + "_profile_family.txt"
  for line in open(profile, "r"):
    if line.startswith("family"):
      continue
    else:
      print(sample + "\t" + line.strip(), file = out_file)

out_file.close()