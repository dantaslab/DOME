#!/usr/bin/env python3

"""Merge normalized HUMANN3 output files


Usage: python3 s24_merge_humann3.py <samples list> <output dir>

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

out_path = outdir + "Merged_pathabundance.txt"
out_file = open(out_path, "w")

print("Sample\tPathwau\tAbundance (CPM)", file = out_file)

for sample in samples:
    profile = outdir + sample + "_pathabundance-cpm.tsv"
    for line in open(profile, "r"):
        if line.startswith("#"):
            continue
        else:
            line = line.strip().split("\t")
            pathway = line[0]
            count = line[1]
            if "|g__" in pathway:
                continue
            elif "|unclassified" in pathway:
                continue
            else:
                print("\t".join([sample, pathway, count]), file = out_file)
            
out_file.close()