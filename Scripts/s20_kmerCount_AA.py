#!/usr/bin/env python3

"""Count kmer occurences in a bakta faa file


Usage: python3 s20_kmerCount_AA.py <bakta faa file> <output file> <sample>

Args:
	bakta faa file: .faa file from bakta which contains all ORFs. Note: each protein is presented in a 2-line format
  output file: path to the file where counts for each kmer will be recorded
  sample: sample ID
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 4:
  print(__doc__)
  exit()

faa_file = sys.argv[1]
output_path = sys.argv[2]
sample = sys.argv[3]

kmer_sizes = [8, 9, 10, 11, 12]

kmers = {}

input_file = open(faa_file, "r")

for line in input_file:
    if not line.startswith(">"):
        line = line.strip()
        for size in kmer_sizes:
            if len(line) > size - 1:
                n = 0
                while (n + size) <= len(line):
                    kmer = line[0+n:size+n]
                    if kmer in kmers.keys():
                        kmers[kmer] += 1
                    else:
                        kmers[kmer] = 1
                    n += 1

input_file.close()

output_file = open(output_path, "w")

print("Sample\tkmer\tCount", file = output_file)

for key, value in kmers.items():
  print(sample + "\t" + key + "\t" + str(value), file = output_file)

output_file.close()