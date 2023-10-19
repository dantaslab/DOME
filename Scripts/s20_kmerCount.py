#!/usr/bin/env python3

"""Count kmer occurences in a fastq file of clean metagenomic reads


Usage: python3 s20_kmerCount.py <fastq file> <output file> <sample>

Args:
	fastq file: fastq file of clean (post trimmomatic, deconseq, bbmap) shotgun metagenomic reads
  output file: path to the file where counts for each kmer will be recorded
  sample: sample ID
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 4:
  print(__doc__)
  exit()

fastq_file = sys.argv[1]
output_path = sys.argv[2]
sample = sys.argv[3]

line_counter = 0

kmers = {}

for line in open(fastq_file, "r"):
    line_counter += 1
    if line_counter%4 == 2:
        line = line.strip()
        if len(line) > 16:
            n = 0
            while (n + 17) <= len(line):
                kmer = line[0+n:17+n]
                if "N" not in kmer:
                  if kmer in kmers.keys():
                      kmers[kmer] += 1
                  else:
                      kmers[kmer] = 1
                n += 1

output_file = open(output_path, "w")

for key, value in kmers.items():
  print(sample + "\t" + key + "\t" + str(value), file = output_file)

output_file.close()