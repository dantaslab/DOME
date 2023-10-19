#!/usr/bin/env python3

"""Fix the megahit assembly fasta headers


Usage: python3 s06_getContigLen_samples.py

Args:
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
	print(__doc__)
	exit()

contig_len_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit_samples/ContigLengths.txt"
contig_len_file = open(contig_len_path, "w")

mapping_file = open("/scratch/gdlab/bmahmud/DOME/shotgun/Mapping_v4_allfc_mod.txt", "r")

for sample in mapping_file:

	sample = sample.strip()

	input_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit_samples/" + sample + "/final.contigs.fa"
	input_file = open(input_path, "r")

	for line in input_file:
		if line.startswith(">"):
			line = line.strip()
			contig = line.split(" ")[0]
			contig = sample + "_" + contig[1:]
			length = line.split(" ")[3]
			length = length.split("=")[1]
			print(contig + "\t" + length, file = contig_len_file)

	input_file.close()

mapping_file.close()
contig_len_file.close()