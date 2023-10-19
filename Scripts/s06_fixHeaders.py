#!/usr/bin/env python3

"""Fix the megahit assembly fasta headers


Usage: python3 s06_fixHeaders.py

Args:
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
	print(__doc__)
	exit()

mapping_file = open("/scratch/gdlab/bmahmud/DOME/shotgun/Mapping_v5_subjects_mod.txt", "r")

for sample in mapping_file:

	sample = sample.strip()

	input_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit/" + sample + "/final.contigs.fa"
	input_file = open(input_path, "r")

	output_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit/" + sample + "/final.contigs.modhead.fa"
	output_file = open(output_path, "w")

	for line in input_file:
		line = line.strip()
		if line.startswith(">"):
			line = line.split(" ")[0]
			line = ">" + sample + "_" + line[1:]
			print(line, file = output_file)
		else:
			print(line, file = output_file)

	input_file.close()
	output_file.close()

mapping_file.close()