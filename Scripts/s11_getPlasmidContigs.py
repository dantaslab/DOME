#!/usr/bin/env python3

"""This 


Usage: python3 s11_getPlasmidContigs.py

Args:

"""

import sys
import os.path

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
	print(__doc__)
	exit()

plas_contings = open("/scratch/gdlab/bmahmud/DOME/shotgun/d11_plasclass/220720_PlasContigs_09_5K.txt", "r")

for contig in plas_contings:
	contig = contig.strip()
	sample = contig.split("_")[0]

	output_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d11_plasclass/PlasmidContigs_09_5K/" + contig + ".fasta"
	output_file = open(output_path, "w")
	print(">" + contig, file = output_file)

	assembly_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit_samples/" + sample + "/final.contigs.modhead.2line.fa"
	assembly_1 = open(assembly_path, "r")
	assembly_2 = open(assembly_path, "r")
	assembly_2 = assembly_2.readlines()

	counter = 0

	for line in assembly_1:
		counter = counter + 1
		if line.strip()[1:] == contig:
			print(assembly_2[counter], file = output_file)

	assembly_1.close()
	output_file.close()

plas_contings.close()