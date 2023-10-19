#!/usr/bin/env python3

"""Parse through the platon outputs and summarize


Usage: python3 parsePlaton.py

Args:

"""

import sys
import os.path

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
	print(__doc__)
	exit()

summary_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d10_platon/AllPlasmid.txt"
summary_file = open(summary_path, "w")

mapping = open("/scratch/gdlab/bmahmud/DOME/shotgun/Mapping_v4_allfc_mod.txt", "r")

for sample in mapping:
	platon_output_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d10_platon/" + sample.strip() + "/" + sample.strip() + ".plasmid.fasta"
	if os.path.exists(platon_output_path):
		for line in open(platon_output_path, "r"):
			if line.startswith(">"):
				print(line.strip()[1:], file = summary_file)
	else:
		print(platon_output_path + " does not exist")

summary_file.close()
mapping.close()