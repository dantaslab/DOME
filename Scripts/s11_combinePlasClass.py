#!/usr/bin/env python3

"""Parse through the plasclass outputs and summarize


Usage: python3 s11_combinePlasClass.py

Args:

"""

import sys
import os.path

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
	print(__doc__)
	exit()

summary_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d11_plasclass/CombinedOutput.txt"
summary_file = open(summary_path, "w")

mapping = open("/scratch/gdlab/bmahmud/DOME/shotgun/Mapping_v4_allfc_mod.txt", "r")

for sample in mapping:
	plasclass_output_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d11_plasclass/" + sample.strip() + "/PlasClass_probabilities.txt"
	if os.path.exists(plasclass_output_path):
		for line in open(plasclass_output_path, "r"):
			line = line.strip()
			print(line, file = summary_file)
	else:
		print(plasclass_output_path + " does not exist")

summary_file.close()
mapping.close()