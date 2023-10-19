#!/usr/bin/env python3

"""Make 2line foramts of the assembly fasta files


Usage: python3 s06_make2line_samples.py

Args:
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 1:
    print(__doc__)
    exit()

mapping_file = open("/scratch/gdlab/bmahmud/DOME/shotgun/Mapping_v4_allfc_mod.txt", "r")

#The for loop below will take in the assemly.fasta file for each sample and convert it to the 2line format
for sample in mapping_file:
    sample = sample.strip()
    fasta_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit_samples/" + sample + "/final.contigs.modhead.fa"
    two_line = []
    seq = ""
    with open(fasta_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                if seq != "":
                    two_line.append(seq)
                    seq = ""
                two_line.append(line.strip())
            else:
                seq = seq + line.strip()
        two_line.append(seq)
    fasta_2line_path = "/scratch/gdlab/bmahmud/DOME/shotgun/d06_megahit_samples/" + sample + "/final.contigs.modhead.2line.fa"
    fasta_2line = open(fasta_2line_path, "w")
    for line in two_line:
        print(line, file = fasta_2line)
    fasta_2line.close()