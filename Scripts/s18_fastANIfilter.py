#!/usr/bin/env python3

"""Filter the pairwise fastANI output


Usage: python3 fastANIfilter.py <fastANI output>

Args:
	fastANI output: Tab-separated txt file contaning the output from fastANI
"""

import sys
import pandas as pd

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 2:
  print(__doc__)
  exit()

fastANI_output = sys.argv[1]

samples = []

for line in open(fastANI_output):
    line = line.strip().split("\t")
    sample1 = line[0]
    sample2 = line[1]
    if sample1 not in samples:
        samples.append(sample1)
    if sample2 not in samples:
        samples.append(sample2)

samples_num = {}

count = 0

for sample in samples:
    count += 1
    samples_num[sample] = count

samples_rev_num = {}

for key, value in samples_num.items():
    samples_rev_num[value] = key

fastANI_num = []

for line in open(fastANI_output):
    line = line.strip().split("\t")
    sample1 = line[0]
    sample2 = line[1]
    val = line[2]
    if sample1 != sample2:
        sample1 = samples_num[sample1]
        sample2 = samples_num[sample2]
        combo = str(min(sample1, sample2)) + "-" + str(max(sample1, sample2))
        fastANI_num.append([combo, float(val)])
    else:
        continue

fastANI_num_df = pd.DataFrame(fastANI_num, columns = ['combo', 'ANI'])

fastANI_num_df = fastANI_num_df.groupby('combo').mean().reset_index()

fastANI_filt = []

for item in fastANI_num_df.values.tolist():
    combo = item[0]
    value = item[1]
    sample1 = combo.split("-")[0]
    sample1 = samples_rev_num[int(sample1)]
    sample2 = combo.split("-")[1]
    sample2 = samples_rev_num[int(sample2)]
    entry = sample1 + "\t" + sample2 + "\t" + str(value)
    fastANI_filt.append(entry)

fastANI_filt_path = fastANI_output.split(".txt")[0] + "_filtered.txt"
fastANI_filt_file = open(fastANI_filt_path, "w")

for line in fastANI_filt:
    print(line.strip(), file = fastANI_filt_file)
    
fastANI_filt_file.close()