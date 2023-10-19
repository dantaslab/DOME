#!/bin/bash
#===============================================================================
# File Name    : s28_fastani_MQHQMAGs.sh
# Description  : Fast all vs. all ANI calculation of many genomes
# Usage        : sbatch s28_fastani_MQHQMAGs.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.2
# Modified     : Wed Oct 12 13:30:06 CDT 2022
# Created      : Wed Oct 12 13:30:06 CDT 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=fastani
#SBATCH --mem=80G
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm_out/fastani_MQHQ/x_fastani_%A.out
#SBATCH --error=slurm_out/fastani_MQHQ/y_fastani_%A.err

eval $( spack load --sh fastani )

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/allHQandMQ"
outdir="${basedir}/d28_fastani_HQMQ"

# Make output directory
mkdir -p ${outdir}

# Find assemblies and write full path to a file
# for an all v all analysis, that file will be ref and query
find ${indir} -type f -name "*.fa" >> ${outdir}/MAG_paths.txt

set -x
time fastANI --ql ${outdir}/MAG_paths.txt \
            --rl ${outdir}/MAG_paths.txt \
            -o ${outdir}/fastani_out.txt \
            -t ${SLURM_CPUS_PER_TASK}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occurred!"
  exit $RC
fi