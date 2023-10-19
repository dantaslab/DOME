#!/bin/bash

#===============================================================================
# File Name    : s00_catrawreads.sh
# Description  : This script will cat sequencing files for each subject
# Usage        : sbatch s00_catrawreads.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2021-12-20
# Last Modified: 2020-12-20
#===============================================================================
# Submission script for HTCF

#SBATCH --job-name=hacking
#SBATCH --array=1-712%100
#SBATCH --mem=8G
#SBATCH --output=slurm_out/catrawseq/x_catrawseq_%a.out
#SBATCH --error=slurm_out/catrawseq/y_catrawseq_%a.err

#store the base directory
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
#store the output path
outdir="${basedir}/d00_rawreads"

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v1_allfc.txt`

set -x

cat ${basedir}/seq*/${sample}_R1.fastq > ${outdir}/${sample}_R1.fastq
cat ${basedir}/seq*/${sample}_R2.fastq > ${outdir}/${sample}_R2.fastq

RC=$?

set +x

#output if job was successful
if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured!"
    exit $RC
fi