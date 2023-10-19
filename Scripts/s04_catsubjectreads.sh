#!/bin/bash

#===============================================================================
# File Name    : s01_catsubjectreads.sh
# Description  : This script will cat sequencing files for each subject
# Usage        : sbatch s01_catsubjectreads.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2021-12-20
# Last Modified: 2020-12-20
#===============================================================================
# Submission script for HTCF

#SBATCH --job-name=hacking
#SBATCH --array=1-284%50
#SBATCH --mem=16G
#SBATCH --output=slurm_out/catsubseqs/x_catsubreads_%a.out
#SBATCH --error=slurm_out/catsubseqs/y_catsubreads_%a.err

#store the base directory
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
#store the input path
indir="${basedir}/d03_clean_reads"
#store the output path
outdir="${basedir}/d04_subjectreads"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v3_subjects.txt`

set -x

cat ${indir}/*${sample}*fwd* > ${outdir}/${sample}_R1.fastq
cat ${indir}/*${sample}*rev* > ${outdir}/${sample}_R2.fastq
cat ${indir}/*${sample}*singletons* > ${outdir}/${sample}_singletons.fastq

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