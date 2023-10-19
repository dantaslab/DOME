#!/bin/bash
#===============================================================================
#
# File Name    : s05_metaspades.sh
# Description  : This script will run spades in parallel
# Usage        : sbatch s05_metaspades.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.2
# Created On   : 2018_01_21
# Modified On  : Wed Aug 28 09:14:08 CDT 2019
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=metaspades
#SBATCH --array=1-284%50
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --output=slurm_out/metaspades/x_metaspades_%a.out
#SBATCH --error=slurm_out/metaspades/y_metaspades_%a.err

module load spades/3.14.0-python-2.7.15

# Basedir
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d04_subjectreads"
outdir="${basedir}/d05_metaspades"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v3_subjects.txt`

# modifiied to utilize /tmp/ on the node

set -x
time spades.py --meta -t ${SLURM_CPUS_PER_TASK} \
    --tmp-dir ${basedir}/tmp/metaspades/${sample}_01_mSpades_tmp \
    --pe1-1 ${indir}/${sample}_R1.fastq \
    --pe1-2 ${indir}/${sample}_R2.fastq \
    --pe1-s ${indir}/${sample}_singletons.fastq \
    -o ${outdir}/${sample}/

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
