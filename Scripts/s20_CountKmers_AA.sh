#!/bin/bash
#===============================================================================
# File Name    : s20_CountKmers_AA.sh
# Description  : This script will count the kmer prevalences in a bakta faa file
# Usage        : sbatch s20_CountKmers_AA.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Sep 20 2022
# Last Modified: Tue Sep 20 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=kmerCount
#SBATCH --array=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm_out/kmerCount/x_kmerCount_%a.out
#SBATCH --error=slurm_out/kmerCount/y_kmerCount_%a.err

basedir="$PWD"

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v1_allfc.txt`

indir="${basedir}/d25_bakta_samples"
outdir="${basedir}/d20_kmerCount_AA"

#make output directory
mkdir -p ${outdir}

set -x

python3 s20_kmerCount_AA.py ${indir}/${sample}_fwd_clean.fastq ${outdir}/${sample}_kmerCounts.txt ${sample}

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