#!/bin/bash

#===============================================================================
# Description  : This script will run bbtools repair (bbmap) to repair deconseq reads
# Usage        : sbatch s03_bbmap.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2020_04_28
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bbtools
#SBATCH --array=1-712%50
#SBATCH --mem=16G
#SBATCH --output=slurm_out/bbmap/x_bbtools_large_%a.out
#SBATCH --error=slurm_out/bbmap/y_bbtools_large_%a.err

module load bbmap

basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d02_deconseq"
outdir="${basedir}/d03_clean_reads"

#make the output directory
mkdir -p ${outdir}

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v1_allfc.txt`

set -x

repair.sh ow=t in=${indir}/${sample}_fwd_clean.fq in2=${indir}/${sample}_rev_clean.fq out=${outdir}/${sample}_fwd_clean.fastq out2=${outdir}/${sample}_rev_clean.fastq outs=${outdir}/${sample}_singletons_CLEAN.fastq repair=t

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