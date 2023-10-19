#!/bin/bash

#===============================================================================
# File Name    : s02_deconseq.sh
# Description  : This script will remove human and cow read contamination from fastq reads in parallel
# Usage        : sbatch s02_deconseq.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2020-04-14
# Last Modified: 2020-04-14
#===============================================================================
# Submission script for HTCF

#SBATCH --job-name=deconseq
#SBATCH --array=1-712%50
#SBATCH --mem=32G
#SBATCH --output=slurm_out/decon/x_deconseq_%a.out
#SBATCH --error=slurm_out/decon/y_deconseq_%a.err

#load the module
module purge
module load deconseq/0.4.3-chr38

#store the base directory
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
#store the input path
indir="${basedir}/d01_trimReads"
#store the output path
outdir="${basedir}/d02_deconseq"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v1_allfc.txt`

#deconseq command

set -x

time deconseq.pl \
     -f ${indir}/${sample}_FW_trimmed_paired.fastq \
     -out_dir ${outdir} \
     -id ${sample}_fwd \
     -dbs hsref38,cow

RC_FW=$?

time deconseq.pl \
     -f ${indir}/${sample}_RV_trimmed_paired.fastq \
     -out_dir ${outdir} \
     -id ${sample}_rev \
     -dbs hsref38,cow

RC_RV=$?

set +x

RC=$((RC_FW + RC_RV))

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured!"
    exit $RC
fi