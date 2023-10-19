#!/bin/bash
#===============================================================================
#
# File Name    : s21_simka_k21.sh
# Description  : This script will run simka on all clean forward reads with kmer = 21
# Usage        : sbatch s17_simka_k21.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2021-09-26
# Last Modified: 2021-09-26
# Modified By  : Bejan Mahmud
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=simka
#SBATCH --mem=200G
#SBATCH --output=slurm_out/simka_k21/x_simka_%a.out
#SBATCH --error=slurm_out/simka_k21/y_simka_%a.err
#SBATCH --cpus-per-task=8

#Load dependencies


# samplelist.txt contains a list of all file names
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d21_Simka/k21"

#make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/tmp

set -x

time /scratch/gdlab/bmahmud/simka-v1.5.3-bin-Linux/bin/simka -in ${basedir}/s21_simka_input.txt \
    -out ${outdir} \
    -out-tmp ${outdir}/tmp \
    -nb-cores ${SLURM_CPUS_PER_TASK} \
    -max-reads 0 \
    -abundance-min 2 \
    -kmer-size 21

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi