#!/bin/bash
#===============================================================================
# File Name    : s19_EMBOSS_FMG-NCBInr.sh
# Description  : This script will run the needle alignment of FMG hits with their top NCBInr matches
# Usage        : s19_EMBOSS_FMG-NCBInr.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=needle
#SBATCH --array=1-2049%50
#SBATCH --mem=100MB
#SBATCH --output=slurm_out/emboss_NCBInr/x_emboss_NCBInr_%a.out
#SBATCH --error=slurm_out/emboss_NCBInr/y_emboss_NCBInr_%a.err

module load emboss/6.6.0-python-3.6.5

basedir="$PWD"
ncbidir="${basedir}/d15_BLAST/FMG-NCBInr_TopUniqueHits"
fmgdir="${basedir}/d19_needle/FMGhits"
outdir="${basedir}/d19_needle/FMG-NCBInr"

mkdir -p ${outdir}

ncbi=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d19_needle/220818_NCBInr_NCBIhits.txt`
fmg=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d19_needle/220818_NCBInr_FMGhits_fixedHeaders.txt`

set -x

time needle ${fmgdir}/${fmg}.fasta ${ncbidir}/${ncbi}.fasta -gapopen 10 -gapextend 0.5 -outfile ${outdir}/${fmg}.txt

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi