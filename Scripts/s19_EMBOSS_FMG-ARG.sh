#!/bin/bash
#===============================================================================
# File Name    : s19_EMBOSS_FMG-ARG.sh
# Description  : This script will run the needle alignment of FMG hits with their top CARD/AMRProt matches
# Usage        : s19_EMBOSS_FMG-ARG.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=needle
#SBATCH --array=1-1992%50
#SBATCH --mem=100MB
#SBATCH --output=slurm_out/emboss_ARG/x_emboss_ARG_%a.out
#SBATCH --error=slurm_out/emboss_ARG/y_emboss_ARG_%a.err

module load emboss/6.6.0-python-3.6.5

basedir="$PWD"
argdir="${basedir}/d19_needle/ARGdb"
fmgdir="${basedir}/d19_needle/FMGhits"
outdir="${basedir}/d19_needle/FMG-ARG"

mkdir -p ${outdir}

arg=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d19_needle/220805_ARGdb_ARGdbMatches.txt`
fmg=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d19_needle/220805_ARGdb_FMGhits.txt`

set -x

time needle ${fmgdir}/${fmg}.fasta ${argdir}/${arg}.fasta -gapopen 10 -gapextend 0.5 -outfile ${outdir}/${fmg}.txt

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi