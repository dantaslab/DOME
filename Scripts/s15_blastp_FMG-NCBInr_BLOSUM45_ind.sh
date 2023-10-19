#!/bin/bash
#===============================================================================
# File Name    : s15_blastp_FMG-NCBInr_BLOSUM45_ind.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s15_blastp_FMG-NCBInr_BLOSUM45_ind.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=blastp
#SBATCH --array=1-2049%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=slurm_out/blastp_NCBInr/x_blastp_NCBInr_BLOSUM45_%a.out
#SBATCH --error=slurm_out/blastp_NCBInr/y_blastp_NCBInr_BLOSUM45_%a.err

#module load blast-plus/2.11.0-python-3.6.5
eval $( spack load --sh blast-plus )

basedir="$PWD"
indir="${basedir}/d19_needle/FMGhits"
outdir="${basedir}/d15_BLAST/FMG-NCBInr"

query=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d19_needle/220806_NCBInr_allFMGhits.txt`

set -x

time blastp \
  -query ${indir}/${query}.fasta \
  -db nr \
  -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid qcovs pident evalue" \
  -matrix BLOSUM45 \
  -out ${outdir}/${query}.txt \
  -remote

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi