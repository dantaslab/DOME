#!/bin/bash
#===============================================================================
# File Name    : s15_blastp_FMG-NCBInr_BLOSUM45.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s15_blastp_FMG-NCBInr_BLOSUM45.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=blastp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=slurm_out/blastp_NCBInr/x_blastp_NCBInr_BLOSUM45.out
#SBATCH --error=slurm_out/blastp_NCBInr/y_blastp_NCBInr_BLOSUM45.err

#module load blast-plus/2.11.0-python-3.6.5
eval $( spack load --sh blast-plus )

basedir="$PWD"
outdir="${basedir}/d15_BLAST"

set -x

time blastp \
  -query /scratch/gdlab/bmahmud/parfums/fmg_DOME_All_filtered_contigs_r05_merged_seqs_filtered.fasta \
  -db /ref/gdlab/data/blast_db/nr_2022-08-04/nr \
  -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid qcovs pident evalue" \
  -matrix BLOSUM45 \
  -out ${outdir}/DOME_FMG-NCBInr-BLOSUM45.txt

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi