#!/bin/bash
#===============================================================================
# File Name    : s15_blastp.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s15_blastp.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=blastp
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --output=slurm_out/blastp/x_blastp_DOMEFMG.out
#SBATCH --error=slurm_out/blastp/y_blastp_DOMEFMG.err

module load blast-plus/2.11.0-python-3.6.5
#eval $( spack load --sh blast-plus )

basedir="$PWD"
dbdir="${basedir}/BLASTdb"
indir="/scratch/gdlab/bmahmud/parfums"
outdir="${basedir}/d15_BLAST"

mkdir -p ${outdir}

set -x

time blastp \
  -query ${indir}/fmg_DOME_All_filtered_contigs_r05_merged_seqs_filtered.fasta \
  -db ${dbdir}/e01_AR-databases_protseqs_FixTrimIDs.fasta \
  -max_target_seqs 5 \
  -outfmt "6 qseqid sseqid qcovs pident evalue" \
  -out ${outdir}/DOME_FMG-NCBI_CARD-blastp.txt

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi