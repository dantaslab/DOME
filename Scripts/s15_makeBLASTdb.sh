#!/bin/bash
#===============================================================================
# File Name    : s15_makeBLASTdb.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s15_makeBLASTdb.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTdb
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=slurm_out/blastdb/x_blastdb_%a.out
#SBATCH --error=slurm_out/blastdb/y_blastdb_%a.err

module load blast-plus/2.11.0-python-3.6.5
#eval $( spack load --sh blast-plus )

basedir="$PWD"
indir="${basedir}/BLASTdb"

set -x

time makeblastdb -in ${indir}/e01_AR-databases_protseqs_FixTrimIDs.fasta \
        -parse_seqids \
        -blastdb_version 5 \
        -title "AMRProt_BLASTdb" \
        -dbtype prot

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi