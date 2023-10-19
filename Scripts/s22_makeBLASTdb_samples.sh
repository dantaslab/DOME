#!/bin/bash
#===============================================================================
# File Name    : s22_makeBLASTdb_samples.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s22_makeBLASTdb_samples.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTdb
#SBATCH --array=2-712%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=slurm_out/blastdb_samples/x_blastdb_%a.out
#SBATCH --error=slurm_out/blastdb_samples/y_blastdb_%a.err

#module load blast-plus/2.11.0-python-3.6.5
eval $( spack load --sh blast-plus@2.12.0 )

basedir="$PWD"
indir="${basedir}/d06_megahit_samples"
outdir="${basedir}/d22_BLASTn_samples/db"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

#make output directory
mkdir -p ${outdir}/${sample}

set -x

time makeblastdb -in ${indir}/${sample}/final.contigs.modhead.2line.fa \
        -parse_seqids \
        -blastdb_version 5 \
        -title "${sample}_BLASTdb" \
        -dbtype nucl \
        -out ${outdir}/${sample}/${sample}

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi