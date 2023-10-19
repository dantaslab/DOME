#!/bin/bash
#===============================================================================
# File Name    : s22_makeBLASTdb_HQMQMAGs.sh
# Description  : This script will make a BLAST database from the input set of seq files
# Usage        : sbatch s22_makeBLASTdb_HQMQMAGs.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTdb
#SBATCH --array=1-4884%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=slurm_out/blastdb_HQMQ/x_blastdb_%a.out
#SBATCH --error=slurm_out/blastdb_HQMQ/y_blastdb_%a.err

#module load blast-plus/2.11.0-python-3.6.5
#eval $( spack load --sh blast-plus@2.12.0 )
eval $( spack load --sh /u7ssbm4 )

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/allHQandMQ"
outdir="${basedir}/d22_BLASTn_HQMQ/db"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v9_allHQandMQ-MAGs.txt`

#make output directory
mkdir -p ${outdir}/${sample}

set -x

time makeblastdb -in ${indir}/${sample}.fa \
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