#!/bin/bash
#===============================================================================
# File Name    : s22_BLASTn_filter_HQMQMAGs.sh
# Description  : This script will filter the BLASTn outputs by size
# Usage        : sbatch s22_BLASTn_filter.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTn
#SBATCH --array=1-4884%300
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_out/blastn_filterHQMQ/x_blastn_%a.out
#SBATCH --error=slurm_out/blastn_filterHQMQ/y_blastn_%a.err

basedir="$PWD"
indir="${basedir}/d22_BLASTn_HQMQ/BLASTn"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v9_allHQandMQ-MAGs.txt`

set -x

python3 s22_filterBLASTn.py ${indir}/merged_${sample}.txt ${indir}/filtered_${sample}.txt

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi