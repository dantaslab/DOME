#!/bin/bash
#===============================================================================
# File Name    : s22_BLASTn_nodup.sh
# Description  : This script will filter out the duplicates (A B vs B A) from the joint BLASTn output file
# Usage        : sbatch s22_BLASTn_nodup.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTn
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm_out/blastn_nodup_HQMQ/x_blastn_BK_%A.out
#SBATCH --error=slurm_out/blastn_nodup_HQMQ/y_blastn_BK_%A.err

basedir="$PWD"
indir="${basedir}/d22_BLASTn_HQMQ"

set -x

python3 s22_nodupBLASTn_BK.py ${indir}/joint_BLASTn_allHQMQMAGs.txt ${indir}/joint_BLASTn_allHQMQMAGs_nodup_BK.txt

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi