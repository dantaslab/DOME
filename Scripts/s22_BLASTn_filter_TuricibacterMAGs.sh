#!/bin/bash
#===============================================================================
# File Name    : s22_BLASTn_filter_TuricibacterMAGs.sh
# Description  : This script will filter the BLASTn outputs by size
# Usage        : sbatch s22_BLASTn_filter_TuricibacterMAGs.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTn_filter_Turicibacter
#SBATCH --array=1-61%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_out/blastn_filterMAGs/x_blastfilt_Turicibacter_%a.out
#SBATCH --error=slurm_out/blastn_filterMAGs/y_blastfilt_Turicibacter_%a.err

basedir="$PWD"
indir="${basedir}/d22_BLASTn_MAGs/Turicibacter/BLASTn"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v7_TuricibacterMAGs.txt`

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