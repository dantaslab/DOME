#!/bin/bash
#===============================================================================
#
# File Name    : s17_5_MAGs_inStrainCompare.sh
# Description  : This script will run instrain compare of sam files
# Usage        : sbatch s17_5_MAGs_inStrainProfile.sh
# Author       : Kimberley Sukhum, kvsukhum@wustl.edu
# Version      : 1.0
# Created On   : 2021-08-09
# Last Modified: 2021-09-09
# Modified By  : Sanjam Sawhney, ssawhney@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=instrainCompare
#SBATCH --array=8
#SBATCH --mem=100G
#SBATCH --output=slurm_out/inStrain_Compare_MAGs/x_inStrain_Compare_%a.out
#SBATCH --error=slurm_out/inStrain_Compare_MAGs/y_inStrain_Compare_%a.err
#SBATCH --cpus-per-task=8

#Load dependencies
eval $( spack load --sh py-instrain@1.5.7 )
#module load python
#module load samtools/1.9
#. /opt/apps/labs/gdlab/envs/instrain/1.5.4/bin/activate

# samplelist.txt contains a list of all file names
basedir="$PWD"

mag=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v6_inStrainMAGs.txt`

indir="${basedir}/d17_inStrain_MAGs/${mag}/Profile"
outdir="${basedir}/d17_inStrain_MAGs/${mag}/Compare"

#make output directory
mkdir -p ${outdir}

set -x

time inStrain compare -i ${indir}/*.IS/ \
    -o ${outdir}/${mag}.IS.COMPARE \
    -p ${SLURM_CPUS_PER_TASK} \
    -s ${basedir}/d17_inStrain_MAGs/${mag}/MAGs/MAGs_cat.stb \
    --database_mode

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi