#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_dasTool.sh
# Description  : Run dasTool
# Usage        : sbatch s07_bin_dasTool.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 5/10/2022
# Last Modified: 5/10/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=dasTool
#SBATCH --array=2-283%50
#SBATCH --mem=4G
#SBATCH --output=slurm_out/bin_dasTool/x_dasTool_%a.out
#SBATCH --error=slurm_out/bin_dasTool/y_dasTool_%a.err

module load das_tool/1.1.4

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
outdir="${basedir}/d07_bin/dastool"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

#make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/${sample}

set -x

DAS_Tool -i ${basedir}/d07_bin/concoct/${sample}/fasta_bins/${sample}_step5.txt,${basedir}/d07_bin/maxbin2/${sample}/${sample}_step5.txt,${basedir}/d07_bin/metabat2/${sample}/${sample}_step5.txt -l concoct,maxbin,metabat -c ${basedir}/d06_megahit/${sample}/final.contigs.modhead.fa -o ${outdir}/${sample}/${sample} --search_engine diamond --write_bins --threads 8 --write_unbinned --score_threshold 0.1

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi