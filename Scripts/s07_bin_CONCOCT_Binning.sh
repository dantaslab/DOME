#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_CONCOCT_Binning.sh
# Description  : Run binning script of CONCOCT
# Usage        : sbatch s07_bin_CONCOCT_Binning.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 05/09/2022
# Last Modified: 05/09/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=concoct_binning
#SBATCH --cpus-per-task=8
#SBATCH --array=98
#SBATCH --mem=8G
#SBATCH --output=slurm_out/bin_concoct_binning/x_binning_%a.out
#SBATCH --error=slurm_out/bin_concoct_binning/y_binning_%a.err

module load miniconda3
source activate /opt/apps/labs/gdlab/envs/concoct/1.1.0

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
indir="${basedir}/d07_bin/bowtie2"
outdir="${basedir}/d07_bin/concoct"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

set -x

concoct --threads 8 --composition_file ${outdir}/${sample}/${sample}_contigs_10K.fa --coverage_file ${outdir}/${sample}/${sample}_coverage_table.tsv -b ${outdir}/${sample}/${sample}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi