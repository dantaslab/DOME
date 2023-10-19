#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_CONCOCT_Extract.sh
# Description  : Run extract script of CONCOCT
# Usage        : sbatch s07_bin_CONCOCT_Extract.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 05/09/2022
# Last Modified: 05/09/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=concoct_extract
#SBATCH --array=2-283%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm_out/bin_concoct_extract/x_extract_%a.out
#SBATCH --error=slurm_out/bin_concoct_extract/y_extract_%a.err

module load biopython/1.65
module load python/2.7.15
module load py-pandas/0.24.1-python-2.7.15
module load py-numpy/1.16.2-python-2.7.15

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
indir="${basedir}//d06_megahit"
outdir="${basedir}/d07_bin/concoct"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

mkdir ${outdir}/${sample}/fasta_bins

set -x

python CONCOCT_scripts/extract_fasta_bins.py ${indir}/${sample}/final.contigs.modhead.fa ${outdir}/${sample}/${sample}_clustering_merged.csv --output_path ${outdir}/${sample}/fasta_bins

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi