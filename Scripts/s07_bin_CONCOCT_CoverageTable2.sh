#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_CONCOCT_CoverageTable.sh
# Description  : Run converge table script of CONCOCT
# Usage        : sbatch s07_bin_CONCOCT_CoverageTable.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 05/04/2022
# Last Modified: 05/04/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=concoct_coverage
#SBATCH --cpus-per-task=8
#SBATCH --array=1
#SBATCH --mem=16G
#SBATCH --output=slurm_out/bin_concoct_coverage/x_coverage2_%a.out
#SBATCH --error=slurm_out/bin_concoct_coverage/y_coverage2_%a.err

module load concoct
module load python/2.7.15
module load py-pandas/0.24.1-python-2.7.15
module load py-numpy/1.16.2-python-2.7.15
module load samtools/1.9
module load py-pybedtools/0.8.0-python-2.7.15
#module load numpy/1.22

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
indir="${basedir}/d07_bin/bowtie2"
outdir="${basedir}/d07_bin/concoct"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

set -x

python CONCOCT_scripts/concoct_coverage_table.py ${outdir}/${sample}/${sample}_contigs_10K.bed ${indir}/${sample}/${sample}_mappedreads_sorted.bam > ${outdir}/${sample}/${sample}_2_coverage_table.tsv

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi