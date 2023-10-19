#!/bin/bash
#===============================================================================
# File Name    : s08_prokka_samples.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s08_prokka_samples.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.3
# Modified     : Thu Jun 11 08:54:19 CDT 2020
# Created      : Mon Jun 10 17:17:41 CDT 2019
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=prokka
#SBATCH --array=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/prokka_samples/x_prokka_%a.out
#SBATCH --error=slurm_out/prokka_samples/y_prokka_%a.err

module load openmpi/3.1.3-python-2.7.15
module load rnammer
module load prokka/1.14.5-python-2.7.15-java-11

basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d06_megahit_samples"
outdir="${basedir}/d08_prokka_samples"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

set -x
time prokka ${indir}/${sample}/final.contigs.fa \
            --outdir ${outdir}/${sample} \
            --prefix ${sample} \
            --mincontiglen 200 \
            --force \
            --rnammer \
            --cpus ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi