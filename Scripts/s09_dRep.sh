#!/bin/bash
#===============================================================================
# File Name    : s09_dRep.sh
# Description  : Dereplicate a set of highly similar genomes
# Usage        : sbatch s09_dRep.sh
# Author       : Luke Diorio-Toth
# Version      : 1.0
# Created On   : Wed Jan 12 09:24:34 CST 2022
# Last Modified: Wed Jan 12 09:24:34 CST 2022
#===============================================================================
# Submission script for HTCF
#SBATCH --job-name=dRep
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --output=slurm_out/dRep/x_dRep_%a.out
#SBATCH --error=slurm_out/dRep/y_dRep_%a.err

#load dependencies
module load prodigal
module load mash
module load MUMmer
module load fastANI
module load ANIcalculator
module load centrifuge
module load openmpi/3.1.3-python-2.7.15
module load hmmer/3.1b1-python-2.7.15
module load pplacer/1.1.alpha19

. /opt/apps/labs/gdlab/envs/drep/3.2.2/bin/activate

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/allbins_noUnbin"
outdir="${basedir}/d09_dRep/dRep_99"

mkdir -p ${outdir}

set -x
time dRep dereplicate \
    -g ${indir}/*.fa \
    -d \
    -p ${SLURM_CPUS_PER_TASK} \
    -sa 0.99 \
    -nc 0.25 \
    ${outdir}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi