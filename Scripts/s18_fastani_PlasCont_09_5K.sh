#!/bin/bash
#===============================================================================
# File Name    : s18_fastani_PlasCont_09_5K.sh
# Description  : This script will run fastANI
# Usage        : sbatch s18_fastani_PlasCont_09_5K.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Modified     : Wed Apr 15 14:42:30 CDT 2020
# Created      : Wed Apr 15 14:42:30 CDT 2020
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=fastani
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm_out/fastani_PlasCont_09_5K/x_fastani_%a.out
#SBATCH --error=slurm_out/fastani_PlasCont_09_5K/y_fastani_%a.err

module load fastani

basedir="$PWD"
# /scratch/gdlab/bmahmud/DOME/shotgun
outdir="${basedir}/d18_fastani_PlasCont_09_5K"
plasmids="${basedir}/PlasCont_09_5K_paths.txt"


# Make output directory
mkdir -p ${outdir}/array


set -x
time fastANI --ql ${plasmids} \
            --rl ${plasmids} \
            --fragLen 250 \
            --minFraction 0.5 \
            -o ${outdir}/array/fastani_PlasCont_09_5K_250bp_0.5min.txt \
            -t ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi