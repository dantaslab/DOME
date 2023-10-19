#!/bin/bash
#===============================================================================
# File Name    : s11_plasclass.sh
# Description  : This script will run find plasmids in short-read assemblies
# Usage        : sbatch s11_plasclass.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Created On   : Tue Jul  7 21:45:13 CDT 2020
# Modified On  : Tue Jul  7 21:45:13 CDT 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=plasclass
#SBATCH --array=1-712%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=slurm_out/plasclass/x_plasclass_%a.out
#SBATCH --error=slurm_out/plasclass/y_plasclass_%a.err

eval $( spack load --sh py-plasclass )

# Basedir
basedir="$PWD"
indir="${basedir}/d06_megahit_samples"
outdir="${basedir}/d11_plasclass"

# read in sample name
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

# Make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/${sample}

set -x

time classify_fasta.py \
		-f ${indir}/${sample}/final.contigs.modhead.fa \
		-o ${outdir}/${sample}/PlasClass_probabilities.txt \
		-p ${SLURM_CPUS_PER_TASK}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi