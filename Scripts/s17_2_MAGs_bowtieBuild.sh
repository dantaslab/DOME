#!/bin/bash

#===============================================================================
# File Name    : s17_2_MAGs_bowtieBuild.sh
# Description  : This script will use quast to do quality control on assemblies
# Usage        : sbatch s17_2_MAGs_bowtieBuild.sh
# Author       : Sanjam Sawhney
# Version      : 1.0
# Created On   : Wed Jul 24 15:46:45 CDT 2019
# Last Modified: Tue Aug 13 19:24:00 CDT 2019
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=8
#SBATCH --job-name=bowtie_build
#SBATCH --array=89
#SBATCH --mem=4G
#SBATCH --output=slurm_out/inStrain_bowtieBuild_MAGs/x_inStrain_bowtieBuild_%a.out
#SBATCH --error=slurm_out/inStrain_bowtieBuild_MAGs/y_inStrain_bowtieBuild_%a.err

#Load dependencies
eval $( spack load --sh bowtie2@2.4.2 )

basedir="$PWD"

mag=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v6_inStrainMAGs.txt`

indir="${basedir}/d17_inStrain_MAGs/${mag}/MAGs"
outdir="${basedir}/d17_inStrain_MAGs/${mag}/BowtieBuild"

#make output directory
mkdir -p ${outdir}

#enter debug mode
set -x

time bowtie2-build ${indir}/MAGs_cat.fasta ${outdir}/MAGs_cat

#store exit code
RC=$?
#exit debug mode
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi