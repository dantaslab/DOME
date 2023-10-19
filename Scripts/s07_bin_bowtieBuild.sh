#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_bowtieBuild.sh
# Description  : This script will index a ref for Burrows-Wheeler alignment
# Usage        : sbatch s07_bin_bowtieBuild.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Created On   : Tue Jul  9 16:13:01 CDT 2019
# Last Modified: Tue Jul  9 16:13:01 CDT 2019
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtieBuild
#SBATCH --array=2-283%50
#SBATCH --mem=8G
#SBATCH --output=slurm_out/bin_bowtieBuild/x_bowtieBuild_%a.out
#SBATCH --error=slurm_out/bin_bowtieBuild/y_bowtieBuild_%a.err

module load bowtie2/2.3.4.1

# Basedir
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

indir="${basedir}/d06_megahit"
outdir="${basedir}/d07_bin/bowtieIndexes"

#make output directory
mkdir -p ${outdir}

bowtie2-build ${indir}/${sample}/final.contigs.modhead.fa ${outdir}/${sample}

RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured in ${sample}!"
    exit $RC
fi