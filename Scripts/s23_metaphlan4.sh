#!/bin/bash
#===============================================================================
# File Name    : s23_metaphlan4.sh
# Description  : This script will run metaphlan4 on paired reads
# Usage        : sbatch s23_metaphlan4.sh
# Author       : Luke Diorio-Toth
# Version      : 1.0
# Created On   : Wed Jan  6 11:49:47 CST 2021
# Last Modified: Wed Jan  6 11:49:47 CST 2021
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=metaphlan4
#SBATCH --array=2-712%50
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm_out/metaphlan4/x_metaphlan_%a.out
#SBATCH --error=slurm_out/metaphlan4/y_metaphlan_%a.err

# load python and activate metaphlan3 venv
eval $( /scratch/gdlab/ldioriot/spack/bin/spack load --sh py-metaphlan )
eval $( spack load --sh python@3.9.12 )
eval $( spack load --sh bowtie2@2.4.2 )

basedir="$PWD"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d23_metaphlan4"

#make output directory
mkdir -p ${outdir}

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

set -x # Start debug mode (will send commands to outfile)

time metaphlan \
	${indir}/${sample}_fwd_clean.fastq,${indir}/${sample}_rev_clean.fastq \
	--bowtie2out ${outdir}/${sample}.bowtie2.bz2 \
	-o ${outdir}/${sample}_profile.txt \
	--nproc ${SLURM_CPUS_PER_TASK} \
	--input_type fastq

RC=$? # Save error code for command

set +x # End debug mode (will stop sending commands to outfile)

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi