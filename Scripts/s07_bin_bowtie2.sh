#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_bowtie2.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s07_bin_bowtie2.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Created On   : Tue Jul  9 16:13:01 CDT 2019
# Last Modified: Tue Jul  9 16:13:01 CDT 2019
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtie2
#SBATCH --array=2-283%50
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/bin_bowtie2/x_bowtie_%a.out
#SBATCH --error=slurm_out/bin_bowtie2/y_bowtie_%a.err

module load bowtie2/2.3.5-python-3.6.5
module load samtools/1.9

# Basedir
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

index="${basedir}/d07_bin/bowtieIndexes/${sample}"
indir="${basedir}/d04_subjectreads"
outdir="${basedir}/d07_bin/bowtie2/${sample}"

#make output directory
mkdir -p ${outdir}

set -x

bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${index} \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 ${indir}/${sample}_R1.fastq \
    -2 ${indir}/${sample}_R2.fastq \

time samtools view -bS ${outdir}/${sample}_mappedreads.sam > ${outdir}/${sample}_mappedreads.bam
rm ${outdir}/${sample}_mappedreads.sam
time samtools sort ${outdir}/${sample}_mappedreads.bam -o ${outdir}/${sample}_mappedreads_sorted.bam
rm ${outdir}/${sample}_mappedreads.bam
time samtools index ${outdir}/${sample}_mappedreads_sorted.bam > ${outdir}/${sample}_mappedreads_sorted.bam.bai

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi