#!/bin/bash

#===============================================================================
# File Name    : s17_3_bowtieAllign.sh
# Description  : This script will use quast to do quality control on assemblies
# Usage        : sbatch s17_3_bowtieAllign.sh
# Author       : Sanjam Sawhney
# Version      : 1.0
# Created On   : Wed Jul 24 15:46:45 CDT 2019
# Last Modified: Tue Aug 13 19:24:00 CDT 2019
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=4
#SBATCH --job-name=bowtie_align_samtools
#SBATCH --mem=30G
#SBATCH --array=136
#SBATCH --output=slurm_out/inStrain_bowtieAlign/x_inStrain_bowtieAlign_%a.out
#SBATCH --error=slurm_out/inStrain_bowtieAlign/y_inStrain_bowtieAlign_%a.err

#Load dependencies
eval $( spack load --sh bowtie2@2.4.2 )
#eval $( spack load --sh samtools@1.14 )
eval $( spack load --sh /6p5wlkk )

basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d17_inStrain/BowtieBuild"
filedir="${basedir}/d03_clean_reads"
outdir="${basedir}/d17_inStrain/BowtieAlign"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

#make output directory
mkdir -p ${outdir}

set -x

#Note mothers have a hyphen between ID and sample instead of an underscore. Need to revise bowtie2 -1 and -2 syntax for array jobs 177-214
time bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${indir}/MAGs_cat \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 ${filedir}/${sample}_fwd_clean.fastq \
    -2 ${filedir}/${sample}_rev_clean.fastq \

time samtools view -bS ${outdir}/${sample}_mappedreads.sam > ${outdir}/${sample}_mappedreads.bam
rm ${outdir}/${sample}_mappedreads.sam
time samtools sort ${outdir}/${sample}_mappedreads.bam -o ${outdir}/${sample}_mappedreads_sorted.bam
rm ${outdir}/${sample}_mappedreads.bam

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi