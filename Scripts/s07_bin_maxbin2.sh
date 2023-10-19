#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_maxbin2.sh
# Description  : Run MaxBin on metagenomic assemblies.
# Usage        : sbatch s07_bin_maxbin2.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 05/03/2022
# Last Modified: 05/03/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=maxbin
#SBATCH --cpus-per-task=1
#SBATCH --array=108
#SBATCH --mem=8G
#SBATCH --output=slurm_out/bin_maxbin2/x_maxbin2_%a.out
#SBATCH --error=slurm_out/bin_maxbin2/y_maxbin2_%a.err

module load fraggenescan/1.31
module load bowtie2
module load hmmer
module load idba-full/1.1.1
module load maxbin/2.2.7

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
indir="${basedir}/d06_megahit"
indir2="${basedir}/d04_subjectreads"
outdir="${basedir}/d07_bin/maxbin2"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

#make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/${sample}

set -x

perl /opt/apps/labs/gdlab/software/maxbin/2.2.7/run_MaxBin.pl -min_contig_length 1500 -thread 8 -out ${outdir}/${sample}/${sample} -contig ${indir}/${sample}/final.contigs.modhead.fa -reads ${indir2}/${sample}_R1.fastq -reads2 ${indir2}/${sample}_R2.fastq

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi