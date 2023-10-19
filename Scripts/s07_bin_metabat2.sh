#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_metabat2.sh
# Description  : Run MetaBat on metagenomic assemblies.
# Usage        : sbatch s07_bin_metabat2.sh
# Author       : Bejan Mahmud (based on Robert Thaenert)
# Version      : 1.0
# Created On   : 05/03/2022
# Last Modified: 05/03/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=metabat2
#SBATCH --cpus-per-task=8
#SBATCH --array=2-283%50
#SBATCH --mem=8G
#SBATCH --output=slurm_out/bin_metabat2/x_metabat2_%a.out
#SBATCH --error=slurm_out/bin_metabat2/y_metabat2_%a.err

module load metabat/2.11.2

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

indir="${basedir}/d06_megahit"
indir2="${basedir}/d07_bin/bowtie2"
outdir="${basedir}/d07_bin/metabat2"

#make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/${sample}

set -x

jgi_summarize_bam_contig_depths --outputDepth ${outdir}/${sample}/final.contigs.modhead.fa.depth.txt --pairedContigs ${outdir}/${sample}/final.contigs.modhead.fa.paired.txt --minContigLength 1000 --minContigDepth 1  ${indir2}/${sample}/${sample}_mappedreads_sorted.bam

metabat2 --minContig 1500 --unbinned --inFile ${indir}/${sample}/final.contigs.modhead.fa --outFile ${outdir}/${sample}/bin --abdFile ${outdir}/${sample}/final.contigs.modhead.fa.depth.txt

#runMetaBat.sh --minContig 1500 --unbinned ${indir}/${sample}/final.contigs.fa ${indir2}/${sample}/${sample}_mappedreads_sorted.bam

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi