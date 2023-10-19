#!/bin/bash
#===============================================================================
# File Name    : s12_shortbred_reads.sh
# Description  : This script will profile sequence data and isolate genomes. 
# Usage        : sbatch s12_shortbred_reads.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2020-05-31
# Last Modified: 2020-05-31
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=shortbred
#SBATCH --array=1-712%50
#SBATCH --mem=4G
#SBATCH --output=slurm_out/shortbred_reads/x_shrtbrd_%a.out
#SBATCH --error=slurm_out/shortbred_reads/y_shrtbrd_%a.out

module load shortbred

# Basedir
basedir="$PWD"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d12_shortbred_reads"

# Make output directory
mkdir -p ${outdir}

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

# Start debug mode (will send commands to outfile)
set -x

time shortbred_quantify.py --markers /scratch/gdlab/bmahmud/shortBRED/shortBRED-identify/220529/220529_shortbred_FMG-CARD-NCBIAMR_clustId-95pc.faa\
 --wgs ${indir}/${sample}_fwd_clean.fastq ${indir}/${sample}_rev_clean.fastq\
 --results ${outdir}/${sample}_isolated.txt\
 --tmp tmp_shorbred_reads/${sample}_quantify

# Save error code for command
RC=$?

# End debug mode (will stop sending commands to outfile)
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occured for Sample ${sample}"
  exit $RC
fi