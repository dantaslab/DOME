#!/bin/bash
#===============================================================================
# File Name    : shortbred_quantify.sh
# Description  : This script will qunatify ARGs using pre-processed shotgun sequencing data. 
# Usage        : sbatch shortbred_quantify.sh
# Author       : Bejan Mahmud
# Version      : 1.1
# Created On   : 2022-06-05
# Last Modified: 2023-04-24
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=shortbred
#SBATCH --array=1-712%50
#SBATCH --mem=4G
#SBATCH --output=slurm_out/shortbred_quantify/x_shrtbrd_%a.out
#SBATCH --error=slurm_out/shortbred_quantify/y_shrtbrd_%a.out

module load shortbred

# Basedir
basedir="$PWD"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d12_shortbred_quantify"

# Make output directory
mkdir -p ${outdir}

# Read in the sample name from the Mapping file
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping.txt`

# Start debug mode (will send commands to outfile)
set -x

time shortbred_quantify.py --markers 220529_shortbred_FMG-CARD-NCBIAMR_clustId-95pc.faa\
 --wgs ${indir}/${sample}_fwd.fastq ${indir}/${sample}_rev.fastq\
 --results ${outdir}/${sample}_ShortBRED.txt\
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