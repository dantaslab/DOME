#!/bin/bash
#===============================================================================
# File Name    : s12_shortbred_ORFs.sh
# Description  : This script will profile sequence data and isolate genomes. 
# Usage        : sbatch s12_shortbred_ORFs.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : 2020-05-31
# Last Modified: 2020-05-31
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=shortbred
#SBATCH --array=1
#SBATCH --mem=4G
#SBATCH --output=slurm_out/shortbred_ORFs/x_shrtbrd_%a.out
#SBATCH --error=slurm_out/shortbred_ORFs/y_shrtbrd_%a.out

module load shortbred

# Basedir
basedir="$PWD"
indir="${basedir}/d08_prokka_samples"
outdir="${basedir}/d12_shortbred_ORFs"

# Make output directory
mkdir -p ${outdir}

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

# Start debug mode (will send commands to outfile)
set -x

time shortbred_quantify.py --markers /scratch/gdlab/bmahmud/shortBRED/shortBRED-identify/220529/220529_shortbred_FMG-CARD-NCBIAMR_clustId-95pc.faa\
 --genome ${indir}/${sample}/${sample}.faa\
 --results ${outdir}/${sample}_isolated.txt\
 --usearch /opt/apps/usearch/11.0.667/usearch\
 --maxhits 0\
 --maxrejects 0\
 --tmp tmp_shorbred_ORFs/${sample}_quantify

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