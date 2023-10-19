#!/bin/bash

#===============================================================================
# File Name    : s17_1_dRep_parse_stb.sh
# Description  : This script will use quast to do quality control on assemblies
# Usage        : sbatch s17_1_dRep_parse_stb.sh
# Author       : Sanjam Sawhney
# Version      : 1.0
# Created On   : Wed Jul 24 15:46:45 CDT 2019
# Last Modified: Tue Aug 13 19:24:00 CDT 2019
#===============================================================================
# Submission script for HTCF
#SBATCH --cpus-per-task=8
#SBATCH --job-name=parse_stb
#SBATCH --mem=5G
#SBATCH --output=slurm_out/inStrain_stb/x_inStrain_stb_%a.out
#SBATCH --error=slurm_out/inStrain_stb/y_inStrain_stb_%a.err

#load dependencies
module load prodigal
module load mash
module load MUMmer
module load fastANI
module load ANIcalculator
module load centrifuge
module load openmpi/3.1.3-python-2.7.15
module load hmmer/3.1b1-python-2.7.15
module load pplacer/1.1.alpha19

. /opt/apps/labs/gdlab/envs/drep/3.2.2/bin/activate

basedir="$PWD"
indir="${basedir}/d09_dRep/dRep_98_nc01"
outdir="${basedir}/d17_inStrain/MAGs"

#make output directory
mkdir -p ${outdir}

#enter debug mode
set -x

#Create concatenated multifasta file for all dREP MAG winners for downstream Bowtie2
time cat ${indir}/dereplicated_genomes/MAGs_final_renamed/*.fa > ${outdir}/MAGs_cat.fasta

#Run stb script
time parse_stb.py --reverse -f ${indir}/dereplicated_genomes/MAGs_final_renamed/*.fa  -o ${outdir}/MAGs_cat.stb

#Create concatenated multi-fna file for all dREP MAG winner Prodigal annotations for downstream inStrain profile
time cat ${indir}/data/prodigal/MAGs_final_renamed/*.fna > ${outdir}/MAGs_cat.fna

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