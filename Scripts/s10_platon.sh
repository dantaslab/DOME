#!/bin/bash
#===============================================================================
# File Name    : s10_platon.sh
# Description  : This script will run find plasmids in short-read assemblies
# Usage        : sbatch s10_platon.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Created On   : Tue Jul  7 21:45:13 CDT 2020
# Modified On  : Tue Jul  7 21:45:13 CDT 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=platon
#SBATCH --array=648
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --output=slurm_out/platon/x_platon_%a.out
#SBATCH --error=slurm_out/platon/y_platon_%a.err

#load dependencies
module load python
module load prodigal
module load diamond
module load ncbi-blast
module load mummer
module load hmmer
module load infernal
#load python venv
. /opt/apps/labs/gdlab/envs/platon/bin/activate

# workaround if there is error loading shared libraries
export LD_LIBRARY_PATH=/opt/slurm/lib:$LD_LIBRARY_PATH

# Basedir
basedir="$PWD"
indir="${basedir}/d06_megahit_samples"
outdir="${basedir}/d10_platon"

# Make output directory
mkdir -p ${outdir}

# read in sample name
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

set -x
time platon \
        --db /scratch/ref/gdlab/platon_db/latest \
        --output ${outdir}/${sample} \
        --verbose \
        --mode sensitivity \
        --threads ${SLURM_CPUS_PER_TASK} \
        --prefix ${sample} \
        ${indir}/${sample}/final.contigs.modhead.fa
RC=$?
set +x

deactivate

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi