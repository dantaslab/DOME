#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_CONCOCT_CutContigs.sh
# Description  : Run cut contigs script of CONCOCT
# Usage        : sbatch s07_bin_CONCOCT_CutContigs.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 05/03/2022
# Last Modified: 05/02/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=concoct_cut
#SBATCH --cpus-per-task=1
#SBATCH --array=2-283%50
#SBATCH --mem=1G
#SBATCH --output=slurm_out/bin_concoct_cut/x_cut_%a.out
#SBATCH --error=slurm_out/bin_concoct_cut/y_cut_%a.err

module load concoct
module load python/3.6.5

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
indir="${basedir}/d06_megahit"
outdir="${basedir}/d07_bin/concoct"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

#make output directory
mkdir -p ${outdir}
mkdir -p ${outdir}/${sample}

set -x

python3 CONCOCT_scripts/cut_up_fasta.py ${indir}/${sample}/final.contigs.modhead.fa -c 10000 -o 0 --merge_last -b ${outdir}/${sample}/${sample}_contigs_10K.bed > ${outdir}/${sample}/${sample}_contigs_10K.fa

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi