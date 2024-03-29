#!/bin/bash
#===============================================================================
#
# File Name    : s07_bin_dasToolPrep_maxbin.sh
# Description  : Prep the maxbin2 output for dasTool
# Usage        : sbatch s07_bin_dasToolPrep_maxbin.sh
# Author       : Robert Thaenert
# Version      : 1.0
# Created On   : 5/10/2022
# Last Modified: 5/10/2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=dTprep_maxbin
#SBATCH --array=2-283%50
#SBATCH --mem=1G
#SBATCH --output=slurm_out/bin_dasToolPrep_maxbin/x_dTprepMaxbin_%a.out
#SBATCH --error=slurm_out/bin_dasToolPrep_maxbin/y_dTprepMaxbin_%a.err

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"

#set output directory
outdir="${basedir}/d07_bin/maxbin2"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v5_subjects_mod.txt`

set -x

grep '>' ${outdir}/${sample}/${sample}*fasta > ${outdir}/${sample}/${sample}_step1.txt

sed 's/:>/\t/g' ${outdir}/${sample}/${sample}_step1.txt > ${outdir}/${sample}/${sample}_step2.txt
sed 's/.fasta//g' ${outdir}/${sample}/${sample}_step2.txt > ${outdir}/${sample}/${sample}_step3.txt
sed 's,.*/,,g' ${outdir}/${sample}/${sample}_step3.txt > ${outdir}/${sample}/${sample}_step4.txt

awk -F'\t' -v OFS='\t' '{print $2,$1}' ${outdir}/${sample}/${sample}_step4.txt > ${outdir}/${sample}/${sample}_step5.txt

rm ${outdir}/${sample}/${sample}_step1.txt
rm ${outdir}/${sample}/${sample}_step2.txt
rm ${outdir}/${sample}/${sample}_step3.txt
rm ${outdir}/${sample}/${sample}_step4.txt

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi