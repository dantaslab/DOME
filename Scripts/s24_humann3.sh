#!/bin/bash
#===============================================================================
# File Name    : s24_humann3.sh
# Description  : This script will run humann3 on concatenated FW/RV reads
# Usage        : sbatch s24_humann3.sh
# Author       : Luke Diorio-Toth
# Version      : 1.0
# Created On   : Wed Jan  6 11:49:47 CST 2021
# Last Modified: Wed Jan  6 11:49:47 CST 2021
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=humann3
#SBATCH --array=177
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --output=slurm_out/humann3/x_humann3_%a.out
#SBATCH --error=slurm_out/humann3/y_humann3_%a.err

# load python and activate humann venv
module load python
. /opt/apps/labs/gdlab/envs/humann/3.0.0a4/bin/activate

basedir="$PWD"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d24_humann3"
tmpdir="${basedir}/tmpHUMANN"

#make output directory
mkdir -p ${outdir}
mkdir -p ${tmpdir}

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v8_humanfc_mod.txt`

set -x # Start debug mode (will send commands to outfile)

cat ${indir}/${sample}_fwd_clean.fastq ${indir}/${sample}_rev_clean.fastq > ${tmpdir}/${sample}_catted.fastq

time humann \
        --input ${tmpdir}/${sample}_catted.fastq \
        --output ${outdir} \
        --output-basename ${sample} \
        --threads ${SLURM_CPUS_PER_TASK}

rm ${tmpdir}/${sample}_catted.fastq
rm -R ${outdir}/${sample}_humann_temp

RC=$? # Save error code for command
set +x # End debug mode (will stop sending commands to outfile)

deactivate # deactivate venv

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi