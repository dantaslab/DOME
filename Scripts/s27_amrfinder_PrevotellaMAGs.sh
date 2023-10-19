#!/bin/bash

#===============================================================================
# File Name    : s27_amrfinder_PrevotellaMAGs.sh
# Description  : Annotations ARGs and VFs in assembled genomes or contigs
# Usage        : sbatch s27_amrfinder_PrevotellaMAGs.sh
# Author       : Luke Diorio-Toth
# Version      : 1.3
# Created On   : Tue Oct 11 14:49:35 CDT 2022
# Last Modified: Tue Oct 11 14:49:35 CDT 2022
#===============================================================================

#SBATCH --job-name=amrfinder_Prevotella
#SBATCH --array=1-101%20
#SBATCH --cpus-per-task=2
#SBATCH --mem=500M
#SBATCH --output=slurm_out/amrfinder_Prevotella/x_amrfinder_%a.out
#SBATCH --error=slurm_out/amrfinder_Prevotella/y_amrfinder_%a.err

eval $( spack load --sh amrfinder@3.10.42 )

basedir="$PWD"
indir="${basedir}/d25_bakta/Prevotella"
outdir="${basedir}/d27_amrfinder/Prevotella"

# because there is inconsistent formatting to GFF files, specify the format
# options include prokka, bakta, pgap, etc. (see docs for the full list)
annotation_format="bakta"

#make output directory and read in the slurm array task
mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v7_PrevotellaMAGs.txt`

set -x
time amrfinder --plus \
  -n ${indir}/${sample}/${sample}.fna \
  -p ${indir}/${sample}/${sample}.faa \
  -g ${indir}/${sample}/${sample}.gff3 \
  -a ${annotation_format} \
  --name ${sample} \
  -o ${outdir}/${sample}_out.tsv \
  --threads ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi