#!/bin/bash
#===============================================================================
# File Name    : s25_bakta_RC9MAGs.sh
# Description  : Annotation of assembled bacterial genomes, MAGs, and plasmids
# Usage        : sbatch s25_bakta_RC9MAGs.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.0
# Modified     : Fri Oct 14 08:14:59 CDT 2022
# Created      : Fri Oct 14 08:14:59 CDT 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=bakta_RC9
#SBATCH --array=1-18%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm_out/bakta_RC9/x_bakta_%a.out
#SBATCH --error=slurm_out/bakta_RC9/y_bakta_%a.err

eval $( spack load --sh py-bakta )

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/allHQandMQ"
outdir="${basedir}/d25_bakta/RC9"
dbdir="/ref/gdlab/data/bakta_db/db"
mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v7_RC9MAGs.txt`

set -x
time bakta \
  --db ${dbdir} \
  --min-contig-length 200 \
  --prefix ${sample} \
  --output ${outdir}/${sample}.fasta \
  --threads ${SLURM_CPUS_PER_TASK} \
  --output ${outdir}/${sample} \
  ${indir}/${sample}.fa
RC=$?
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi