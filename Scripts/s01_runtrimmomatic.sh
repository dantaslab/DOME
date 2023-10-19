#!/bin/bash
#==============================================================================
#
# File Name    : s01_runtrimmomatic.sh
# Description  : This script will run trimmomatic in parallel
# Usage        : sbatch s01_runtrimmomatic.sh
# Author       : Bejan Mahmud, bmahmud@wustl.edu
# Version      : 2.0
# Created On   : 2020-04-06
# Last Modified: 2020-04-14
#==============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=trimmomatic
#SBATCH --array=1-712%50
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm_out/trim/x_trimReads_%a.out
#SBATCH --error=slurm_out/trim/y_trimReads_%a.err

module load trimmomatic/0.38

#store adapter specific paths
adapt="/opt/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa"

#basedir if you want
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d00_rawreads"
outdir="${basedir}/d01_trimReads"

#make output directory
mkdir -p ${outdir}

export JAVA_ARGS="-Xmx15000M"

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v1_allfc.txt`

set -x
time java -jar $TRIMMOMATIC_HOME/trimmomatic-0.38.jar PE -phred33 -trimlog ${outdir}/Paired_${sample}_trimlog.txt ${indir}/${sample}_R1.fastq ${indir}/${sample}_R2.fastq ${outdir}/${sample}_FW_trimmed_paired.fastq ${outdir}/${sample}_FW_trimmed_unpaired.fastq ${outdir}/${sample}_RV_trimmed_paired.fastq ${outdir}/${sample}_RV_trimmed_unpaired.fastq ILLUMINACLIP:${adapt}:2:30:10:1:true SLIDINGWINDOW:4:20 LEADING:10 TRAILING:10 MINLEN:60
RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured!"
    exit $RC
fi