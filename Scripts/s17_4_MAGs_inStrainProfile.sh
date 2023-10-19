#!/bin/bash
#===============================================================================
#
# File Name    : s17_4_MAGs_inStrainProfile.sh
# Description  : This script will run instrain profile of sam files
# Usage        : sbatch s17_4_MAGs_inStrainProfile.sh
# Author       : Kimberley Sukhum, kvsukhum@wustl.edu
# Version      : 1.0
# Created On   : 2021-08-09
# Last Modified: 2021-09-09
# Modified By  : Sanjam Sawhney, ssawhney@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=instrain
#SBATCH --array=7261-7458%100
#SBATCH --mem=20G
#SBATCH --output=slurm_out/inStrain_Profile_MAGs/x_inStrain_Profile_%a.out
#SBATCH --error=slurm_out/inStrain_Profile_MAGs/y_inStrain_Profile_%a.err
#SBATCH --cpus-per-task=8

#Load dependencies
eval $( spack load --sh py-instrain@1.5.7 )
#module load python
#module load samtools/1.9
#module load prodigal
#. /opt/apps/labs/gdlab/envs/instrain/1.5.4/bin/activate

#basedir
basedir="$PWD"

mag=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d17_inStrain_MAGs/221008_pairs_MAGs.txt`
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/d17_inStrain_MAGs/221008_pairs_samples.txt`

indir="${basedir}/d17_inStrain_MAGs/${mag}/BowtieAlign"
outdir="${basedir}/d17_inStrain_MAGs/${mag}/Profile"
WD="/tmp/bck-${SLURM_JOBID}-${SLURM_ARRAY_TASK_ID}"

#make output directory
mkdir -p ${outdir}
mkdir ${WD}

set -x

cp ${indir}/${sample}_mappedreads_sorted.bam ${basedir}/d17_inStrain_MAGs/${mag}/MAGs/MAGs_cat.fasta ${basedir}/d17_inStrain_MAGs/${mag}/MAGs/MAGs_cat.fna ${basedir}/d17_inStrain_MAGs/${mag}/MAGs/MAGs_cat.stb ${WD}

pushd ${WD}

time inStrain profile ${sample}_mappedreads_sorted.bam MAGs_cat.fasta \
    -o output/${sample}.IS \
    -g MAGs_cat.fna \
    -s MAGs_cat.stb \
    -p ${SLURM_CPUS_PER_TASK} \
    --skip_plot_generation \
    --database_mode

popd

cp -a ${WD}/output/${sample}.IS ${outdir}
rm -rf ${WD}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi