#!/bin/bash
#===============================================================================
#
# File Name    : s17_5_inStrainCompare.sh
# Description  : This script will run instrain compare of sam files
# Usage        : sbatch s17_5_inStrainProfile.sh
# Author       : Kimberley Sukhum, kvsukhum@wustl.edu
# Version      : 1.0
# Created On   : 2021-08-09
# Last Modified: 2021-09-09
# Modified By  : Sanjam Sawhney, ssawhney@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=instrainCompare
#SBATCH --mem=200G
#SBATCH --output=slurm_out/inStrain_Compare_subset/x_inStrain_Compare_%a.out
#SBATCH --error=slurm_out/inStrain_Compare_subset/y_inStrain_Compare_%a.err
#SBATCH --cpus-per-task=8

#Load dependencies
eval $( spack load --sh py-instrain@1.5.7 )
#module load python
#module load samtools/1.9
#. /opt/apps/labs/gdlab/envs/instrain/1.5.4/bin/activate

# samplelist.txt contains a list of all file names
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d17_inStrain/Profile_subset"
outdir="${basedir}/d17_inStrain/Compare_subset"

#make output directory
mkdir -p ${outdir}

set -x

for i in `less ${basedir}/Mapping_v4_allfc_mod.txt`
do
  [ -f ${indir}/${i}.IS/output/${i}.IS_genome_info.tsv ] || mv ${indir}/${i}.IS ${indir}/${i}.IS.none
done

time inStrain compare -i ${indir}/*.IS/ \
    -o ${outdir}/${ID}.IS.COMPARE \
    -p ${SLURM_CPUS_PER_TASK} \
    --store_mismatch_locations \
    -s ${basedir}/d17_inStrain/MAGs_subset/MAGs_cat.stb \
    --database_mode

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi