#!/bin/bash
#===============================================================================
# Name         : s13_checkm.sh
# Description  : This script will run checkm on a unicycler output dir
# Usage        : sbatch s13_checkm.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.5
# Created On   : 2018_03_25
# Last Modified: Sat Aug 22 12:52:49 CDT 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkm
#SBATCH --array=1-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --output=slurm_out/checkm/x_checkm_%a.out
#SBATCH --error=slurm_out/checkm/y_checkm_%a.err

module load openmpi/3.1.3-python-2.7.15
module load py-checkm-genome/1.0.13-python-2.7.15

# workaround if there is error loading shared libraries
export LD_LIBRARY_PATH=/opt/slurm/lib:$LD_LIBRARY_PATH

basedir="$PWD"
indir="${basedir}/d09_dRep/dRep_99/dereplicated_genomes/set${SLURM_ARRAY_TASK_ID}"
outdir="${basedir}/d13_dRepMAGs_checkm/set${SLURM_ARRAY_TASK_ID}"

mkdir -p ${outdir}

set -x
time  checkm lineage_wf -f ${outdir}/output_set${SLURM_ARRAY_TASK_ID}.txt \
        -t ${SLURM_CPUS_PER_TASK} \
        -x fa \
        --tab_table \
        ${indir} ${outdir}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi