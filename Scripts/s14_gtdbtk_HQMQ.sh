#!/bin/bash
#===============================================================================
# File Name    : s14_gtdbtk_HQMQ.sh
# Description  : This script will run the gtdbtk classify workflow
# Usage        : sbatch s14_gtdbtk_HQMQ.sh
# Author       : Luke Diorio-Toth
# Version      : 1.0
# Created On   : Mon Apr  4 16:23:16 CDT 2022
# Last Modified: Mon Apr  4 16:23:16 CDT 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=gtdbtk_classify
#SBATCH --array=1-49%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=250G
#SBATCH --output=slurm_out/gtdbtk_HQMQ/x_gtdbtk_%a.out
#SBATCH --error=slurm_out/gtdbtk_HQMQ/y_gtdbtk_%a.err


module load miniconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /opt/apps/labs/gdlab/envs/gtdbtk/1.7.0

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/subset${SLURM_ARRAY_TASK_ID}"
outdir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/gtdbtk/subset${SLURM_ARRAY_TASK_ID}"

mkdir -p ${outdir}

set -x
time gtdbtk classify_wf \
        --genome_dir ${indir} \
        -x fa \
        --out_dir ${outdir} \
        --cpus ${SLURM_CPUS_PER_TASK}

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi