#!/bin/bash

#===============================================================================
# File Name    : s26_panaroo_TreponemaMAGs.sh
# Description  : Runs the panaroo pangenome tool on genomes annotated by bakta
# Usage        : sbatch s26_panaroo_TreponemaMAGs.sh
# Author       : Luke Diorio-Toth
# Version      : 1.4
# Created On   : Thu Aug 15 14:44:06 CDT 2019
# Last Modified: Tue Oct 18 15:32:58 CDT 2022
#===============================================================================

#SBATCH --job-name=panaroo_Treponema
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm_out/panaroo/x_panaroo_Treponema_%A.out
#SBATCH --error=slurm_out/panaroo/y_panaroo_Treponema_%A.err

eval $( spack load --sh py-panaroo )

basedir="$PWD"
indir="${basedir}/d25_bakta/Treponema"
outdir="${basedir}/d26_panaroo/Treponema"

mkdir -p ${outdir}

# make text file with paths to gff files
find ${indir} -iname "*.gff3" >> gffs.${SLURM_JOB_ID}.tmp

set -x
time panaroo \
        -i ${basedir}/gffs.${SLURM_JOB_ID}.tmp \
        -o ${outdir} \
        --clean-mode strict \
        --core_threshold 0.98 \
        -a core \
        --aligner mafft \
        --remove-invalid-gene \
        -t ${SLURM_CPUS_PER_TASK}

RC=$?

rm ${basedir}/gffs.${SLURM_JOB_ID}.tmp

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi