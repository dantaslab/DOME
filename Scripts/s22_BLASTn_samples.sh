#!/bin/bash
#===============================================================================
# File Name    : s22_BLASTn_samples.sh
# Description  : This script will run BLASTn on all sample assemblies against all sample databases (all-vs-all)
# Usage        : sbatch s22_BLASTn_samples.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTn
#SBATCH --array=2-712%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_out/blastn_samples/x_blastdb_%a.out
#SBATCH --error=slurm_out/blastn_samples/y_blastdb_%a.err

#module load blast-plus/2.11.0-python-3.6.5
eval $( spack load --sh blast-plus@2.12.0 )

basedir="$PWD"
indir="${basedir}/d06_megahit_samples"
outdir="${basedir}/d22_BLASTn_samples/BLASTn"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

#make output directory
mkdir -p ${outdir}

set -x

for i in {1..712}
do
  [ ${i} -eq ${SLURM_ARRAY_TASK_ID} ] && continue
  subject=`sed -n ${i}p ${basedir}/Mapping_v4_allfc_mod.txt`
  echo "Blasting ${sample} with ${subject}"
  time blastn \
          -query ${indir}/${sample}/final.contigs.modhead.2line.fa \
          -db ${basedir}/d22_BLASTn_samples/db/${subject}/${subject} \
          -outfmt 6 \
          -perc_identity 99 \
          -out ${outdir}/${sample}_${subject}.txt
done

cat ${outdir}/${sample}_*.txt > ${outdir}/merged_${sample}.txt
rm ${outdir}/${sample}_*.txt

RC=$?

set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi