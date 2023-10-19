#!/bin/bash
#===============================================================================
# File Name    : s22_BLASTn_RC9MAGs.sh
# Description  : This script will run BLASTn on all sample assemblies against all sample databases (all-vs-all)
# Usage        : sbatch s22_BLASTn_RC9MAGs.sh
# Author       : Bejan Mahmud
# Version      : 1.0
# Created On   : Tue Jun 21 2022
# Last Modified: Tue Jun 21 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=BLASTn_RC9
#SBATCH --array=1-18
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_out/blastn_MAGs/x_blastn_RC9_%a.out
#SBATCH --error=slurm_out/blastn_MAGs/y_blastn_RC9_%a.err

#module load blast-plus/2.11.0-python-3.6.5
#eval $( spack load --sh blast-plus@2.12.0 )
eval $( spack load --sh /u7ssbm4 )

basedir="$PWD"
indir="${basedir}/d07_bin/dastool/HQandMQ_MAGs/allHQandMQ"
outdir="${basedir}/d22_BLASTn_MAGs/RC9/BLASTn"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v7_RC9MAGs.txt`

#make output directory
mkdir -p ${outdir}

set -x

for i in {1..18}
do
  [ ${i} -eq ${SLURM_ARRAY_TASK_ID} ] && continue
  subject=`sed -n ${i}p ${basedir}/Mapping_v7_RC9MAGs.txt`
  echo "Blasting ${sample} with ${subject}"
  time blastn \
          -query ${indir}/${sample}.fa \
          -db ${basedir}/d22_BLASTn_MAGs/RC9/db/${subject}/${subject} \
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