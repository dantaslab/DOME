#!/bin/bash                                                                                                                                                                                                
#===============================================================================                                                                                                                           
# Name         : s06_megahit.sh                                                                                                                                                                            
# Description  : This script will run MEGAHIT on samples in parallel                                                                                                                                       
# Usage        : sbatch s06_megahit.sh                                                                                                                                                                     
# Author       : Kevin Blake, kevin.blake@wustl.edu                                                                                                                                                        
# Version      : 1.0                                                                                                                                                                                       
# Created On   : 26 Feb 2022                                                                                                                                                                               
# Modified On  : Sat Feb 26 20:55:00 CST 2022                                                                                                                                                              
#===============================================================================                                                                                                                           
#Submission script for HTCF                                                                                                                                                                                
#SBATCH --job-name=megahit_1000                                                                                                                                                                            
#SBATCH --array=1-284%50                                                                                                                                                                                      
#SBATCH --mem=128G                                                                                                                                                                                         
#SBATCH --cpus-per-task=8                                                                                                                                                                                  
#SBATCH --output=slurm_out/megahit/x_megahit_%a.out                                                                                                                                                     
#SBATCH --error=slurm_out/megahit/y_megahit_%a.err                                                                                                                                                                                                                                                                                                                      

module load megahit/1.1.4

# Basedir                                                                                                                                                                                                  
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d04_subjectreads"
outdir="${basedir}/d06_metahit"

# Make output directory                                                                                                                                                                                    
mkdir -p ${outdir}

# Read in the slurm array task                                                                                                                                                                             
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v3_subjects.txt`

# modifiied to utilize /tmp/ on the node                                                                                                                                                                   
set -x
time megahit \
     -1 ${indir}/${sample}_R1.fastq \
     -2 ${indir}/${sample}_R2.fastq \
     --min-contig-len 1000 \
     -t ${SLURM_CPUS_PER_TASK} \
     -o ${outdir}/${sample}

RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured in ${sample}!"
    exit $RC
fi