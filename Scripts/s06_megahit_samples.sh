#!/bin/bash                                                                                                                                                                                                
#===============================================================================                                                                                                                           
# Name         : s06_megahit_samples.sh                                                                                                                                                                            
# Description  : This script will run MEGAHIT on samples in parallel                                                                                                                                       
# Usage        : sbatch s06_megahit_samples.sh                                                                                                                                                                     
# Author       : Kevin Blake, kevin.blake@wustl.edu                                                                                                                                                        
# Version      : 1.0                                                                                                                                                                                       
# Created On   : 26 Feb 2022                                                                                                                                                                               
# Modified On  : Sat Feb 26 20:55:00 CST 2022                                                                                                                                                              
#===============================================================================                                                                                                                           
#Submission script for HTCF                                                                                                                                                                                
#SBATCH --job-name=megahit_1000                                                                                                                                                                            
#SBATCH --array=1-712%50                                                                                                                                                                                    
#SBATCH --mem=32G                                                                                                                                                                                         
#SBATCH --cpus-per-task=8                                                                                                                                                                                  
#SBATCH --output=slurm_out/megahit_samples/x_megahit_%a.out                                                                                                                                                     
#SBATCH --error=slurm_out/megahit_samples/y_megahit_%a.err                                                                                                                                                                                                                                                                                                                      

module load megahit/1.1.4

# Basedir                                                                                                                                                                                                  
basedir="/scratch/gdlab/bmahmud/DOME/shotgun"
indir="${basedir}/d03_clean_reads"
outdir="${basedir}/d06_megahit_samples"

# Make output directory                                                                                                                                                                                    
mkdir -p ${outdir}

# Read in the slurm array task                                                                                                                                                                             
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v4_allfc_mod.txt`

# modifiied to utilize /tmp/ on the node                                                                                                                                                                   
set -x
time megahit \
     -1 ${indir}/${sample}_fwd_clean.fastq \
     -2 ${indir}/${sample}_rev_clean.fastq \
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