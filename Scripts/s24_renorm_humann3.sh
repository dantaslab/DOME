#!/bin/bash
#===============================================================================
# Description  : This script will create metabolic pathway profiles for the input sequences. 
# Usage        : sbatch s24_renorm_humann3.sh
# Author       : Alaric D'Souza, alaric.dsouza@wustl.edu
# Version      : 2.5
# Created On   : 2018-10-19
# Last Modified: Mon Dec 17 23:03:40 CST 2018
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=humann
#SBATCH --array=114,145,160,177,190,225,22,25,278,281,4,73,94
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --output=slurm_out/humann3_renorm/x_humann3_%a.out
#SBATCH --error=slurm_out/humann3_renorm/y_humann3_%a.err

module load python
. /opt/apps/labs/gdlab/envs/humann/3.0.0a4/bin/activate

# Basedir
basedir="$PWD"
indir="${basedir}/d24_humann3"
outdir="${basedir}/d24_humann3/norm"
#nucleotide_db='/scratch/ref/gdlab/humann3/'
#protein_db="${basedir}/uniref90"

# Make output directory
mkdir -p ${outdir}

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/Mapping_v8_humanfc_mod.txt`

# Start debug mode (will send commands to outfile)
set -x

# Run and time command
humann_renorm_table --input ${indir}/${sample}_pathabundance.tsv --output ${outdir}/${sample}_pathabundance-cpm.tsv --units cpm --update-snames
#humann2_rename_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_uniref_names.tsv -n uniref90
#humann2_regroup_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_ko.tsv -g uniref90_ko
#humann2_rename_table -i ${indir}/${sample}_gf_ko.tsv -o ${indir}/${sample}_gf_ko_names.tsv -n kegg-orthology
#humann2_renorm_table --input ${indir}/${sample}_gf_ko_names.tsv --output ${indir}/${sample}_gf_ko_names_ra.tsv --units relab --update-snames
#humann2_regroup_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_eggnog.tsv -g uniref90_eggnog

#humann2_rename_table -i ${indir}/${sample}_gf_eggnog.tsv -o ${indir}/${sample}_gf_eggnog_names.tsv -n eggnog
#humann2_renorm_table --input ${indir}/${sample}_gf_eggnog_names.tsv --output ${indir}/${sample}_gf_eggnog_names_ra.tsv --units relab --update-snames
#humann2_regroup_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_pfam.tsv -g uniref90_pfam
#humann2_rename_table -i ${indir}/${sample}_gf_pfam.tsv -o ${indir}/${sample}_gf_pfam_names.tsv -n pfam
#humann2_renorm_table --input ${indir}/${sample}_gf_pfam_names.tsv --output ${indir}/${sample}_gf_pfam_names_ra.tsv --units relab --update-snames
#humann2_regroup_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_level4ec.tsv -g uniref90_level4ec
#humann2_rename_table -i ${indir}/${sample}_gf_level4ec.tsv -o ${indir}/${sample}_gf_level4ec_names.tsv -n ec
#humann2_renorm_table --input ${indir}/${sample}_gf_level4ec_names.tsv --output ${indir}/${sample}_gf_level4ec_names_ra.tsv --units relab --update-snames
#humann2_regroup_table -i ${indir}/${sample}_genefamilies.tsv -o ${indir}/${sample}_gf_infogo1000.tsv -g uniref90_infogo1000
#humann2_rename_table -i ${indir}/${sample}_gf_infogo1000.tsv -o ${indir}/${sample}_gf_infogo1000_names.tsv -n infogo1000
#humann2_renorm_table --input ${indir}/${sample}_gf_infogo1000_names.tsv --output ${indir}/${sample}_gf_infogo1000_names_ra.tsv --units relab --update-snames
#humann2_renorm_table --input ${indir}/${sample}_gf_uniref_names.tsv --output ${indir}/${sample}_gf_uniref_names_ra.tsv --units relab --update-snames
# Variations
# Instead of "--bypass-translated-search" which skips the protein comparison, you can choose to use the flag "--bypass-nucleotide-search", which bypasses gene comparison
# --gap-fill on will assume that is just one part of the metabolic pathway is missing, then it will give the benefit of doubt

# Save error code for command
RC=$?

# End debug mode (will stop sending commands to outfile)
set +x

# Output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi