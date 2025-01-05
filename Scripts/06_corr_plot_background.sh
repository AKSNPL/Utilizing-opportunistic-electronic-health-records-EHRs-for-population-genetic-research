#!/bin/sh

## script to compute corrrelation between EHR and UKB in background

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=48:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=30

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-1%1

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/01_corr_background/obtain_UKB_EHR_corr-%x-%j.out

## change directory
cd /Aakash/12_GWAS_EHR_UKB/

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"

echo "Node ID: $SLURM_NODELIST"

## run the R script

## This is the container to be used
R_CONTAINER='/programs/all-inclusive-rstudio-apptainer/sif/all_inclusive_rstudio_4.3.1.sif'

# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='/Aakash/12_GWAS_EHR_UKB/scripts/05_correlation_plots.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="/Aakash/"

## The container 
singularity exec \
  --bind $BIND_DIR \
  $R_CONTAINER Rscript $R_SCRIPT
  
echo 'Done'
