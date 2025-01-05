#!/bin/sh

## script to runs econd step of REGENIE

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=12:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=30

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-650%100

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/05_proxies/obtain_reg_sent_only_EHR_qt-%x-%j.out

## change directory
cd /Aakash/12_GWAS_EHR_UKB/

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' output/24_res_sentinels_only_EHR_qt/region.query_all_only_EHR_qt.txt)"
low="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' output/24_res_sentinels_only_EHR_qt/region.query_all_only_EHR_qt.txt)"
upp="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' output/24_res_sentinels_only_EHR_qt/region.query_all_only_EHR_qt.txt)"
var="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' output/24_res_sentinels_only_EHR_qt/region.query_all_only_EHR_qt.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Chromosome : ${chr} | Start : ${low} | End : ${upp} | MarkerName(s) : ${var}"

echo "${chr}:${low}-${upp}"

## set up R environment
eval "$(/opt/conda/bin/conda shell.bash hook)"
conda activate r_env

## run the R script
scripts/42_obtain_regional_proxies.R ${chr} ${low} ${upp} ${var} 

## deactivate conda env
conda deactivate
