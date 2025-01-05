#!/bin/sh

## script to runs econd step of REGENIE

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=47:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=12

#SBATCH --mem-per-cpu=12G

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-650%20

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/05_proxies/all_obtain_reg_sent-%x-%j.out

## change directory
cd /Aakash/12_GWAS_EHR_UKB/

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' output/08_res_sentinels_all/region.query_sukb_UKB.txt)"
low="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' output/08_res_sentinels_all/region.query_sukb_UKB.txt)"
upp="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' output/08_res_sentinels_all/region.query_sukb_UKB.txt)"
var="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' output/08_res_sentinels_all/region.query_sukb_UKB.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Chromosome : ${chr} | Start : ${low} | End : ${upp} | MarkerName(s) : ${var}"

echo "${chr}:${low}-${upp}"


## run the R script
R_CONTAINER='/programs/all-inclusive-rstudio-apptainer/sif/all_inclusive_rstudio_4.3.1.sif'
R_SCRIPT='/Aakash/12_GWAS_EHR_UKB/scripts/14_obtain_regional_proxies.R'
BIND_DIR="/./"
BIND_DIR2="/imputed/bgen_files_44448/"
singularity exec \
  --bind $BIND_DIR,$BIND_DIR2 \
  $R_CONTAINER $R_SCRIPT ${chr} ${low} ${upp} ${var}

echo "Done!"
exit 0
