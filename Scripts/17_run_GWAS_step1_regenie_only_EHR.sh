#!/bin/sh

## script to runs 1st step of REGENIE

#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 120
#SBATCH --array=1-1%1
#SBATCH --mail-type=FAIL
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/06_GWAS_step1_indications_EHR_only/Step1-%j.log.out
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/06_GWAS_step1_indications_EHR_only/Step1-%j.log.err


## general 
source /source_proj.sh

source /source_ukb.sh

## get datasets 

### my own dir
export aakash=${ppl}/Aakash
export current_proj=${aakash}/12_GWAS_EHR_UKB
export input=${current_proj}/input/05_indications_phe_cov
export step1=${current_proj}/output/12_GWAS_step1_ind_EHR_only
export tmpdir=${current_proj}/tmpdir

names="ehr"

## change to relevant directory
cd ${current_proj}

${progs_dir}/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
--step 1 \
--bed ${pruned_arr}/ukb_allChrs.pruned \
--extract ${pruned_arr}/ukb_allChrs.prune.in \
--keep ${input}/EUR_panukbb_regenie_format.id \
--phenoFile ${input}/phenotypes_${names}.txt \
--covarFile ${input}/covariates.txt \
--phenoColList ${names} \
--bt \
--threads 30 \
--bsize 1000 \
--lowmem \
--lowmem-prefix ${tmpdir}/regenie_tmp_preds_${names} \
--use-relative-path --gz \
--out ${step1}/pred_step1_${names}
