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
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/08_GWAS_step1_indications_all/Step1-%j.log.out
#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/08_GWAS_step1_indications_all/Step1-%j.log.err


## general 
source /source_proj.sh

source /source_ukb.sh

## get datasets 

### my own dir
export aakash=${ppl}/Aakash
export current_proj=${aakash}/12_GWAS_EHR_UKB
export input=${current_proj}/input/05_indications_phe_cov
export step1=${current_proj}/output/14_GWAS_step1_ind_all
export tmpdir=${current_proj}/tmpdir

names="all"

## change to relevant directory
cd ${current_proj}

${progs_dir}/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
--step 1 \
--bed ${pruned_arr}/ukb_allChrs.pruned \
--extract ${pruned_arr}/ukb_allChrs.prune.in \
--keep ${input}/EUR_panukbb_regenie_format.id \
--phenoFile ${input}/phenotypes_${names}.txt \
--covarFile ${input}/covariates.txt \
--phenoColList ALP,ALT,Albumin,Basophills,CRP,Calcium,Cholesterol,Creatinine,DBP,Eosinophills,FEV1,FVC,Glucose,HDL,Haematocritperc,Haemoglobinconc,HbA1c,Lymphocytes,MCHbconc,MCV,Monocytes,Neutrophills,Platelets,RBC,SBP,Totalbilirubin,Urea,WBC,Weight,Height,Triglycerides \
--bt \
--threads 30 \
--bsize 1000 \
--lowmem \
--lowmem-prefix ${tmpdir}/regenie_tmp_preds_${names} \
--use-relative-path --gz \
--out ${step1}/pred_step1_${names}
