#!/bin/sh

## script to run the 2nd  Step for REGENIE
#SBATCH --job-name 99_GWAS
#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=120G
#SBATCH --array=1-23%23

#SBATCH --output=/Aakash/12_GWAS_EHR_UKB/logs/11_GWAS_step2_UKB_only/test_Step2-%j.log
#SBATCH --error=/Aakash/12_GWAS_EHR_UKB/logs/11_GWAS_step2_UKB_only/test_Step2-%j.err


## general 
source /source_proj.sh

source /source_ukb.sh

## get datasets:IMPUTED 
export var=/imputed/variant_qc/EUR
### my own dir
export aakash=${ppl}/Aakash
export current_proj=${aakash}/12_GWAS_EHR_UKB
export input=${current_proj}/input/06_UKB_only_phe_cov
export step1=${current_proj}/output/16_GWAS_step1_UKB_only
export tmpdir=${current_proj}/tmpdir
export output=${current_proj}/output/17_GWAS_step2_UKB_only

names="all"

## a file from 1 ...23
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' /Aakash/12_GWAS_EHR_UKB/input/04_others/chr_all.txt)"

cd ${current_proj}

## create variables for loading easily!
 gen_ids="${input}/EUR_panukbb_regenie_format.id"
 cov_file="${input}/covariates.txt"
 pheno_file="${input}/phenotypes_${names}.txt" 
 sample_X="${ukb_imp}/bgen_files/ukb_cX_b0_v3.sample"
 sample_all="${ukb_imp}/bgen_files/ukb_c1_b0_v3.sample"
 step1_file="${step1}/pred_step1_${names}_pred.list"
 step2_file="${output}/GWAS_chr${chr}_${names}"



## get the chromosome
if [[ $chr -eq 23 ]]; then


  export chr="X"

   echo "echo chr $chr .."

  cat ${var}/output/ukb_imp_chr${chr}_qced.txt | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > ${tmpdir}/tmp.ex.chr${chr}_${names}.list

## run REGENIE (needs a remove flag for samples to exclude)
${progs_dir}/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
	--step 2 \
 	--bgen ${ukb_imp}/bgen_files/ukb_c${chr}_b0_v3.bgen \
 	--ref-first \
 	--extract ${tmpdir}/tmp.ex.chr${chr}_${names}.list \
 	--sample ${sample_X} \
 	--keep ${gen_ids} \
 	--phenoFile ${pheno_file} \
 	--covarFile ${cov_file} \
	--phenoColList ALP,ALT,Albumin,Basophills,CRP,Calcium,Cholesterol,Creatinine,DBP,Eosinophills,FEV1,FVC,Glucose,HDL,Haematocritperc,Haemoglobinconc,HbA1c,Lymphocytes,MCHbconc,MCV,Monocytes,Neutrophills,Platelets,RBC,SBP,Totalbilirubin,Urea,WBC,Weight,Height,Triglycerides \
	--qt \
 	--threads 30 \
 	--pred ${step1_file} \
 	--bsize 1000 \
	--use-relative-path --gz \
 	--out ${step2_file}

## remove variant inclusion list by chr number..
rm ${tmpdir}/tmp.ex.chr${chr}_${names}.list

else
echo "echo chromsome in run is.. $chr .."
  
  ## create variant inclusion list for the respective chromosome
cat ${var}/ukb_imp_EUR_chr${chr}_snpstat.out | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' - > ${tmpdir}/tmp.ex.chr${chr}_${names}.list


## run REGENIE (needs a remove flag for samples to exclude)
${progs_dir}/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl \
	--step 2 \
 	--bgen ${ukb_imp}/bgen_files/ukb_c${chr}_b0_v3.bgen \
 	--ref-first \
 	--extract ${tmpdir}/tmp.ex.chr${chr}_${names}.list \
 	--sample ${sample_all} \
 	--keep ${gen_ids} \
 	--phenoFile ${pheno_file} \
 	--covarFile ${cov_file} \
	--phenoColList ALP,ALT,Albumin,Basophills,CRP,Calcium,Cholesterol,Creatinine,DBP,Eosinophills,FEV1,FVC,Glucose,HDL,Haematocritperc,Haemoglobinconc,HbA1c,Lymphocytes,MCHbconc,MCV,Monocytes,Neutrophills,Platelets,RBC,SBP,Totalbilirubin,Urea,WBC,Weight,Height,Triglycerides \
	--qt \
 	--threads 30 \
 	--pred ${step1_file} \
 	--bsize 1000 \
	--use-relative-path --gz \
 	--out ${step2_file}

## remove variant inclusion list
rm ${tmpdir}/tmp.ex.chr${chr}_${names}.list
fi
## change to relevant directory
