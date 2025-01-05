#!/bin/sh

## script to obtain SNP dosages to create an LD-matrix

## export location of files
export dir=/imputed/bgen_files

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}
export pheno=${4}

echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

if [ ${chr} -eq 23 ]; then

## create subset bgen file get only SNPs in the data set
/programs/bgen/build/apps/bgenix \
-g ${dir}/ukb22828_cX_b0_v3.bgen \
-incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen

## create dosage file (restrict to random subset of 100k individuals to be computationally efficient but with enough precision for SuSiE)
/programs/qctool/build/release/qctool_v2.0.7 \
-g tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen \
-s ${dir}/ukb22828_cX_b0_v3.sample \
-incl-samples /Aakash/12_GWAS_EHR_UKB/input/04_others/samples_30k.sample \
-og - \
-ofiletype dosage > tmpdir/tmp.${pheno}.${chr}.${lowpos}.${uppos}.dosage

else
  
## create subset bgen file get only SNPs in the data set
/programs/bgen/build/apps/bgenix \
-g ${dir}/ukb22828_c${chr}_b0_v3.bgen \
-incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen

## create dosage file (restrict to random subset of 100k individuals to be computationally efficient but with enough precision for SuSiE)
/programs/qctool/build/release/qctool_v2.0.7 \
-g tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen \
-s ${dir}/ukb22828_c${chr}_b0_v3.sample \
-incl-samples /Aakash/12_GWAS_EHR_UKB/input/04_others/samples_30k.sample \
-og - \
-ofiletype dosage > tmpdir/tmp.${pheno}.${chr}.${lowpos}.${uppos}.dosage

fi

rm tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen
rm tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst
