### create subset 30Kfor LD clumping################
####################################################

rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/12_GWAS_EHR_UKB")
for (i in c(
  "data.table",
  "lubridate",
  "tidyr",
  "ggplot2",
  "dplyr",
  "stringr",
  "arrow",
  "readxl",
  "Hmisc",
  "tidyverse",
  "here"
)
) {
  suppressPackageStartupMessages(
    library(i, character.only = TRUE
    ))
}

## set root and parallel settings ####
# Knitr should use the project root and not the script location as root
# base.dir refers to the plot location, that should remain with the script
knitr::opts_knit$set(
  root.dir = here()
)

# Give data.table enough threads
writeLines(paste0("Threads available: ", parallel::detectCores()))
writeLines(paste0("Threads given to data.table: ", parallel::detectCores() / 2))
setDTthreads(parallel::detectCores() / 2)

# load the required data ####
new <- read.table("/imputed/bgen_files_44448/ukb22828_c10_b0_v3.sample", header=T, sep = " ")

#length(intersect(old$ID_1, new$ID_1))  114

# also load withrawals
withr <- read.table("/Aakash/12_GWAS_EHR_UKB/input/05_indications_phe_cov/EUR_panukbb_regenie_format.id")

## now extract the withdrawal  and randomly select 30k samples

# load withdr freshly.. draw random 30 samples, do the rest as below given withr now has 30k samples

set.seed(12345)
sub_30 <- withr[sample(nrow(withr), 30000), "V1"] %>% as.data.frame() %>% dplyr::rename(., V1=.)

##### the format for the soft (ldsc?) requires 0s in the first row
new_row <- c(V1 = 0)
sub_30 <- rbind(new_row,sub_30)
length(sort(sub_30$V1))
dim(sub_30)

#now I need to select sub_30 from new ID1 column. It already includes 0


sample_30k <- new %>% filter(new$ID_1 %in% sub_30$V1)
summary(as.factor(sample_30k$sex))# 13.5k 1 vs 16.3k 2

vc<- intersect(sample_30k$ID_1,new$ID_1)
fwrite(sample_30k, "/Aakash/12_GWAS_EHR_UKB/input/04_others/samples_30k.sample", sep = " ", row.names=F, quote=F)

test<- fread("/data/UK_biobank/genotypes/variant_qc/input/UKBB.samples.random.30k.20211119.sample")

jj<- merge(sample_30k,test, by = "ID_1")