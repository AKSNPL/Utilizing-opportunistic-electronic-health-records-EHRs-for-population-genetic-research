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

##################################
#### import basic information ####
##################################

## An example list of columns:
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.20116.0.0", "f.1558.0.0", paste0("f.22009.0.", 1:10), paste0("f.20003.0.", 1:47))

## names to assign
cl.names        <- c("f.eid", "age", "sex", "baseline_date", "centre", "bmi", "smoking", "alcohol", paste0("pc", 1:10), paste0("self_med", 1:47))

## import data from the main release
ukb.dat         <- read_parquet("/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
                                col_select = cl.select)

#select the unique eids in ukb data
ukb.dat1 <- data.frame(eid=unique(ukb.dat$f.eid))

#read all the data from EHR_UKB
ttl_calc <- arrow::read_feather("input/02_ehr_ukbbiom_combined/all_combined_EHR_UKB.feather")

# select all the unique phenotypes
concept_names <- unique(ttl_calc$concept_name)

# select all unique eids in qced ehr
uniq_EHR<- unique(ttl_calc$eid)

ttl_impli1<- data.frame(ehr=ukb.dat1$eid%in%uniq_EHR,eid=ukb.dat1$eid)

concept_names<- concept_names

all_resu<- list()
for(concept in concept_names){
  ttl_impli<- data.frame()
  #select a concept e.g. weight
  ttl_calc_w <- ttl_calc %>%
    filter(concept_name==concept)
  
  # select the unique eids
  ttl_calc_uniq<- unique(ttl_calc_w$eid)
  
  names<- paste0("bin_",concept,"_ehr")
  
  # check for the indications
  ttl_impli <- data.frame(names = (uniq_EHR%in%ttl_calc_uniq)) %>%
    rename(!!names:= names)
  
  ttl_impli<- as.integer(as.logical(ttl_impli[,1]))
  all_resu[[concept]]<- ttl_impli
}

all_resu <- as.data.frame(all_resu)

all_resu$eid<-uniq_EHR

ttl_im<- merge(ttl_impli1,all_resu, by= c("eid"),all.x = TRUE)

ttl_im$ehr <- as.integer(as.logical(ttl_im[,2]))

pheno1<- ttl_im%>%
  select(eid,ehr)

pheno1$FID <- pheno1$eid
pheno1$IID <- pheno1$eid
pheno1$eid<- NULL
pheno1<- pheno1 %>%
  select(FID,IID,ehr)


pheno2<- ttl_im %>%
  select(!ehr) 
pheno2$FID <- pheno2$eid
pheno2$IID <- pheno2$eid
pheno2$eid<- NULL

pheno2<- pheno2 %>%
  select(FID,IID,everything())

fwrite(pheno1,paste0("input/05_indications_phe_cov/","phenotypes_ehr.txt"), sep = "\t", row.names=F, quote=F, na = NA)

fwrite(pheno2,paste0("input/05_indications_phe_cov/","phenotypes_all.txt"), sep = "\t", row.names=F, quote=F, na = NA)

all_resu<- list()
for(concept in concept_names){
  #select a concept e.g. weight
  ttl_calc_w <- ttl_calc %>%
    filter(concept_name==concept)
  
  names<- paste0("bin_",concept,"_ehr")
  
  # select the unique eids
  ttl_calc_uniq<- unique(ttl_calc_w$eid)
  
  # check for the indications
  ttl_impli <- data.frame(eid=uniq_EHR,bin_weight_ehr= uniq_EHR%in%ttl_calc_uniq)
  
  ttl_im<- merge(ttl_impli1,ttl_impli, by= c("eid"),all.x = TRUE)
  
  ttl_imp<- as.integer(as.logical(ttl_im[,3]))
  all_resu[[concept]]<- ttl_imp
}

all_resu <- as.data.frame(all_resu)


