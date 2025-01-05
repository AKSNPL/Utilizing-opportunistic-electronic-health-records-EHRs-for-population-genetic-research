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

## make the phenotypes_files: ####
folder_path <- "input/02_ehr_ukbbiom_combined/"
data0 <- arrow::read_feather(paste0(folder_path,"all_combined_EHR_UKB.feather"))
concept_names<- unique(data0$concept_name)
for (concept in concept_names) {
  cat(concept)
  data <- data0 %>%
    filter(concept_name== concept)
  
  data <- data %>%
    rename(!!paste0(concept, "_ehr_md") := ehr_min_dt,
           !!paste0(concept, "_ukb_md") := ukb_value,
           !!paste0(concept, "_ehr_median") := ehr_median_value,
           !!paste0(concept, "_ehr_mean") := ehr_mean_value,
           !!paste0(concept, "_ehr_max") := ehr_maximum_value,
           !!paste0(concept, "_ehr_min") := ehr_minimum_value,
           FID = eid)
  
  data$IID <- data$FID
  
  data <- data %>%
    select("FID","IID",paste0(concept,"_ehr_md"),paste0(concept,"_ukb_md"),paste0(concept, "_ehr_median")
           ,paste0(concept, "_ehr_mean"),paste0(concept, "_ehr_max"),paste0(concept, "_ehr_min") )
  
  data <- data[complete.cases(data), ]
  
  fwrite(data,paste0("input/03_pheno_cov_files/","phenotypes_",concept,".txt"), sep = "\t", row.names=F, quote=F, na = NA)
  cat("\n")

}


## covariate file: ###################################################

### 0. set paths of folders ##########################################
dir.main <- c('/.')
dir.ukbb <- paste0(dir.main,'/data/UK_biobank')
dir.people <- paste0(dir.main, '/people')
dir.aakash <- paste0(dir.people, '/Aakash')
current.proj <- paste0(dir.aakash, '/12_GWAS_EHR_UKB')


### 1. Set working directory in the project folder #####################
setwd(current.proj)


### 2. Load packages, if they are not installed => install.package('arrow') #######################################
### 3. create file names (fn) with address to read or write  easily #########################################################

ukbb_dict_fn <- paste0(dir.ukbb, '/phenotypes/working_data/Data.dictionary.UKBB.txt')
ukbb_data_fn <- paste0(dir.ukbb, "/phenotypes/working_data/parquet_files/ukb.parquet")
ukbb_EUR_fn <- paste0("/Aakash/11_GWAS_EHR_W_H/input/EUR_panukbb_regenie_format.id")#/regenie_44448/EUR_panukbb_regenie_format.id



#### 4. read files: Metadata information  ################################################################################################
ukbb_dict_tbl <- read.table(ukbb_dict_fn, sep="\t", header=T)
ukbb_EUR <- read.table(ukbb_EUR_fn, header=F)


#### 5. Get data from ukb main table ###########################################

## Columns to read from phenotype data (Refer to ukbb_dict_tbl below for id:label mapping)
col_ids <- c("f.eid","f.53.0.0", "f.21022.0.0", "f.31.0.0", "f.54.0.0", "f.22000.0.0", paste0("f.22009.0.", 1:20), "f.21001.0.0", "f.50.0.0")
cols_labels <- c("f.eid","baseline_date", "age", "sex", "centre", "batch", paste0("pc", 1:20), "bmi", "height")


## Read only the columns of interest (col_ids) from the big data file for UKBB 'ukb45268.parquet' (baseline charateristics, PCs)
ukbb_data_tbl <- arrow::read_parquet(ukbb_data_fn, col_select = all_of(col_ids )) 

## check how many unique individuals are in ukbb_data_tbl, indivdiuals are anonmyised and they are IDs are 'f.eid' column.
length(unique(ukbb_data_tbl$f.eid)) #502,490

## Rename the columns of interest in vector 'col_ids'
names(ukbb_data_tbl) <- cols_labels

## Convert to data table
ukbb_data_tbl <- as.data.table(ukbb_data_tbl)

#### 6. Only EUR ###########################################
## Only keep EURopeans in ukbb_data_tbl
ukbb_eur_tbl <- ukbb_data_tbl[which(ukbb_data_tbl$f.eid%in%ukbb_EUR$V1), ] 
## remove ukbb_data_tbl
rm(ukbb_data_tbl)

#### 7. some additional processing on variables in ukbb_eur_tbl  ###########################################

##  Transform batch to binary
ukbb_eur_tbl$batch  <- ifelse(ukbb_eur_tbl$batch < 0, 1, 0)
#table(ukbb.dat$batch)

## Check non-nulls
colSums(!is.na(ukbb_eur_tbl))
# f.eid    age    sex centre  batch    pc1 ... pc20    bmi height    
# 502490 502489 502489 502489 488248 488248 ... 488248 499385 499951 

## Add FID (means family ID) and IID (individual ID) columns, both are same. 
ukbb_eur_tbl$FID <- ukbb_eur_tbl$f.eid
ukbb_eur_tbl$IID <- ukbb_eur_tbl$f.eid

ukb.cov <- ukbb_eur_tbl[, c("FID", "IID", "age", "sex", "centre", "batch", paste0("pc", 1:20)), with = F]

## drop participants with missing values in covariates
ukb.cov         <- na.omit(ukb.cov)


## make biological sex binary
ukb.cov$sex     <- ifelse(ukb.cov$sex == "Female", 1, 0)

#saving
ukb.cov.fn <- paste0(current.proj, '/input/03_pheno_cov_files/covariates.txt')


fwrite(ukb.cov, ukb.cov.fn , sep = "\t", row.names=F, quote=F)

