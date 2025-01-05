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

## import the required data

ukb_data<- arrow::read_feather( paste0("input/02_ehr_ukbbiom_combined/ukb_data", ".feather"))

# select the unique concept names
concept_names<- unique(ukb_data$concept_name)

# extract the phenotypes and save them
for (concept in concept_names) {
  cat(concept)
  data <- ukb_data %>%
    filter(concept_name== concept)
  
  data <- data %>%
    rename(FID = eid,
           !!concept:= bm_value)
  
  data$IID <- data$FID
  
  data <- data %>%
    select("FID","IID",paste0(concept))
  
  data <- data[complete.cases(data), ]
  
  fwrite(data,paste0("input/06_UKB_only_phe_cov/","phenotypes_",concept,".txt"), sep = "\t", row.names=F, quote=F, na = NA)
  cat("\n")
  
}

file<- fread(paste0("input/06_UKB_only_phe_cov/","phenotypes_","ALT",".txt"))

file1<- fread(paste0("input/06_UKB_only_phe_cov/","phenotypes_","ALP",".txt"))


file2<- merge(file,file1,by=c(("FID"= "FID"),("IID" = "IID")), all.x= TRUE, all.y= TRUE)

files<-data.frame(FID=NA,IID=NA)
for (concept in concept_names){
  file<- fread(paste0("input/06_UKB_only_phe_cov/","phenotypes_",concept,".txt"))
  files<- merge(file,files,by=c(("FID"= "FID"),("IID" = "IID")), all.x= TRUE, all.y= TRUE)
}
files<- files[-1,]
fwrite(files,paste0("input/06_UKB_only_phe_cov/","phenotypes_all",".txt"), sep = "\t", row.names=F, quote=F, na = NA)



