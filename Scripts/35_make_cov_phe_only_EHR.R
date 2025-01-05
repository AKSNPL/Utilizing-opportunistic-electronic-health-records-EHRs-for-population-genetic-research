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

ehr_data<- arrow::read_feather( paste0("input/02_ehr_ukbbiom_combined/ehr_data", ".feather"))

# select the unique concept names
concept_names<- unique(ehr_data$concept_name)

#########################median##########################################################
# extract the phenotypes and save them
for (concept in concept_names) {
  cat(concept)
  data <- ehr_data %>%
    filter(concept_name== concept)
  
  data<- data%>%
    group_by(eid)%>%
    summarise(value_n= median(value_n))
  
  data <- data %>%
    rename(FID = eid,
           !!concept:= value_n)
  
  data$IID <- data$FID
  
  data <- data %>%
    select("FID","IID",paste0(concept))
  
  data <- data[complete.cases(data), ]
  
  fwrite(data,paste0("input/07_EHR_only_phe_cov/","phenotypes_",concept,".txt"), sep = "\t", row.names=F, quote=F, na = NA)
  cat("\n")
  
}

file<- fread(paste0("input/07_EHR_only_phe_cov/","phenotypes_","ALT",".txt"))

file1<- fread(paste0("input/07_EHR_only_phe_cov/","phenotypes_","ALP",".txt"))


file2<- merge(file,file1,by=c(("FID"= "FID"),("IID" = "IID")), all.x= TRUE, all.y= TRUE)

files<-data.frame(FID=NA,IID=NA)
for (concept in concept_names){
  file<- fread(paste0("input/07_EHR_only_phe_cov/","phenotypes_",concept,".txt"))
  files<- merge(file,files,by=c(("FID"= "FID"),("IID" = "IID")), all.x= TRUE, all.y= TRUE)
}
files<- files[-1,]
fwrite(files,paste0("input/07_EHR_only_phe_cov/","phenotypes_all",".txt"), sep = "\t", row.names=F, quote=F, na = NA)


#####################################################################################
## Median ####

folder_path <- "input/02_ehr_ukbbiom_combined/"
ehr_data <- arrow::read_feather(paste0(folder_path,"ehr_data.feather"))
ehr_data <- ehr_data %>%
  distinct(eid, value_n, date_EHR,concept_name, .keep_all = TRUE)

concept_names<- unique(ehr_data$concept_name)
for(concpt in concept_names){
  if(concpt=="Height" | concpt=="Haematocritperc"){
    ehr_data$value_n[ehr_data$concept_name ==concpt] <- ehr_data$value_n[ehr_data$concept_name ==concpt]*100
  }
  
}


ehr_median_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_median_value = median(value_n)) %>%
  ungroup()


arrow::write_feather(ehr_median_df, paste0("input/07_EHR_only_phe_cov/EHR_only_median", ".feather"), compression = "zstd")


## Mean ####
folder_path <- "input/02_ehr_ukbbiom_combined/"
ehr_data <- arrow::read_feather(paste0(folder_path,"ehr_data.feather"))
ehr_data <- ehr_data %>%
  distinct(eid, value_n, date_EHR,concept_name, .keep_all = TRUE)

concept_names<- unique(ehr_data$concept_name)
for(concpt in concept_names){
  if(concpt=="Height" | concpt=="Haematocritperc"){
    ehr_data$value_n[ehr_data$concept_name ==concpt] <- ehr_data$value_n[ehr_data$concept_name ==concpt]*100
  }
  
}


ehr_mean_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_mean_value = mean(value_n)) %>%
  ungroup()


arrow::write_feather(ehr_mean_df, paste0("input/07_EHR_only_phe_cov/EHR_only_mean", ".feather"), compression = "zstd")


## Minimum #### 
folder_path <- "input/02_ehr_ukbbiom_combined/"
ehr_data <- arrow::read_feather(paste0(folder_path,"ehr_data.feather"))
ehr_data <- ehr_data %>%
  distinct(eid, value_n, date_EHR,concept_name, .keep_all = TRUE)

concept_names<- unique(ehr_data$concept_name)
for(concpt in concept_names){
  if(concpt=="Height" | concpt=="Haematocritperc"){
    ehr_data$value_n[ehr_data$concept_name ==concpt] <- ehr_data$value_n[ehr_data$concept_name ==concpt]*100
  }
  
}


ehr_minimum_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_minimum_value = min(value_n)) %>%
  ungroup()


arrow::write_feather(ehr_minimum_df, paste0("input/07_EHR_only_phe_cov/EHR_only_minimum", ".feather"), compression = "zstd")
#####################
## Maximum ####
folder_path <- "input/02_ehr_ukbbiom_combined/"
ehr_data <- arrow::read_feather(paste0(folder_path,"ehr_data.feather"))
ehr_data <- ehr_data %>%
  distinct(eid, value_n, date_EHR,concept_name, .keep_all = TRUE)

concept_names<- unique(ehr_data$concept_name)
for(concpt in concept_names){
  if(concpt=="Height" | concpt=="Haematocritperc"){
    ehr_data$value_n[ehr_data$concept_name ==concpt] <- ehr_data$value_n[ehr_data$concept_name ==concpt]*100
  }
  
}

ehr_maximum_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_maximum_value = max(value_n)) %>%
  ungroup()


arrow::write_feather(ehr_maximum_df, paste0("input/07_EHR_only_phe_cov/EHR_only_maximum", ".feather"), compression = "zstd")


##########combine_data########################
combined0_data<- merge(ehr_minimum_df,ehr_maximum_df, by= c("eid","concept_name"))

combined1_data<- merge(ehr_median_df,ehr_mean_df, by= c("eid","concept_name"))

combined_data<- merge(combined0_data,combined1_data,by= c("eid","concept_name"))


arrow::write_feather(combined_data, paste0("input/07_EHR_only_phe_cov/all_combined_EHR_only", ".feather"), compression = "zstd")

###############################
folder_path <- "input/07_EHR_only_phe_cov/"
mainD<- arrow::read_feather(paste0(folder_path,"all_combined_EHR_only.feather"))

concept_names<- unique(mainD$concept_name)

# extract the phenotypes and save them
for (concept in concept_names) {
  cat(concept)
  data <- mainD %>%
    filter(concept_name== concept)
  
  data <- data %>%
    rename(FID = eid)
  
  data$IID <- data$FID
  
  data <- data %>%
    select("FID","IID",everything())
  
  
  data <- data %>%
    rename(!!paste0(concept, "_ehr_median") := ehr_median_value,
           !!paste0(concept, "_ehr_mean") := ehr_mean_value,
           !!paste0(concept, "_ehr_max") := ehr_maximum_value,
           !!paste0(concept, "_ehr_min") := ehr_minimum_value)
  
  data$concept_name<- NULL
  
  data <- data[complete.cases(data), ]
  
  fwrite(data,paste0("input/07_EHR_only_phe_cov/","phenotypes_",concept,".txt"), sep = "\t", row.names=F, quote=F, na = NA)
  cat("\n")
  
}

