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

## EHR ####
# Define a vector of concept_names
concept_names <- c("ALP", "ALT", "Albumin", "Basophills", "CRP",
                   "Calcium", "Cholesterol", "Creatinine", "DBP", "Eosinophills", "FEV1",
                   "FVC", "Glucose", "HDL", "Haematocritperc", "Haemoglobinconc", "HbA1c",
                   "Lymphocytes", "MCHbconc", "MCV", "Monocytes", "Neutrophills", "Platelets",
                   "RBC", "SBP", "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides"
)
#concept_name<- "ALP"

# Define a list of 'minmaxmid"
minmaxmid_values <- c("all")


# Iterate through concept_names
for (concept_name in concept_names) {
    # List all files in the directory that match the pattern
    file_list <- list.files(path = "output/15_GWAS_step2_ind_all/",pattern = paste0("GWAS_chr\\d+_", "all_", concept_name, ".regenie.gz"))
    
    if(length(file_list)==22){
      print(length(file_list))
      # Read REGENIE files and filter
      filtered_data <- data.frame()
      for (file in file_list) {
        cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
        file<- paste0("output/15_GWAS_step2_ind_all/",file)
        data <- fread(file, fill = TRUE)
        filtered_data <- rbind(filtered_data, data[data$INFO >= 0.4 & data$SE <= 10 & data$A1FREQ >= 0.01 & data$A1FREQ <= 0.99, ])
        cat("Done\n")
      }
      
      # Save the filtered data to a new file
      fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/19_combined_allchrs_all_pheno/",concept_name, "_allchrs_", ".tsv.gz"))
      
      # Remove the original files
      #file_paths <- file.path(getwd(), file_list)
      #unlink(file_paths)
      #print(file_paths)
    }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}
