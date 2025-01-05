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


# Iterate through concept_names
for (concept_name in concept_names) {
    # List all files in the directory that match the pattern
    file_list <- list.files(path = "output/17_GWAS_step2_UKB_only/",pattern = paste0("GWAS_chr\\d+_", concept_name, "_", concept_name,".regenie.gz"))
    
    if(length(file_list)==22){
      print(length(file_list))
      # Read REGENIE files and filter
      filtered_data <- data.frame()
      for (file in file_list) {
        cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
        file<- paste0("output/17_GWAS_step2_UKB_only/",file)
        data <- fread(file, fill = TRUE)
        filtered_data <- rbind(filtered_data, data[data$INFO >= 0.4 & data$SE <= 10 & data$A1FREQ >= 0.01 & data$A1FREQ <= 0.99, ])
        cat("Done\n")
      }
      
      # Save the filtered data to a new file
      fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/20_combined_allchrs_UKB_only/",concept_name, "_allchrs_", ".tsv.gz"))
      
      # Remove the original files
      #file_paths <- file.path(getwd(), file_list)
      #unlink(file_paths)
      #print(file_paths)
    }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}




###new

# Define a vector of concept_names
concept_names <- c("ALP", "ALT", "Albumin", "Basophills", "CRP",
                   "Calcium", "Cholesterol", "Creatinine", "DBP", "Eosinophills", "FEV1",
                   "FVC", "Glucose", "HDL", "Haematocritperc", "Haemoglobinconc", "HbA1c",
                   "Lymphocytes", "MCHbconc", "MCV", "Monocytes", "Neutrophills", "Platelets",
                   "RBC", "SBP", "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides"
)
#height:66056099
#alt:   65610373

# Iterate through concept_names
for (concept_name in concept_names) {
  # List all files in the directory that match the pattern
  file_list <- list.files(path = "output/17_GWAS_step2_UKB_only/",pattern = paste0("GWAS_chr\\d+_", concept_name, "_", concept_name,".regenie.gz"))
  
  if(length(file_list)==22){
    print(length(file_list))
    # Read REGENIE files and filter
    filtered_data <- data.frame()
    for (file in file_list) {
      cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
      file<- paste0("output/17_GWAS_step2_UKB_only/",file)
      data <- fread(file, fill = TRUE)
      filtered_data <- rbind(filtered_data, data)
      cat("Done\n")
    }
    
    # Save the filtered data to a new file
    fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/20_combined_allchrs_UKB_only/",concept_name, "_allchrs_new", ".tsv.gz"))
    
    # Remove the original files
    #file_paths <- file.path(getwd(), file_list)
    #unlink(file_paths)
    #print(file_paths)
  }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}

###new combine:
rm(list = ls())
gc()
setwd("/Aakash/12_GWAS_EHR_UKB/")
# Load required libraries
library(dplyr)
library(readr)
library(data.table)
# Get the current working directory
current_dir <- getwd()

# Define a vector of concept_names
concept_names <- c("Height")
concept_name <- "Height"
#height:66056099
#alt:   65610373
#new: 9764427

# Iterate through concept_names
for (concept_name in concept_names) {
  # List all files in the directory that match the pattern
  file_list <- list.files(path = "output/17_GWAS_step2_UKB_only/",pattern = paste0("GWAS_chr\\d+_", concept_name, "_", concept_name,".regenie.gz"))
  
  if(length(file_list)==22){
    print(length(file_list))
    # Read REGENIE files and filter
    filtered_data <- data.frame()
    for (file in file_list) {
      cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
      file<- paste0("output/17_GWAS_step2_UKB_only/",file)
      data <- fread(file, fill = TRUE)
      filtered_data <- rbind(filtered_data, data)
      cat("Done\n")
    }
    
    # Save the filtered data to a new file
    #fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/20_combined_allchrs_UKB_only/",concept_name, "_allchrs_new", ".tsv.gz"))
    
    # Remove the original files
    #file_paths <- file.path(getwd(), file_list)
    #unlink(file_paths)
    #print(file_paths)
  }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}
# Make the ID column unique
filtered_data <- filtered_data %>%
  distinct(ID, .keep_all = TRUE)
#66056099
# Save the combined data to a new file 
fwrite(filtered_data, "output/17_GWAS_step2_UKB_only/Height_allchrs_combined_unique.tsv.gz", sep = "\t", compress = "gzip")

# List of file names
file_names <- c(
  "f1to5_allchrs_combined_unique.tsv.gz",
  "f6to10_allchrs_combined_unique.tsv.gz",
  "f11to15_allchrs_combined_unique.tsv.gz",
  "f16to20_allchrs_combined_unique.tsv.gz",
  "f21to25_allchrs_combined_unique.tsv.gz",
  "f26to31_allchrs_combined_unique.tsv.gz"
)

# Read in all files and combine them into one data frame
combined_data <- file_names %>%
  lapply(read_tsv) %>%
  bind_rows()

# Make the ID column unique
combined_data <- combined_data %>%
  distinct(ID, .keep_all = TRUE)

# Save the combined data to a new file
fwrite(combined_data, "allchrs_combined_unique_o_ukb.tsv.gz", sep = "\t", compress = "gzip")

file_names <- c(
  "/Aakash/12_GWAS_EHR_UKB/output/20_combined_allchrs_UKB_only/allchrs_combined_unique_o_ukb.tsv.gz",
  "/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/allchrs_combined_unique.tsv.gz"
)

# Read in all files and combine them into one data frame
combined_data <- file_names %>%
  lapply(read_tsv) %>%
  bind_rows()

# Make the ID column unique
combined_data <- combined_data %>%
  distinct(ID, .keep_all = TRUE)

combined_data <- combined_data %>%
  select(!EXTRA)

combined_data <- combined_data %>%
  drop_na()

# Save the combined data to a new file
fwrite(combined_data, "allchrs_combined_uniqueplus.tsv.gz", sep = "\t", compress = "gzip")

####
aa<- fread("/Aakash/12_GWAS_EHR_UKB/output/20_combined_allchrs_UKB_only/Basophills_allchrs_.tsv.gz")
