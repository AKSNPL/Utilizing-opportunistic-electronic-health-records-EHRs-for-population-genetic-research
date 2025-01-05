rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("Aakash/12_GWAS_EHR_UKB")
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

# Define a list of 'minmaxmid"
minmaxmid_values <- c("min", "max", "median", "mean", "md")


# Iterate through concept_names
for (concept_name in concept_names) {
  # Iterate through 'minmaxmid' values
  for (minmaxmid in minmaxmid_values) {
    # List all files in the directory that match the pattern
    file_list <- list.files(path = "output/06_output_step2/",pattern = paste0("GWAS_chr\\d+_", concept_name, "_", concept_name, "_ehr_", minmaxmid, ".regenie.gz"))
    
    if(length(file_list)==22){
      print(length(file_list))
      # Read REGENIE files and filter
      filtered_data <- data.frame()
      for (file in file_list) {
        cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
        file<- paste0("output/06_output_step2/",file)
        data <- fread(file, fill = TRUE)
        filtered_data <- rbind(filtered_data, data[data$INFO >= 0.4 & data$SE <= 10 & data$A1FREQ >= 0.01 & data$A1FREQ <= 0.99, ])
        cat("Done\n")
      }
      
      # Save the filtered data to a new file
      fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/07_combined_allchrs/",concept_name, "_allchrs_", minmaxmid, ".tsv.gz"))
      
      # Remove the original files
      #file_paths <- file.path(getwd(), file_list)
      #unlink(file_paths)
      #print(file_paths)
    }
  }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}

## UKBB ####
concept_names <- c("ALP", "ALT", "Albumin", "Basophills", "CRP",
"Calcium", "Cholesterol", "Creatinine", "DBP", "Eosinophills", "FEV1",
"FVC", "Glucose", "HDL", "Haematocritperc", "Haemoglobinconc", "HbA1c",
"Lymphocytes", "MCHbconc", "MCV", "Monocytes", "Neutrophills", "Platelets",
"RBC", "SBP", "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides"
)


# Define a list of 'minmaxmid' values

minmaxmid_values<- "_ukb_md"

# Iterate through concept_names
for (concept_name in concept_names) {
  # Iterate through 'minmaxmid' values
  for (minmaxmid in minmaxmid_values) {
    # List all files in the directory that match the pattern
    file_list <- list.files(path = "output/06_output_step2/", pattern = paste0("GWAS_chr\\d+_", concept_name, "_", concept_name, minmaxmid, ".regenie.gz"))
    
    if(length(file_list)==22){
      print(length(file_list))
      # Read REGENIE files and filter
      filtered_data <- data.frame()
      for (file in file_list) {
        cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
        file<- paste0("output/06_output_step2/",file)
        data <- fread(file, fill = TRUE)
        filtered_data <- rbind(filtered_data, data[data$INFO >= 0.4 & data$SE <= 10 & data$A1FREQ >= 0.01 & data$A1FREQ <= 0.99, ])
        cat("Done\n")
      }
      
      # Optionally add your calculations here (e.g., LOG10P calculation)
      
      # Save the filtered data to a new file
      fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/07_combined_allchrs/",concept_name, "_allchrs_", minmaxmid, ".tsv.gz"))
      
      # Remove the original files
      #file_paths <- file.path(getwd(), file_list)
      #unlink(file_paths)
      #print(file_paths)
    }
  }
  cat("done:",concept_name)
  cat(" : ",length(file_list),"\n")
}

###a new check ####
# Load necessary packages
library(data.table)
library(dplyr)
library(purrr)

# Define a vector of concept_names
concept_names <- c("ALP", "ALT", "Albumin", "Basophills", "CRP",
                   "Calcium", "Cholesterol", "Creatinine", "DBP", "Eosinophills", "FEV1",
                   "FVC", "Glucose", "HDL", "Haematocritperc", "Haemoglobinconc", "HbA1c",
                   "Lymphocytes", "MCHbconc", "MCV", "Monocytes", "Neutrophills", "Platelets",
                   "RBC", "SBP", "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides"
)

# Define lists of 'minmaxmid' values
ehr_minmaxmid_values <- c("mean")
ukbb_minmaxmid_values <- c("")

# Combine both minmaxmid value lists
all_minmaxmid_values <- list(ehr = ehr_minmaxmid_values, ukbb = ukbb_minmaxmid_values)

# Iterate through concept_names
for (concept_name in concept_names) {
  # Iterate through all minmaxmid value types
  for (type in names(all_minmaxmid_values)) {
    minmaxmid_values <- all_minmaxmid_values[[type]]
    
    for (minmaxmid in minmaxmid_values) {
      # Define the pattern based on the type
      if (type == "ehr") {
        pattern <- paste0("GWAS_chr\\d+_", concept_name, "_", concept_name, "_ehr_", minmaxmid, ".regenie.gz")
      } else if (type == "ukbb") {
        pattern <- paste0("GWAS_chr\\d+_", concept_name, "_", concept_name, minmaxmid, ".regenie.gz")
      }
      
      # List all files in the directory that match the pattern
      file_list <- list.files(path = "output/06_output_step2/", pattern = pattern)
      
      if (length(file_list) == 22) {
        print(length(file_list))
        # Read REGENIE files and filter
        filtered_data <- data.frame()
        for (file in file_list) {
          cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
          file <- paste0("output/06_output_step2/", file)
          data <- fread(file, fill = TRUE)
          filtered_data <- rbind(filtered_data, data)
          cat("Done\n")
        }
        
        # Save the filtered data to a new file
        fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/07_combined_allchrs/", concept_name, "_allchrs_new_", minmaxmid, ".tsv.gz"))
        
        # Remove the original files if necessary
        # file_paths <- file.path(getwd(), file_list)
        # unlink(file_paths)
        # print(file_paths)
      }
    }
  }
  cat("done:", concept_name)
  cat(" : ", length(file_list), "\n")
}




####
setwd("/people/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs")
# Load required libraries
library(dplyr)
library(readr)
library(data.table)
# Get the current working directory
current_dir <- getwd()

# List all files in the current directory
all_files <- list.files(current_dir)

# Filter files based on the pattern
filtered_files <- all_files[grepl("_allchrs_new_mean.tsv.gz$", all_files)][26:31]

# Read in all files and combine them into one data frame
combined_data <- filtered_files %>%
  lapply(read_tsv) %>%
  bind_rows()

# Make the ID column unique
combined_data <- combined_data %>%
  distinct(ID, .keep_all = TRUE)

head(combined_data) 
# Display the combined data
print(combined_data)

# Save the combined data to a new file
fwrite(combined_data, "f26to31_allchrs_combined_unique.tsv.gz", sep = "\t", compress = "gzip")
#1-5: 56642629
#6-10 : 57482105
#11-15: 56411764
#16-20: 57232811
#21-25: 56496776
#26-31: 57835459
#ukb:   66065638

setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs")
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
fwrite(combined_data, "allchrs_combined_unique.tsv.gz", sep = "\t", compress = "gzip")

#check:
aa<- fread("Height_allchrs_all_ukb_ukb_md.tsv.gz")

com<- combined_data%>%
  filter(GENPOS == "141148854")


############test
setwd("/sc-projects/sc-proj-computational-medicine/people/Aakash/12_GWAS_EHR_UKB/")
# Load necessary packages
library(data.table)
library(dplyr)
library(purrr)

# Define a vector of concept_names
concept_names <- c('Height')


# Define lists of 'minmaxmid' values
# Define a list of 'minmaxmid"
minmaxmid_values <- c("")

ukbb_minmaxmid_values <- c("1")

# Combine both minmaxmid value lists
all_minmaxmid_values <- list(ukbb = ukbb_minmaxmid_values)

# Iterate through concept_names
for (concept_name in concept_names) {
      # Define the pattern based on the type
      pattern <- paste0("GWAS_chr\\d+_", concept_name, "_", concept_name, ".regenie.gz")
      
      # List all files in the directory that match the pattern
      file_list <- list.files(path = "output/8999_test_only_ukb_step2/", pattern = pattern)
      
      if (length(file_list) == 22) {
        print(length(file_list))
        # Read REGENIE files and filter
        filtered_data <- data.frame()
        for (file in file_list) {
          cat("Filtering: ", file, " for low quality MAF(A1FREQ) ...\n")
          file <- paste0("output/8999_test_only_ukb_step2/", file)
          data <- fread(file, fill = TRUE)
          filtered_data <- rbind(filtered_data, data)
          cat("Done\n")
        }
        
        # Save the filtered data to a new file
        fwrite(filtered_data, quote = FALSE, sep = "\t", row.names = FALSE, compress = "gzip", file = paste0("output/8999_test_only_ukb_step2/", concept_name, "_allchrs_new_", minmaxmid, ".tsv.gz"))
        
        # Remove the original files if necessary
        # file_paths <- file.path(getwd(), file_list)
        # unlink(file_paths)
        # print(file_paths)
      }
  cat("done:", concept_name)
  cat(" : ", length(file_list), "\n")
}
