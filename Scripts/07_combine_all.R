# Clear the workspace and garbage collect
rm(list = ls()) # Remove all objects from the workspace
ay <- gc() # Trigger garbage collection to free up memory

# Set display options
options(scipen = 5) # Prevent scientific notation for numbers

# Load session information for reproducibility
sessioninfo::session_info()

# Set the working directory to the specified path
setwd("/Aakash/12_GWAS_EHR_UKB")

# Load required libraries
for (i in c(
  "data.table", "lubridate", "tidyr", "ggplot2", "dplyr", 
  "stringr", "arrow", "readxl", "Hmisc", "tidyverse", "here"
)) {
  suppressPackageStartupMessages(
    library(i, character.only = TRUE) # Load libraries without showing messages
  )
}

## Set root and parallel processing settings ####
# Configure Knitr to use the project root for paths
knitr::opts_knit$set(
  root.dir = here()
)

# Set the number of threads for data.table operations
writeLines(paste0("Threads available: ", parallel::detectCores())) # Show total cores available
writeLines(paste0("Threads given to data.table: ", parallel::detectCores() / 2)) # Show cores allocated
setDTthreads(parallel::detectCores() / 2) # Allocate half the available threads to data.table

# Define the folder path where input data is located
folder_path <- "input/02_ehr_ukbbiom_combined/"

# Load the dataset containing minimum date distances
min_date_data <- arrow::read_feather(paste0(folder_path, "combined_EHR_UKB_min_date.feather"))

# Rename columns for clarity
min_date_data <- min_date_data %>%
  rename(
    ehr_min_dt = value_n,          # Rename EHR value column
    ukb_min_dt = bm_value,         # Rename UKB value column
    concept_name = concept_name_ehr # Rename concept name column
  )

# Load datasets for minimum, maximum, median, and mean values
minimum_data <- arrow::read_feather(paste0(folder_path, "combined_EHR_UKB_minimum.feather"))
maximum_data <- arrow::read_feather(paste0(folder_path, "combined_EHR_UKB_maximum.feather"))
median_data <- arrow::read_feather(paste0(folder_path, "combined_EHR_UKB_median.feather"))
mean_data <- arrow::read_feather(paste0(folder_path, "combined_EHR_UKB_mean.feather"))

# Merge the minimum and maximum datasets by `eid` and `concept_name`
combined0_data <- merge(minimum_data, maximum_data, by = c("eid", "concept_name"))

# Merge the median and mean datasets by `eid` and `concept_name`
combined1_data <- merge(median_data, mean_data, by = c("eid", "concept_name"))

# Merge the combined datasets with minimum/maximum and median/mean data
combined2_data <- merge(combined0_data, combined1_data, by = c("eid", "concept_name"))

# Merge the result with the minimum date distance data
combined_data <- merge(combined2_data, min_date_data, by = c("eid", "concept_name"))

# Select and rename relevant columns for the final dataset
combined_data <- combined_data %>%
  select(
    eid, concept_name, ehr_min_dt, ehr_minimum_value, ehr_maximum_value, 
    ehr_median_value, ehr_mean_value, bm_minimum_value
  ) %>%
  rename(ukb_value = bm_minimum_value) # Rename UKB value column for clarity

# Save the combined dataset as a Feather file with compression
arrow::write_feather(
  combined_data, 
  paste0("input/02_ehr_ukbbiom_combined/all_combined_EHR_UKB", ".feather"), 
  compression = "zstd"
)
