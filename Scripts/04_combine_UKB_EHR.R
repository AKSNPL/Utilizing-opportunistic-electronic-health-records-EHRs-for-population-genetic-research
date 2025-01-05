# Clear the workspace and garbage collect
rm(list = ls())  # Removes all objects from the workspace
ay <- gc()  # Triggers garbage collection and stores the result

# Set options
options(scipen = 5)  # Avoid scientific notation for numbers

# Load required packages
sessioninfo::session_info()  # Prints the current session information for reproducibility
setwd("/Aakash/12_GWAS_EHR_UKB")  # Set working directory

# Load necessary libraries dynamically
for (i in c(
  "data.table", "lubridate", "tidyr", "ggplot2", "dplyr", 
  "stringr", "arrow", "readxl", "Hmisc", "tidyverse", "here"
)) {
  suppressPackageStartupMessages(
    library(i, character.only = TRUE)  # Dynamically loads each package without startup messages
  )
}

# Set project root for knitr and adjust settings for parallel processing
knitr::opts_knit$set(root.dir = here())  # Sets project root as the working directory
writeLines(paste0("Threads available: ", parallel::detectCores()))  # Displays available CPU cores
writeLines(paste0("Threads given to data.table: ", parallel::detectCores() / 2))  # Allocates half of cores
setDTthreads(parallel::detectCores() / 2)  # Configures threads for data.table

# Data preparation: Load and process EHR and UKB data
folder_path <- "input/02_ehr_ukbbiom_combined/"  # Define folder path
ehr_data <- arrow::read_feather(paste0(folder_path, "ehr_data.feather")) %>%
  distinct(eid, value_n, date_EHR, concept_name, .keep_all = TRUE)  # Remove duplicates

ukb_data <- arrow::read_feather(paste0(folder_path, "ukb_data.feather")) %>%
  distinct(eid, bm_value, baseline_date, concept_name, .keep_all = TRUE)  # Remove duplicates

# Summarize data by concept_name
data_summary <- ukb_data %>%
  group_by(concept_name) %>%
  dplyr::summarise(
    median_bm = median(bm_value, na.rm = TRUE),
    min_bm = min(bm_value, na.rm = TRUE),
    max_bm = max(bm_value, na.rm = TRUE)
  )

data_summary.data <- ehr_data %>%
  group_by(concept_name) %>%
  dplyr::summarise(
    median_ehr = median(value_n, na.rm = TRUE),
    min_ehr = min(value_n, na.rm = TRUE),
    max_ehr = max(value_n, na.rm = TRUE)
  )

# Merge summaries from EHR and UKB
metadata <- data_summary %>%
  left_join(data_summary.data, by = "concept_name")

# Convert date columns to Date type
ukb_data$baseline_date <- as.Date(ukb_data$baseline_date)
ehr_data$date_EHR <- as.Date(ehr_data$date_EHR)

# Filter and process data by concept_name
concept_names <- unique(ehr_data$concept_name)
filtered_data_list <- list()  # Initialize list for storing filtered data

for (concept_nam in concept_names) {
  print(concept_nam)  # Display current concept name
  ehr_filtered <- ehr_data %>%
    filter(concept_name == concept_nam)  # Filter EHR data for current concept
  
  # Scale values for specific concepts
  if (concept_nam == "Height" | concept_nam == "Haematocritperc") {
    ehr_filtered$value_n <- ehr_filtered$value_n * 100
  }

  ukb_filtered <- ukb_data %>%
    filter(concept_name == concept_nam)  # Filter UKB data for current concept

  # Match EHR and UKB data by closest date
  matched_data <- ehr_filtered %>%
    inner_join(ukb_filtered, by = "eid", suffix = c("_ehr", "_ukb")) %>%
    mutate(
      date_distance = date_EHR - baseline_date,
      bm_value_new = ifelse(bm_value == 0, lag(na.omit(bm_value), default = bm_value[1]), bm_value),
      value_n_new = ifelse(value_n == 0, lag(na.omit(value_n), default = value_n[1]), value_n)
    ) %>%
    group_by(eid) %>%
    mutate(min_abs_date_distance = min(abs(date_distance))) %>%
    filter(abs(date_distance) == min_abs_date_distance) %>%
    arrange(eid, min_abs_date_distance, date_distance) %>%
    slice(1) %>%
    select(-min_abs_date_distance) %>%
    mutate(
      bm_value = ifelse(bm_value == 0, bm_value_new, bm_value),
      value_n = ifelse(value_n == 0, value_n_new, value_n)
    ) %>%
    select(-bm_value_new, -value_n_new)

  # Store matched data
  filtered_data_list[[concept_nam]] <- matched_data
}

# Combine all filtered data into a single dataframe
combined_data <- bind_rows(filtered_data_list)
arrow::write_feather(
  combined_data, 
  paste0(folder_path, "combined_EHR_UKB_min_date.feather"), 
  compression = "zstd"
)

arrow::write_feather(combined_data, paste0("input/02_ehr_ukbbiom_combined/combined_EHR_UKB_min_date", ".feather"), compression = "zstd")


## do similar for Mean ####
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

ukb_data<- arrow::read_feather(paste0(folder_path,"ukb_data.feather"))
ukb_data <- ukb_data %>%
  distinct(eid, bm_value, baseline_date,concept_name, .keep_all = TRUE)


ehr_mean_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_mean_value = mean(value_n)) %>%
  ungroup()


ukb_mean_df <- ukb_data %>%
  group_by(eid, concept_name) %>%
  summarise(bm_mean_value = mean(bm_value)) %>%
  ungroup()

matched_data_mean <- ehr_mean_df %>%
  left_join(ukb_mean_df, by = c("eid" = "eid","concept_name"="concept_name")) %>%
  filter(!is.na(bm_mean_value))

arrow::write_feather(matched_data_mean, paste0("input/02_ehr_ukbbiom_combined/combined_EHR_UKB_mean", ".feather"), compression = "zstd")

## do similar for Median ####

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

ukb_data<- arrow::read_feather(paste0(folder_path,"ukb_data.feather"))
ukb_data <- ukb_data %>%
  distinct(eid, bm_value, baseline_date,concept_name, .keep_all = TRUE)


ehr_median_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_median_value = median(value_n)) %>%
  ungroup()


ukb_median_df <- ukb_data %>%
  group_by(eid, concept_name) %>%
  summarise(bm_median_value = median(bm_value)) %>%
  ungroup()

matched_data_median <- ehr_median_df %>%
  left_join(ukb_median_df, by = c("eid" = "eid","concept_name"="concept_name")) %>%
  filter(!is.na(bm_median_value))

arrow::write_feather(matched_data_median, paste0("input/02_ehr_ukbbiom_combined/combined_EHR_UKB_median", ".feather"), compression = "zstd")

## do similar for Minimum #### 
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

ukb_data<- arrow::read_feather(paste0(folder_path,"ukb_data.feather"))
ukb_data <- ukb_data %>%
  distinct(eid, bm_value, baseline_date,concept_name, .keep_all = TRUE)


ehr_minimum_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_minimum_value = min(value_n)) %>%
  ungroup()


ukb_minimum_df <- ukb_data %>%
  group_by(eid, concept_name) %>%
  summarise(bm_minimum_value = min(bm_value)) %>%
  ungroup()

matched_data_minimum <- ehr_minimum_df %>%
  left_join(ukb_minimum_df, by = c("eid" = "eid","concept_name"="concept_name")) %>%
  filter(!is.na(bm_minimum_value))

arrow::write_feather(matched_data_minimum, paste0("input/02_ehr_ukbbiom_combined/combined_EHR_UKB_minimum", ".feather"), compression = "zstd")

## do similar for Maximum ####
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

ukb_data<- arrow::read_feather(paste0(folder_path,"ukb_data.feather"))
ukb_data <- ukb_data %>%
  distinct(eid, bm_value, baseline_date,concept_name, .keep_all = TRUE)


ehr_maximum_df <- ehr_data %>%
  group_by(eid, concept_name) %>%
  summarise(ehr_maximum_value = max(value_n)) %>%
  ungroup()


ukb_maximum_df <- ukb_data %>%
  group_by(eid, concept_name) %>%
  summarise(bm_maximum_value = max(bm_value)) %>%
  ungroup()

matched_data_maximum <- ehr_maximum_df %>%
  left_join(ukb_maximum_df, by = c("eid" = "eid","concept_name"="concept_name")) %>%
  filter(!is.na(bm_maximum_value))

arrow::write_feather(matched_data_maximum, paste0("input/02_ehr_ukbbiom_combined/combined_EHR_UKB_maximum", ".feather"), compression = "zstd")

# check ####
matched_data_max<- matched_data_maximum %>%
  filter(concept_name == "ALP")
