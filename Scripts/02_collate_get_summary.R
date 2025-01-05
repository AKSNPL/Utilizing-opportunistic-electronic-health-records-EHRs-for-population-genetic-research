# By Aakash Nepal

# Clear workspace and memory
rm(list = ls())
gc() # Garbage collection
options(scipen = 5) # Disable scientific notation for small numbers

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
folder_path <- "input/01_extracted_raw_data/"
feather_files <- list.files(path = folder_path, pattern = "\\.feather$", full.names = TRUE)
#feather_file<- feather_files[4]

## Boxplots Data providers vs Values for the phenotypes ####
for(feather_file in feather_files){
  #read the file
  ttl_calc <- arrow::read_feather(feather_file)
  # Extract the file name without extension
  file_name <- basename(feather_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)
  cat("Processing:", phenotype_name, "\n")
  #phenotype_name <- unique(ttl_calc$concept_name)
  # Convert the dataprovider column to a factor
  ttl_calc$data_provider <- as.factor(ttl_calc$data_provider)
  # Create a custom color palette for the data provider groups
  my_colors <- c("#FF5733", "#33FF57", "#5733FF", "#FF33D1")
  
  # Define custom labels for the legend
  legend_labels <- c("England Vision", "Scotland EMIS and Vision", "England TPP", "Wales")
  
  
  # Create the boxplot
  p1<- ggplot(ttl_calc, aes(x = data_provider, y = value_n, fill = data_provider)) +
    geom_boxplot() +
    
    # Customize the appearance
    labs(x = "Data Provider", y = "Values") +
    ggtitle(phenotype_name) +
    scale_fill_manual(values = my_colors) +  # Use the custom color palette
    theme_minimal() +  # Apply a minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for readability +
    scale_fill_discrete(name = "DataProvider", labels = legend_labels)
   #p1
  
  # Save the plot with the name of the current concept_name
  ggsave(paste("output/01_boxplots_DP_vs_Val/",phenotype_name, ".png", sep = ""), plot = p1, width = 6, height = 4, dpi = 300)
  
}

## Histograms Counts per Year for the phenotypes ####
for(feather_file in feather_files){
  #read the file
  ttl_calc <- arrow::read_feather(feather_file)
  # Extract the file name without extension
  file_name <- basename(feather_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)
  cat("Processing:", phenotype_name, "\n")
  
  ttl_calc$year <- year(ttl_calc$date_EHR)
  
  count_data1 <- ttl_calc %>%
    group_by(year) %>%
    dplyr::summarise(count = n_distinct(eid))%>%
    drop_na()
  
  
  # Create a histogram using ggplot2 with nicer aesthetics
  p00<- ggplot(count_data1, aes(x = as.factor(year), y = count)) +
    geom_bar(stat = "identity", color = "white") +
    geom_text(aes(label = count), vjust = -0.5, size = 3) +  # Add count numbers at the top of each bar
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Year", y = "Count", title = paste0("Counts per Year for ",phenotype_name)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),  # Do not rotate x-axis labels and center them
      legend.position = "right",   # Move the legend to the right
      legend.title = element_blank(),  # Remove the legend title
      plot.title = element_text(size = 16, hjust = 0.5)  # Customize plot title size and alignment
    )
  # Save the plot with the name of the current concept_name
  ggsave(paste("output/02_histograms_per_year/",phenotype_name, ".png", sep = ""), plot = p00, width = 10, height = 4, dpi = 300)
  #p00
}


## summary stats#### #####################
summary_total<-c()
## Histograms Counts per Year for the phenotypes ####
for(feather_file in feather_files){
  summary_stats<- c()
  #read the file
  ttl_calc <- arrow::read_feather(feather_file)
  # Extract the file name without extension
  file_name <- basename(feather_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)
  cat("Processing:", phenotype_name, "\n")
  summary_stats <- ttl_calc %>%
    mutate(date = as.Date(date_EHR)) %>%
    mutate(value = as.numeric(value_n)) %>%
    filter(date >= as.Date("1904-01-01") & date <= as.Date("2023-12-31")) %>%
    group_by(concept_name) %>%
    dplyr::summarize(
      mean_value = round(mean(value, na.rm = TRUE), 2),
      median_value = round(median(value, na.rm = TRUE), 2),
      lower_bound = round(min(value, na.rm = TRUE), 2),
      upper_bound = round(max(value, na.rm = TRUE), 2),
      start_date = min(date),
      end_date = max(date),
      duration_years = year(max(date)) - year(min(date)) + 1,
      last_n_eids = n_distinct(eid),
      missing_or_empty_unit_proportion = sum(units %in% c(NA,"", " ")) / n()
    )
  
  # bind for saving
  summary_total<- rbind(summary_total, summary_stats)
  
}

## histograms DISTRIBUTION PLOTS###########

for(feather_file in feather_files){
  summary_stats<- c()
  #read the file
  ttl_calc <- arrow::read_feather(feather_file)
  # Extract the file name without extension
  file_name <- basename(feather_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)
  cat("Processing:", phenotype_name, "\n")
  num_bins <- 1 + log2(length(ttl_calc$value_n))
  binwidth <- (max(ttl_calc$value_n) - min(ttl_calc$value_n)) / num_bins
  
  # Create histograms for each concept_name
  p<- ttl_calc %>%
    ggplot(aes(x = value_n)) +
    geom_histogram(binwidth = binwidth , color = "black", fill = "lightblue") +
    facet_wrap(~ concept_name, scales = "free") +
    labs(title = paste0("Distribution plot for ",phenotype_name),
         x = "Value",
         y = "Frequency")
  #p
  # Save the plot with the name of the current concept_name
  ggsave(paste("output/03_distribution_plots/",phenotype_name, ".png", sep = ""), plot = p, width = 6, height = 4, dpi = 300)
  
}


