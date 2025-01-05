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
  "here",
  "corrplot",
  "ggcorrplot"
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

folder_path <- "input/02_ehr_ukbbiom_combined/"
combined_data<- arrow::read_feather(paste0(folder_path,"combined_EHR_UKB_min_date.feather"))


# Get unique biomarker names from ehr_data and ukb_data
biomarker_names_ehr <- unique(combined_data$concept_name_ehr)
biomarker_names_ukb <- unique(combined_data$concept_name_ukb)

# Create an empty matrix to store the correlation results
correlation_matrix <- matrix(NA, nrow = length(biomarker_names_ehr), ncol = length(biomarker_names_ukb))

# Create row and column names for the correlation matrix
rownames(correlation_matrix) <- biomarker_names_ehr
colnames(correlation_matrix) <- biomarker_names_ukb
correlation_matrix<- correlation_matrix[,order(colnames(correlation_matrix))]
correlation_matrix<- correlation_matrix[order(row.names(correlation_matrix)),]
correlation_matrix[upper.tri(correlation_matrix)] <- 0

# Loop through biomarker pairs and calculate correlations
for (i in 1:length(biomarker_names_ehr)) {
  for (j in 1:length(biomarker_names_ukb)) {
    if (is.na(correlation_matrix[i, j])) {
      cat(biomarker_names_ehr[i],"->",biomarker_names_ukb[j]," : ")
      # Filter ehr_data and ukb_long for the current pair of biomarkers
      ehr_subset <- combined_data %>% 
        filter(concept_name_ehr == biomarker_names_ehr[i]) %>%
        select(eid,date_EHR,value_n,data_provider)
      ukb_subset <- combined_data %>% 
        filter(concept_name_ukb == biomarker_names_ukb[j]) %>%
        select(eid,baseline_date,bm_value,data_provider)
      
      matched_data <- ehr_subset %>%
        inner_join(ukb_subset, by = "eid", suffix = c("_ehr", "_ukb")) %>%
        mutate(
          date_distance = date_EHR - baseline_date,
          bm_value_new = ifelse(bm_value == 0, 
                                lag(na.omit(bm_value), default = bm_value[1]), 
                                bm_value),
          value_n_new = ifelse(value_n == 0, 
                               lag(na.omit(value_n), default = value_n[1]), 
                               value_n)
        ) %>%
        group_by(eid) %>%
        mutate(min_abs_date_distance = min(abs(date_distance))) %>%
        filter(abs(date_distance) == min_abs_date_distance) %>%
        arrange(eid, min_abs_date_distance, date_distance) %>%
        slice(1) %>%
        select( -min_abs_date_distance) %>%
        mutate(bm_value = ifelse(bm_value == 0, bm_value_new, bm_value),
               value_n = ifelse(value_n == 0, value_n_new, value_n)) %>%
        select(-bm_value_new, -value_n_new)
      
      # Calculate the pearson correlation
      correlation_value <- cor(matched_data$value_n, matched_data$bm_value, method = "pearson")
      cat(correlation_value,"\n")
      # Store the correlation value in the matrix
      correlation_matrix[i, j] <- correlation_value
    }
  }
}

# make the correlation matrix symmetrical
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

correlation_matrix<- makeSymm(correlation_matrix)

# Convert the correlation matrix to a data frame
correlation_df <- as.data.frame(correlation_matrix)

write.table(correlation_matrix, "output/04_corr_plots/correlation_matrix_min_date.txt", sep="\t")

## Mean ####

# load the required data ####

folder_path <- "input/02_ehr_ukbbiom_combined/"
combined_data<- arrow::read_feather(paste0(folder_path,"combined_EHR_UKB_mean.feather"))


# Get unique biomarker names from ehr_data and ukb_long
biomarker_names_ehr <- unique(combined_data$concept_name)
biomarker_names_ukb <- unique(combined_data$concept_name)

# Create an empty matrix to store the correlation results
mean_correlation_matrix <- matrix(NA, nrow = length(biomarker_names_ehr), ncol = length(biomarker_names_ukb))
# Create row and column names for the correlation matrix
rownames(mean_correlation_matrix) <- biomarker_names_ehr
colnames(mean_correlation_matrix) <- biomarker_names_ukb

#mean_correlation_matrix<- mean_correlation_matrix[,order(colnames(mean_correlation_matrix))]
#mean_correlation_matrix<- mean_correlation_matrix[order(row.names(mean_correlation_matrix)),]
#mean_correlation_matrix[upper.tri(mean_correlation_matrix)] <- 0

# Loop through biomarker pairs and calculate correlations
for (i in 1:length(biomarker_names_ehr)) {
  for (j in 1:length(biomarker_names_ukb)) {
    if (is.na(mean_correlation_matrix[i, j])) {
      cat(biomarker_names_ehr[i],"->",biomarker_names_ukb[j],"\n")
      # Filter ehr_data and ukb_long for the current pair of biomarkers
      ehr_subset <- combined_data %>%
        filter(concept_name == biomarker_names_ehr[i]) %>%
        select(eid,concept_name,ehr_mean_value)
      
      ukb_subset <- combined_data %>%
        filter(concept_name ==  biomarker_names_ukb[j]) %>%
        select(eid,concept_name,bm_mean_value)
      
      # Join the data based on 'eid' and calculate correlations
      matched_data <- ehr_subset %>%
        left_join(ukb_subset, by = c("eid" = "eid")) %>%
        filter(!is.na(bm_mean_value))
      
      # Calculate the Spearman correlation
      correlation_value <- cor(matched_data$ehr_mean_value, matched_data$bm_mean_value, method = "pearson")
      print(correlation_value)
      # Store the correlation value in the matrix
      mean_correlation_matrix[i, j] <- correlation_value
    }
  }
}


mean_correlation_matrix <- as.matrix(mean_correlation_matrix)#
#mean_correlation_matrix<- makeSymm(mean_correlation_matrix)

write.table(mean_correlation_matrix, "output/04_corr_plots/correlation_matrix_mean.txt", sep="\t")

## Maximum ####
# load the required data ####

folder_path <- "input/02_ehr_ukbbiom_combined/"
combined_data<- arrow::read_feather(paste0(folder_path,"combined_EHR_UKB_maximum.feather"))


# Get unique biomarker names from ehr_data and ukb_long
biomarker_names_ehr <- unique(combined_data$concept_name)
biomarker_names_ukb <- unique(combined_data$concept_name)

# Create an empty matrix to store the correlation results
maximum_correlation_matrix <- matrix(NA, nrow = length(biomarker_names_ehr), ncol = length(biomarker_names_ukb))
# Create row and column names for the correlation matrix
rownames(maximum_correlation_matrix) <- biomarker_names_ehr
colnames(maximum_correlation_matrix) <- biomarker_names_ukb

#maximum_correlation_matrix<- maximum_correlation_matrix[,order(colnames(maximum_correlation_matrix))]
#maximum_correlation_matrix<- maximum_correlation_matrix[order(row.names(maximum_correlation_matrix)),]
#maximum_correlation_matrix[upper.tri(maximum_correlation_matrix)] <- 0

# Loop through biomarker pairs and calculate correlations
for (i in 1:length(biomarker_names_ehr)) {
  for (j in 1:length(biomarker_names_ukb)) {
    if (is.na(maximum_correlation_matrix[i, j])) {
      cat(biomarker_names_ehr[i],"->",biomarker_names_ukb[j],"\n")
      # Filter ehr_data and ukb_long for the current pair of biomarkers
      ehr_subset <- combined_data %>%
        filter(concept_name == biomarker_names_ehr[i]) %>%
        select(eid,concept_name,ehr_maximum_value)
      
      ukb_subset <- combined_data %>%
        filter(concept_name ==  biomarker_names_ukb[j]) %>%
        select(eid,concept_name,bm_maximum_value)
      
      # Join the data based on 'eid' and calculate correlations
      matched_data <- ehr_subset %>%
        left_join(ukb_subset, by = c("eid" = "eid")) %>%
        filter(!is.na(bm_maximum_value))
      
      # Calculate the Spearman correlation
      correlation_value <- cor(matched_data$ehr_maximum_value, matched_data$bm_maximum_value, method = "pearson")
      print(correlation_value)
      # Store the correlation value in the matrix
      maximum_correlation_matrix[i, j] <- correlation_value
    }
  }
}


maximum_correlation_matrix <- as.matrix(maximum_correlation_matrix)#

#maximum_correlation_matrix<- makeSymm(maximum_correlation_matrix)

write.table(maximum_correlation_matrix, "output/04_corr_plots/correlation_matrix_max.txt", sep="\t")

## Minimum ####

#### load the required data ####

folder_path <- "input/02_ehr_ukbbiom_combined/"
combined_data<- arrow::read_feather(paste0(folder_path,"combined_EHR_UKB_minimum.feather"))


# Get unique biomarker names from ehr_data and ukb_long
biomarker_names_ehr <- unique(combined_data$concept_name)
biomarker_names_ukb <- unique(combined_data$concept_name)

# Create an empty matrix to store the correlation results
minimum_correlation_matrix <- matrix(NA, nrow = length(biomarker_names_ehr), ncol = length(biomarker_names_ukb))
# Create row and column names for the correlation matrix
rownames(minimum_correlation_matrix) <- biomarker_names_ehr
colnames(minimum_correlation_matrix) <- biomarker_names_ukb

#minimum_correlation_matrix<- minimum_correlation_matrix[,order(colnames(minimum_correlation_matrix))]
#minimum_correlation_matrix<- minimum_correlation_matrix[order(row.names(minimum_correlation_matrix)),]
#minimum_correlation_matrix[upper.tri(minimum_correlation_matrix)] <- 0
# Loop through biomarker pairs and calculate correlations
for (i in 1:length(biomarker_names_ehr)) {
  for (j in 1:length(biomarker_names_ukb)) {
    if (is.na(minimum_correlation_matrix[i, j])) {
      cat(biomarker_names_ehr[i],"->",biomarker_names_ukb[j],"\n")
      # Filter ehr_data and ukb_long for the current pair of biomarkers
      ehr_subset <- combined_data %>%
        filter(concept_name == biomarker_names_ehr[i]) %>%
        select(eid,concept_name,ehr_minimum_value)
      
      ukb_subset <- combined_data %>%
        filter(concept_name ==  biomarker_names_ukb[j]) %>%
        select(eid,concept_name,bm_minimum_value)
      
      # Join the data based on 'eid' and calculate correlations
      matched_data <- ehr_subset %>%
        left_join(ukb_subset, by = c("eid" = "eid")) %>%
        filter(!is.na(bm_minimum_value))
      
      # Calculate the Spearman correlation
      correlation_value <- cor(matched_data$ehr_minimum_value, matched_data$bm_minimum_value, method = "pearson")
      print(correlation_value)
      # Store the correlation value in the matrix
      minimum_correlation_matrix[i, j] <- correlation_value
    }
  }
}

minimum_correlation_matrix <- as.matrix(minimum_correlation_matrix)#

#minimum_correlation_matrix<- makeSymm(minimum_correlation_matrix)

write.table(minimum_correlation_matrix, "output/04_corr_plots/correlation_matrix_minimum.txt", sep="\t")

## Median ####

### load the required data ####

folder_path <- "input/02_ehr_ukbbiom_combined/"
combined_data<- arrow::read_feather(paste0(folder_path,"combined_EHR_UKB_median.feather"))


# Get unique biomarker names from ehr_data and ukb_long
biomarker_names_ehr <- unique(combined_data$concept_name)
biomarker_names_ukb <- unique(combined_data$concept_name)

# Create an empty matrix to store the correlation results
median_correlation_matrix <- matrix(NA, nrow = length(biomarker_names_ehr), ncol = length(biomarker_names_ukb))
# Create row and column names for the correlation matrix
rownames(median_correlation_matrix) <- biomarker_names_ehr
colnames(median_correlation_matrix) <- biomarker_names_ukb

#median_correlation_matrix<- median_correlation_matrix[,order(colnames(median_correlation_matrix))]
#median_correlation_matrix<- median_correlation_matrix[order(row.names(median_correlation_matrix)),]
#median_correlation_matrix[upper.tri(median_correlation_matrix)] <- 0
# Loop through biomarker pairs and calculate correlations
for (i in 1:length(biomarker_names_ehr)) {
  for (j in 1:length(biomarker_names_ukb)) {
    if (is.na(median_correlation_matrix[i, j])) {
      cat(biomarker_names_ehr[i],"->",biomarker_names_ukb[j],"\n")
      # Filter ehr_data and ukb_long for the current pair of biomarkers
      ehr_subset <- combined_data %>%
        filter(concept_name == biomarker_names_ehr[i]) %>%
        select(eid,concept_name,ehr_median_value)
      
      ukb_subset <- combined_data %>%
        filter(concept_name ==  biomarker_names_ukb[j]) %>%
        select(eid,concept_name,bm_median_value)
      
      # Join the data based on 'eid' and calculate correlations
      matched_data <- ehr_subset %>%
        left_join(ukb_subset, by = c("eid" = "eid")) %>%
        filter(!is.na(bm_median_value))
      
      # Calculate the Spearman correlation
      correlation_value <- cor(matched_data$ehr_median_value, matched_data$bm_median_value, method = "pearson")
      print(correlation_value)
      # Store the correlation value in the matrix
      median_correlation_matrix[i, j] <- correlation_value
    }
  }
}


median_correlation_matrix <- as.matrix(median_correlation_matrix)#

#median_correlation_matrix<- makeSymm(median_correlation_matrix)

write.table(median_correlation_matrix, "output/04_corr_plots/correlation_matrix_median.txt", sep="\t")

## plots ####

###maximum ####
maximum_correlation_matrix<- read.table("output/04_corr_plots/correlation_matrix_max.txt")
# Convert the data to a correlation matrix
maximum_correlation_matrix <- as.matrix(maximum_correlation_matrix)#
#correlation_matrix<- t(correlation_matrix)
# Fill missing correlations with NA for "height" and "triglycerides"
# Add the new column to the correlation matrix
maximum_correlation_matrix <- cbind(maximum_correlation_matrix)

# 2. Find the indices of the highest correlation values in the matrix
max_cor <- which(maximum_correlation_matrix == max(maximum_correlation_matrix, na.rm = TRUE), arr.ind = TRUE)
max_corr_value <- max(maximum_correlation_matrix, na.rm = TRUE)

# Get the order of column names
col_order <- order(colnames(maximum_correlation_matrix))

# Rearrange the columns
maximum_correlation_matrix <- maximum_correlation_matrix[sort(rownames(maximum_correlation_matrix)),sort(colnames(maximum_correlation_matrix))]

# Set up the correlation plot
p0<- ggcorrplot(maximum_correlation_matrix[,31:1],hc.order = FALSE,lab = TRUE,tl.cex = 15,lab_size = 3,type = "full" ) +
  scale_fill_gradient2(limit = c(-1,1), low = "red", high =  "blue", mid = "white", midpoint = 0)

ggsave(paste("output/04_corr_plots/", "maximum_corr_plot.png", sep = ""), plot = p0, width = 20, height = 15, dpi = 500)

###minimum ####
minimum_correlation_matrix<- read.table("output/04_corr_plots/correlation_matrix_minimum.txt")
# Convert the data to a correlation matrix
minimum_correlation_matrix <- as.matrix(minimum_correlation_matrix)#
#correlation_matrix<- t(correlation_matrix)
# Fill missing correlations with NA for "height" and "triglycerides"
# Add the new column to the correlation matrix
minimum_correlation_matrix <- cbind(minimum_correlation_matrix)

# 2. Find the indices of the highest correlation values in the matrix
max_cor <- which(minimum_correlation_matrix == max(minimum_correlation_matrix, na.rm = TRUE), arr.ind = TRUE)
max_corr_value <- max(minimum_correlation_matrix, na.rm = TRUE)

# Get the order of column names
col_order <- order(colnames(minimum_correlation_matrix))

# Rearrange the columns
minimum_correlation_matrix <- minimum_correlation_matrix[sort(rownames(minimum_correlation_matrix)),sort(colnames(minimum_correlation_matrix))]

# Set up the correlation plot
p1<- ggcorrplot(minimum_correlation_matrix[,31:1],hc.order = FALSE,lab = TRUE,tl.cex = 15,lab_size = 3,type = "full" ) +
  scale_fill_gradient2(limit = c(-1,1), low = "red", high =  "blue", mid = "white", midpoint = 0)

ggsave(paste("output/04_corr_plots/", "minimum_corr_plot.png", sep = ""), plot = p1, width = 20, height = 15, dpi = 500)

###median ####
median_correlation_matrix<- read.table("output/04_corr_plots/correlation_matrix_median.txt")
# Convert the data to a correlation matrix
median_correlation_matrix <- as.matrix(median_correlation_matrix)#
#correlation_matrix<- t(correlation_matrix)
# Fill missing correlations with NA for "height" and "triglycerides"
# Add the new column to the correlation matrix
median_correlation_matrix <- cbind(median_correlation_matrix)

# 2. Find the indices of the highest correlation values in the matrix
max_cor <- which(median_correlation_matrix == max(median_correlation_matrix, na.rm = TRUE), arr.ind = TRUE)
max_corr_value <- max(median_correlation_matrix, na.rm = TRUE)

# Get the order of column names
col_order <- order(colnames(median_correlation_matrix))

# Rearrange the columns
median_correlation_matrix <- median_correlation_matrix[sort(rownames(median_correlation_matrix)),sort(colnames(median_correlation_matrix))]

# Set up the correlation plot
p2<- ggcorrplot(median_correlation_matrix[,31:1],hc.order = FALSE,lab = TRUE,tl.cex = 15,lab_size = 3,type = "full" ) +
  scale_fill_gradient2(limit = c(-1,1), low = "red", high =  "blue", mid = "white", midpoint = 0)

ggsave(paste("output/04_corr_plots/", "maximum_corr_plot.png", sep = ""), plot = p2, width = 20, height = 15, dpi = 500)

###mean ####
mean_correlation_matrix<- read.table("output/04_corr_plots/correlation_matrix_mean.txt")
# Convert the data to a correlation matrix
mean_correlation_matrix <- as.matrix(mean_correlation_matrix)#
#correlation_matrix<- t(correlation_matrix)

# Fill missing correlations with NA for "height" and "triglycerides"
# Add the new column to the correlation matrix
mean_correlation_matrix <- cbind(mean_correlation_matrix)

# 2. Find the indices of the highest correlation values in the matrix
max_cor <- which(mean_correlation_matrix == max(mean_correlation_matrix, na.rm = TRUE), arr.ind = TRUE)
max_corr_value <- max(mean_correlation_matrix, na.rm = TRUE)

# Get the order of column names
col_order <- order(colnames(mean_correlation_matrix))

# Rearrange the columns
mean_correlation_matrix <- mean_correlation_matrix[sort(rownames(mean_correlation_matrix)),sort(colnames(mean_correlation_matrix))]

# Set up the correlation plot
p3<-ggcorrplot(mean_correlation_matrix[,31:1],hc.order = FALSE,lab = TRUE,tl.cex = 7,lab_size = 2,type = "full" ) +
  scale_fill_gradient2(limit = c(-1,1), low = "red", high =  "blue", mid = "white", midpoint = 0)

ggsave(paste("output/04_corr_plots/", "mean_corr_plot.png", sep = ""), plot = p3, width = 20, height = 15, dpi = 500)

###minimum_date ####
min_date_correlation_matrix<- read.table("output/04_corr_plots/correlation_matrix_min_date.txt")
# Convert the data to a correlation matrix
min_date_correlation_matrix <- as.matrix(min_date_correlation_matrix)#
#correlation_matrix<- t(correlation_matrix)

# Fill missing correlations with NA for "height" and "triglycerides"
# Add the new column to the correlation matrix
min_date_correlation_matrix <- cbind(min_date_correlation_matrix)

# 2. Find the indices of the highest correlation values in the matrix
max_cor <- which(min_date_correlation_matrix == max(min_date_correlation_matrix, na.rm = TRUE), arr.ind = TRUE)
max_corr_value <- max(min_date_correlation_matrix, na.rm = TRUE)

# Get the order of column names
col_order <- order(colnames(min_date_correlation_matrix))

# Rearrange the columns
min_date_correlation_matrix <- min_date_correlation_matrix[sort(rownames(min_date_correlation_matrix)),sort(colnames(min_date_correlation_matrix))]

# Set up the correlation plot
p4<- ggcorrplot(min_date_correlation_matrix[,31:1],hc.order = FALSE,lab = TRUE,tl.cex = 7,lab_size = 2,type = "full" ) +
  scale_fill_gradient2(limit = c(-1,1), low = "red", high =  "blue", mid = "white", midpoint = 0)


ggsave(paste("output/04_corr_plots/", "min_date_corr_plot.png", sep = ""), plot = p4, width = 20, height = 15, dpi = 500)



