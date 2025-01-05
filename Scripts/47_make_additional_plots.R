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


# list of phenotypes:
phenotypes <- c("ALP", "ALT", "Albumin", "Basophills", "CRP", "Calcium", "Cholesterol", 
                "Creatinine", "DBP", "Eosinophills", "Glucose", "HDL", 
                "Haematocritperc", "Haemoglobinconc", "HbA1c", "Lymphocytes", "MCHbconc", 
                "MCV", "Monocytes", "Neutrophills", "Platelets", "RBC", "SBP", 
                "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides")

# load the required data ####

corr_mat<-(read.table("output/04_corr_plots/correlation_matrix_min_date.txt", sep="\t"))
corr_mat<- as.matrix(corr_mat)
diag_corr <- diag(corr_mat)


diag_min_date_df <- as.data.frame(diag_corr)

colnames(diag_min_date_df)<- ("min_date")
diag_min_date_df$id <- row.names(diag_min_date_df)


# load the required data ####

corr_mat<-(read.table("output/04_corr_plots/correlation_matrix_median.txt", sep="\t"))
corr_mat<- as.matrix(corr_mat)
diag_corr <- diag(corr_mat)

diag(corr_mat) <- 0

diag_median_df <- as.data.frame(diag_corr)

colnames(diag_median_df)<- "median"
diag_median_df$id <- row.names(diag_median_df)
# load the required data ####

corr_mat<-(read.table("output/04_corr_plots/correlation_matrix_mean.txt", sep="\t"))
corr_mat<- as.matrix(corr_mat)
diag_corr <- diag(corr_mat)

diag_mean_df <- as.data.frame(diag_corr)

colnames(diag_mean_df)<- "mean"
diag_mean_df$id <- row.names(diag_mean_df)
# load the required data ####

corr_mat<-(read.table("output/04_corr_plots/correlation_matrix_minimum.txt", sep="\t"))
corr_mat<- as.matrix(corr_mat)
diag_corr <- diag(corr_mat)

diag_min_df <- as.data.frame(diag_corr)

colnames(diag_min_df)<- "min"
diag_min_df$id <- row.names(diag_min_df)
# load the required data ####

corr_mat<-(read.table("output/04_corr_plots/correlation_matrix_max.txt", sep="\t"))
corr_mat<- as.matrix(corr_mat)
diag_corr <- diag(corr_mat)

diag_max_df <- as.data.frame(diag_corr)

colnames(diag_max_df)<- "max"
diag_max_df$id <- row.names(diag_max_df)
L <- list(diag_min_date_df,diag_median_df,diag_mean_df,diag_min_df,diag_max_df)

#merge all data frames in list
fin<- L %>% reduce(full_join, by="id")

row.names(fin)<- fin$id

fin$id <- NULL

#write.table(fin, "output/04_corr_plots/correlation_matrix_all_diag.txt", sep="\t")

# Assuming your data frame is 'fin'
fin$id <- row.names(fin)

# Reshape the data to long format and retain the 'id' column
fin_data <- fin %>%
  pivot_longer(cols = c(median, mean, min, max, min_date), names_to = "statistic", values_to = "value") %>%
  mutate(id = rep(fin$id, each = 5))  # Replicate 'id' to match the reshaped data

# Filter the data to include only values within the specified range
fin_data <- fin_data %>%
  filter(value >= 0 & value <= 0.99)

# Verify the structure of the reshaped data
str(fin_data)

# Specify the ids to highlight
highlight_ids <- c("Triglycerides", "Urea", "WBC") 

# Filter the data for the specific ids
highlight_data <- fin_data %>%
  filter(id %in% highlight_ids)

# Create a new column to indicate the original order
highlight_data <- highlight_data %>%
  mutate(order = row_number())

# Reorder with "md" values at the end, preserving the original order
highlight_data <- highlight_data %>%
  arrange(statistic == "min_date", order) %>%
  select(-order)  # Remove the helper column

# Create the violin plot with highlighted ids
ggplot(fin_data, aes(x = statistic, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, color = "darkgray", alpha = 0.6) +
  geom_line(data = highlight_data, aes(group = id, color = id), size = 1.0) +
  geom_point(data = highlight_data, aes(color = id), size = 0.5) +
  geom_text(data = highlight_data, aes(label = id), check_overlap = TRUE, vjust = -1, size = 3, fontface = "bold") +
  labs(x = element_blank(), y = "correlation: EHR-Bio") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_color_viridis_d(option = "plasma")

# Create the violin plot with enhanced styling for all IDs and connecting lines
ggplot(fin_data, aes(x = statistic, y = value, group = id)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(aes(color = id), alpha = 0.6, size = 0.8) +
  geom_jitter(width = 0.2, size = 2, color = "darkgray", alpha = 0.6) +
  geom_point(aes(color = id), size = 2, alpha = 0.8) +
  geom_text(aes(label = id), check_overlap = TRUE, vjust = -1, size = 3, fontface = "bold") +
  labs(x = "Data", y = "Value", title = "Enhanced Violin Plot with Lines for All IDs") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_color_viridis_d(option = "plasma")


ggplot(fin_data, aes(x = statistic, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_jitter(width = 0, size = 2, color = "darkgray", alpha = 0.6) +  # Set width to 0 for vertical alignment
  geom_line(data = highlight_data, aes(group = id, color = id), size = 1.0) +
  geom_point(data = highlight_data, aes(color = id), size = 1) +
  geom_text(data = highlight_data, aes(label = id), check_overlap = TRUE, vjust = -1, size = 3, fontface = "bold") +
  labs(x = "Statistic", y = "Value", title = "Enhanced Violin Plot") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_color_viridis_d(option = "plasma")

###################################################################################################
Results_max <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_max.txt")

Results_max<- Results_max %>%
  mutate(id1= paste0("max"))

Results_min <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_min.txt")
Results_min<- Results_min %>%
  mutate(id1= paste0("min"))

Results_mean <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_mean.txt")
Results_mean<- Results_mean %>%
  mutate(id1= paste0("mean"))

Results_median <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_median.txt")
Results_median<- Results_median %>%
  mutate(id1= paste0("median"))

Results_md <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_md.txt")
Results_md<- Results_md %>%
  mutate(id1= paste0("min_date"))

Results_ukb <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_ukb.txt")
Results_ukb<- Results_ukb %>%
  mutate(id1= paste0("ukb"))

ewa<- Results_ukb %>%
  filter(id== "Height")

results_all <- rbind(Results_max,Results_min,Results_mean,Results_median,Results_md,Results_ukb)

Results_oukb <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.UKB_only.202404m08.txt")
Results_oukb<- Results_oukb %>%
  mutate(id1= paste0("ukb_only"))

count_data0 <- Results_oukb %>%
  group_by(id1,id) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

# Count the number of unique ids for each signal_rank
count_data <- results_all %>%
  group_by(id1,id) %>%
  dplyr::summarize(n = n()) %>%
  ungroup()

colnames(count_data)<- c("statistic","id","value")

# Specify the IDs to highlight
highlight_ids <- c("MCV", "Monocytes", "Neutrophills", "Platelets", "RBC", "SBP" )

# Filter the data for the specific ids
highlight_data <- count_data %>%
  filter(id %in% highlight_ids)

# Create a new column to indicate the original order
highlight_data <- highlight_data %>%
  mutate(order = row_number())

# Reorder with "md" values at the end, preserving the original order
highlight_data <- highlight_data %>%
  arrange(statistic == "min_date", order) %>%
  dplyr::select(-order)  # Remove the helper column

# Create a new column to indicate the original order
count_data <- count_data %>%
  mutate(order = row_number())

count_data <- count_data %>%
  arrange(statistic == "min_date", order) %>%
  dplyr::select(-order)  # Remove the helper column

# Create the violin plot with highlighted ids
ggplot(count_data, aes(x = statistic, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, color = "darkgray", alpha = 0.6) +
  geom_line(data = highlight_data, aes(group = id, color = id), size = 1.0) +
  geom_point(data = highlight_data, aes(color = id), size = 1) +
  geom_text(data = highlight_data, aes(label = id), check_overlap = TRUE, vjust = -1, size = 3, fontface = "bold") +
  labs(x = element_blank(), y = "number of regional sentinels") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  )  +
  scale_color_viridis_d(option = "plasma")



Results_GWAS <- read_delim("output/10_result_regions/Results.GWAS_ukb_md_genes2024_ukb_sukb.txt", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE)

#Results_GWAS <- Results_GWAS %>%
#  drop_na()

Results_GWAS<- Results_GWAS %>%
  dplyr::select(MarkerName,ID,CHROM,GENPOS,region_start,region_end,id,r2,R2.group)
concept_names<- c("ALP", "ALT", "Albumin", "Basophills", "CRP", "Calcium", "Cholesterol", 
                  "Creatinine", "DBP", "Eosinophills", "Glucose", "HDL", 
                  "Haematocritperc", "Haemoglobinconc", "HbA1c", "Lymphocytes", "MCHbconc", 
                  "MCV", "Monocytes", "Neutrophills", "Platelets", "RBC", "SBP", 
                  "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides")
dats<- c("_md","_mean","_median","_min","_max")

concept_name <- "Height"
dat<- "_sukb"
for(concept_name in concept_names){
  for(dat in dats){
  # Filter data for _ukb and _max
  ukb_data <- Results_GWAS %>% filter(str_detect(id, "_ukb"))%>% filter(id==paste0(concept_name,"_ukb"))
  max_data <- Results_GWAS %>% filter(str_detect(id, dat))%>% filter(id==paste0(concept_name,dat))
  
  # Load required libraries
  library(VennDiagram)
  
  # Extract R2.group values from both datasets
  ukb_R2_groups <- unique(ukb_data$R2.group)
  max_R2_groups <- unique(max_data$R2.group)
  
  # Find the intersection of the two sets of unique values
  intersection <- intersect(ukb_R2_groups, max_R2_groups)
  
  # Create a Venn diagram
  venn.plot <- venn.diagram(
    x = list(UKB = ukb_R2_groups, Max = max_R2_groups),
    category.names = c(paste0(concept_name,"_UKBB"), paste0(concept_name,dat)),
    filename = NULL,
    output = TRUE,
    fill = c("red", "blue"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.05,
    cat.fontface = "bold",
    main = paste0("Overlap of R2.groups between UKBB and EHR",dat),
    main.cex = 2,
    margin = 0.1
  )
  # Add annotations for the totals
  png(paste0("output/27_venndiagrams/",concept_name,dat,".png"), width = 3000, height = 3000, res = 300)
  # Plot the Venn diagram
  grid.draw(venn.plot)
  
  dev.off()
  }
  
}

#/sc-projects/sc-proj-computational-medicine/people/Aakash/12_GWAS_EHR_UKB/



Results_GWAS <- read_delim("output/10_result_regions/Results.GWAS_only_UKB_genes20240408_new.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

#Res
Results_GWAS<- as.data.frame(Results_GWAS)

Results_GWAS<- Results_GWAS %>%
  dplyr::select(MarkerName,ID,CHROM,GENPOS,region_start,region_end,id,r2,R2.group)
concept_names<- c("ALP", "ALT", "Albumin", "Basophills", "CRP", "Calcium", "Cholesterol", 
                  "Creatinine", "DBP", "Eosinophills", "Glucose", "HDL", 
                  "Haematocritperc", "Haemoglobinconc", "HbA1c", "Lymphocytes", "MCHbconc", 
                  "MCV", "Monocytes", "Neutrophills", "Platelets", "RBC", "SBP", 
                  "Totalbilirubin", "Urea", "WBC", "Weight", "Height", "Triglycerides")
dats<- c("_md","_mean","_median","_min","_max")
concept_name<- "Height"
dat<- "_md"
for(concept_name in concept_names){
  for(dat in dats){
    # Filter data for _ukb and _max
    ukb_data <- Results_GWAS %>% filter(str_detect(id, "_ukb"))%>% filter(id==paste0(concept_name,"_ukb"))
    max_data <- Results_GWAS %>% filter(str_detect(id, dat))%>% filter(id==paste0(concept_name,dat))
    
    # Load required libraries
    library(VennDiagram)
    
    # Extract R2.group values from both datasets
    ukb_R2_groups <- unique(ukb_data$R2.group)
    max_R2_groups <- unique(max_data$R2.group)
    
    # Create a Venn diagram
    venn.plot <- venn.diagram(
      x = list(UKB = ukb_R2_groups, Max = max_R2_groups),
      category.names = c(paste0(concept_name,"_UKBB"), paste0(concept_name,dat)),
      filename = NULL,
      output = TRUE,
      fill = c("red", "blue"),
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.5,
      cat.pos = 0,
      cat.dist = 0.05,
      cat.fontface = "bold",
      main = paste0("Overlap of R2.groups between UKBB and EHR",dat),
      main.cex = 2,
      margin = 0.1
    )
    # Add annotations for the totals
    png(paste0("output/28_venndiagramm_wo_allUKB/",concept_name,dat,".png"), width = 3000, height = 3000, res = 300)
    # Plot the Venn diagram
    grid.draw(venn.plot)
    
    dev.off()
  }
  
}
