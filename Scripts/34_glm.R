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

mdata<- fread("/Aakash/12_GWAS_EHR_UKB/input/05_indications_phe_cov/covariates.txt")

m1data<- fread("/Aakash/12_GWAS_EHR_UKB/input/05_indications_phe_cov/phenotypes_all.txt")

##################################
#### import basic information ####
##################################

## An example list of columns:
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", "f.54.0.0", "f.21001.0.0", "f.20116.0.0", "f.1558.0.0", paste0("f.22009.0.", 1:10), paste0("f.20003.0.", 1:47),"f.189.0.0")
## names to assign
cl.names        <- c("f.eid", "age", "sex", "baseline_date", "centre", "bmi", "smoking", "alcohol", paste0("pc", 1:10), paste0("self_med", 1:47),"SES")

## import data from the main release
ukb.dat         <- read_parquet("/data/UK_biobank/phenotypes/working_data/parquet_files/ukb.parquet",
                                col_select = all_of(cl.select))
## change names
names(ukb.dat)  <- cl.names
concept_names <- colnames(m1data)[3:33]
#concept_name <- "FEV1"
for(concept_name in concept_names){
  m1data_w<- m1data %>%
    select(FID,!!(concept_name))
  
  merged_data <- merge(m1data_w,ukb.dat, by.x = "FID" , by.y="f.eid" )
  
  merged_data <- merged_data %>%
    drop_na(!!concept_name)
  
  merged_data<- merged_data %>%
    select(FID,!!concept_name,age,sex,bmi,smoking,SES,alcohol)
  
  ## make biological sex binary
  merged_data$sex     <- ifelse(merged_data$sex == "Female", 1, 0)
  merged_data$smoking     <- ifelse(merged_data$smoking == "Never", 0, 1)
  merged_data$alcohol     <- ifelse(merged_data$alcohol == "Never", 0, 1)
  
  g_mod<- glm(paste0(concept_name,"~age + sex + bmi + smoking + alcohol + SES") ,merged_data,family="binomial")
  #a one-unit increase in age is associated with a 0.030093 increase in the log-odds of the outcome variable (Weight), holding other variables constant.
  
  summary(g_mod)
  
  ################################################################################
  #plots
  or_CI <- round(exp(cbind(coef(g_mod), confint(g_mod))), digits=3) %>% 
    as.data.frame()
  
  or_CI <- or_CI %>% 
    mutate(variable=rownames(or_CI)) # extract the variables from rownames
  
  or_CI <- rename(or_CI, c("AOR" = "V1",
                           "lower_bound" = "2.5 %",
                           "upper_bound" = "97.5 %"))
  
  
  # Reorder variables
  col_order <- c("variable", "AOR", "lower_bound", "upper_bound")
  or_CI <- or_CI[, col_order] #reorder variables in the data frame
  
  
  plot_logit_model <- or_CI[-1,] %>%  # Remove row number 1 (The intercept) 
    ggplot(aes(x = variable, y = AOR)) +
    geom_point(shape = 15,
               size  = 1.5,
               position = position_dodge(width = 0.1),
               color = "black") + 
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),
                  width = 0.2,
                  size  = 0.7,
                  position = position_dodge(width = 0.1),
                  color = "turquoise4") +
    theme(axis.title = element_text(face = "bold")) +
    xlab("Variables") +
    ylab(paste0("Adjusted OR with 95% CI: ", concept_name)) +
    coord_flip(ylim = c(0.5, 2.5)) + 
    geom_hline(yintercept = 1, color = "red", size = 0.5) +
    theme(axis.title = element_text(size = 17)) + 
    theme(axis.text = element_text(size = 14))
  
  
  ggsave(paste("output/26_OR_plots/",concept_name, ".png", sep = ""), plot = plot_logit_model, width = 6, height = 4, dpi = 300)
}

#plot2.0:
concept_name <- "ALT"
library(dplyr)
library(ggplot2)
library(broom)
# Prepare an empty data frame to store results
results_df <- data.frame()

for(concept_name in concept_names){
  m1data_w <- m1data %>%
    select(FID, !!sym(concept_name))
  
  merged_data <- merge(m1data_w, ukb.dat, by.x = "FID", by.y = "f.eid")
  
  merged_data <- merged_data %>%
    drop_na(!!sym(concept_name)) %>%
    select(FID, !!sym(concept_name), age, sex, bmi, smoking, SES, alcohol)
  
  # Convert variables to binary where needed
  merged_data$sex <- ifelse(merged_data$sex == "Female", 1, 0)
  merged_data$smoking <- ifelse(merged_data$smoking == "Never", 0, 1)
  merged_data$alcohol <- ifelse(merged_data$alcohol == "Never", 0, 1)
  
  # Fit the logistic regression model
  formula <- as.formula(paste0(concept_name, " ~ age + sex + bmi + smoking + alcohol + SES"))
  g_mod <- glm(formula, merged_data, family = "binomial")
  
  # Extract the model summary
  model_summary <- tidy(g_mod)
  
  # Set a lower limit for p-values to handle p-values of 0
  model_summary <- model_summary %>%
    mutate(p.value = ifelse(p.value == 0, 1e-300, p.value)) %>%
    mutate(log10_p = -log10(p.value + 1e-300)) %>%
    filter(term != "(Intercept)") %>%
    mutate(direction = ifelse(estimate > 0, "up", "down"),
           concept_name = concept_name)
  
  # Bind to the results data frame
  #results_df <- bind_rows(results_df, model_summary)
  
  # Define colors for directions
  direction_colors <- c("up" = "red", "down" = "blue")
  
  # Plotting
  plot <- ggplot(model_summary, aes(x = term, y = log10_p, color = direction)) +
    geom_point(data = subset(model_summary, direction == "up"), shape = 24, size = 3) +
    geom_point(data = subset(model_summary, direction == "down"), shape = 25, size = 3) +
    scale_color_manual(values = direction_colors, labels = c("up" = "OR > 1", "down" = "OR < 1")) +  # Updated labels
    labs(x = element_blank(), y = "-log10(p-value)",color = concept_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 350)+
    guides(color = guide_legend(override.aes = list(shape = c(25, 24))))  
  
  
  # Save the plot
  ggsave(paste("output/26_OR_plots/l_",concept_name, ".png", sep = ""), plot = plot, width = 6, height = 4, dpi = 300)
}
