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

## compute meta-analysis across severe outcomes
# Load necessary libraries
library(metafor)
library(data.table)
Results_max <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_max.txt")

Results_max<- Results_max %>%
  mutate(id1= paste0("max"))

Results_md <- read.delim("/Aakash/12_GWAS_EHR_UKB/output/08_res_sentinels_all/Results.GWAS.all_chrs.20240317_md.txt")
Results_md<- Results_md %>%
  mutate(id1= paste0("min_date"))

data1<- Results_max %>%
  select(CHROM,GENPOS,ID,BETA,SE,id) %>%
  filter(id=="Albumin")
data2 <- Results_md %>%
  select(CHROM,GENPOS,ID,BETA,SE,id) %>%
  filter(id=="Albumin")

## compute meta-analysis
meta_data <- lapply(1:nrow(data1), function(x){
  ## compute MA across COVID-19 outcomes
  ma_result <- rma(yi=c(data1$BETA[x], data2$BETA[x]),
                   sei=c(data1$SE[x], data2$SE[x]),
                   method="FE")
  ## return results
  return(data.frame(data1[x,], data2[x,], heterogeneity=ma_result$I2, pval.hetero=ma_result$QEp))
})

## combine results
meta_data <- do.call(rbind, meta_data)

## view results
# Calculate the confidence intervals
meta_data$lower_ci <- meta_data$BETA - 1.96 * meta_data$SE
meta_data$upper_ci <- meta_data$BETA + 1.96 * meta_data$SE

# Create the forest plot
library(ggplot2)

ggplot(meta_data, aes(x = ID, y = BETA)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  coord_flip() +
  xlab("SNP ID") +
  ylab("Effect Size (BETA)") +
  ggtitle("Forest Plot of Meta-Analysis") +
  theme_minimal()

ggplot(meta_data, aes(x = BETA, y = 1 / SE)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  xlab("Effect Size (BETA)") +
  ylab("Precision (1 / SE)") +
  ggtitle("Funnel Plot of Meta-Analysis") +
  theme_minimal()

library(meta)
# Create meta-analysis object
meta_analysis <- metagen(TE = BETA, seTE = SE, data = meta_data, sm = "SMD")

# Print the meta-analysis results
print(meta_analysis)

# Heterogeneity statistics
Q <- meta_analysis$Q
I2 <- meta_analysis$I2
tau2 <- meta_analysis$tau^2

cat("Cochran's Q: ", Q, "\n")
cat("I² statistic: ", I2, "%\n")
cat("Tau²: ", tau2, "\n")
