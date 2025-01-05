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

names<- read.table("/Aakash/12_GWAS_EHR_UKB/input/04_others/01_names.txt")

chr<- read.table("/Aakash/12_GWAS_EHR_UKB/input/04_others/chr_all.txt")

names_df <- data.frame(names=rep(names$V1,23), chr=rep(chr$V1,31))

write.table(names_df,"/Aakash/12_GWAS_EHR_UKB/input/04_others/99_names.txt", sep = "\t",col.names = FALSE,row.names = FALSE,quote =FALSE )


