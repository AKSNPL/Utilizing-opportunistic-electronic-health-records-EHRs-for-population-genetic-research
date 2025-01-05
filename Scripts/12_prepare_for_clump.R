rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
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
## minimum date #### 
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd <- fread("zcat Weight_allchrs_md.tsv.gz | head")
# Reorder columns
hd <- hd %>%
  select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)

endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")

#names:
x<- endpoint_selection$concept_name[1]

## import regional sentinels for each outcome
res.sentinel <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_md.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel     <- do.call(rbind, res.sentinel)


#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel <- res.sentinel[order(res.sentinel$id, -res.sentinel$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel <- res.sentinel %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_md.txt", sep="\t", row.names = F)

## ukb ####

setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd               <- fread("zcat Weight_allchrs_md.tsv.gz | head")
#hd$EXTRA<- NULL
# Reorder columns
hd <- hd %>%
  select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)
#bb<- names(res.region)
#hd<- names(res.region)#[1:13]
endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")
x<- endpoint_selection$concept_name


## import regional sentinels for each outcome
res.sentinel_ukb <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_ukb.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel_ukb     <- do.call(rbind, res.sentinel_ukb)

#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel_ukb <- res.sentinel_ukb[order(res.sentinel_ukb$id, -res.sentinel_ukb$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel_ukb <- res.sentinel_ukb %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel_ukb, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_ukb.txt", sep="\t", row.names = F)


#merged_df <- merge(res.sentinel, res.sentinel_ukb, by = c("ID", "id"))
#merged_df1<- merged_df%>%
#  filter(id=="Height")

#res.sentinel_ukb1<- res.sentinel_ukb %>%
#  filter(id=="Height")

#res.sentinel_md1<- res.sentinel %>%
#  filter(id=="Height")
#################################################
# List all files with the specified extension in the folder
files <- list.files(pattern = "_ukb_md.tsv.gz", full.names = TRUE)

# Initialize an empty data.table to store unique RSIDs
unique_ids <- data.table(CHROM=numeric(),GENPOS=numeric(),ID = character())
#file<- files[1]
# Loop through each file, read the RSIDs, and append to the unique_rsids data.table
for (file in files) {
  current_data <- fread(file, header = TRUE, fill=TRUE)  # Assuming the RSID column name is "RSID"
  # Extract unique RSIDs and append to the unique_rsids_df
  unique_ids <- unique(rbind(unique_ids, data.table(CHROM=current_data$CHROM,GENPOS=current_data$GENPOS,ID = current_data$ID)))
  cat("file done:",file, "\n")
}
unique_ids1 <- unique(unique_ids)

fwrite(unique_ids1,quote=FALSE, sep = "\t", row.names=FALSE,compress ="gzip", file=(paste0("../../input/04_others/","allchrs_RSIDS.tsv.gz")))

## minimum ####
rm(list = ls())
setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd               <- fread("zcat Weight_allchrs_md.tsv.gz | head")
#hd$EXTRA<- NULL
# Reorder columns
hd <- hd %>%
  dplyr::select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)
#bb<- names(res.region)
#hd<- names(res.region)#[1:13]
endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")
x<- endpoint_selection$concept_name


## import regional sentinels for each outcome
res.sentinel <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_min.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel     <- do.call(rbind, res.sentinel)

#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel <- res.sentinel[order(res.sentinel$id, -res.sentinel$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel <- res.sentinel %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_min.txt", sep="\t", row.names = F)

## mean ####
rm(list = ls())
setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd               <- fread("zcat Weight_allchrs_md.tsv.gz | head")
#hd$EXTRA<- NULL
# Reorder columns
hd <- hd %>%
  dplyr::select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)
#bb<- names(res.region)
#hd<- names(res.region)#[1:13]
endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")
x<- endpoint_selection$concept_name


## import regional sentinels for each outcome
res.sentinel <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_mean.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel     <- do.call(rbind, res.sentinel)

#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel <- res.sentinel[order(res.sentinel$id, -res.sentinel$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel <- res.sentinel %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_mean.txt", sep="\t", row.names = F)

## median ####
rm(list = ls())
setwd("/Aakash/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd               <- fread("zcat Weight_allchrs_md.tsv.gz | head")
#hd$EXTRA<- NULL
# Reorder columns
hd <- hd %>%
  dplyr::select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)
#bb<- names(res.region)
#hd<- names(res.region)#[1:13]
endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")
x<- endpoint_selection$concept_name


## import regional sentinels for each outcome
res.sentinel <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_median.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel     <- do.call(rbind, res.sentinel)

#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel <- res.sentinel[order(res.sentinel$id, -res.sentinel$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel <- res.sentinel %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_median.txt", sep="\t", row.names = F)

## maximum ####

rm(list = ls())
setwd("/12_GWAS_EHR_UKB/output/07_combined_allchrs/")
##################################################
####         import regional sentinels        ####
##################################################

## import header
hd               <- fread("zcat Weight_allchrs_md.tsv.gz | head")
#hd$EXTRA<- NULL
# Reorder columns
hd <- hd %>%
  dplyr::select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P)


hd               <- names(hd)
#bb<- names(res.region)
#hd<- names(res.region)#[1:13]
endpoint_selection <- read_excel("/Aakash/07_primary_care/data/01_relevant_data/ehr_bm_metadata.xlsx")
x<- endpoint_selection$concept_name


## import regional sentinels for each outcome
res.sentinel <- lapply(as.character(endpoint_selection$concept_name), function(x) {
  tryCatch({
    file_path <- paste0(x, "_regional_sentinels_max.txt")
    # Read the file using read.table with tab delimiter
    tmp_tab <- read.table(file_path, header = FALSE, sep = "\t")
    tmp_tab$V14 <- NULL
    # Read the file using read.table with space delimiter
    tmp_space <- read.table(file_path, header = FALSE, sep = " ")
    tmp_space$V1 <- NULL
    #tmp_tab$V14 <- NA
    tmp_tab$V15 <- tmp_space$V2
    tmp_tab$V16 <- tmp_space$V3
    # Combine the two datasets if needed
    tmp <- tmp_tab
    ## add header
    if (nrow(tmp) > 0) {
      names(tmp) <- c(hd, "region_start", "region_end")
      ## add outcome
      tmp$id <- x
      return(tmp)
    }
  }, error = function(e) {
    # Handle the error here (e.g., print an error message)
    cat("Error processing", x, ":", conditionMessage(e), "\n")
    return(NULL)  # Return NULL or any default value in case of an error
  })
})
## combine all into one file
res.sentinel     <- do.call(rbind, res.sentinel)

#-----------------------------------------------#
##-- order and look at top variants for each --##
#-----------------------------------------------#

## create vector to sort by phenotype
# Sort the data frame by 'id' and descending 'Log10P'
res.sentinel <- res.sentinel[order(res.sentinel$id, -res.sentinel$LOG10P), ]


# Use group_by and mutate to add a new column by group
res.sentinel <- res.sentinel %>%
  group_by(id) %>%
  mutate(signal_rank = row_number(dplyr::desc(LOG10P)))

write.table(res.sentinel, "../08_res_sentinels_all/Results.GWAS.all_chrs.20240317_max.txt", sep="\t", row.names = F)



