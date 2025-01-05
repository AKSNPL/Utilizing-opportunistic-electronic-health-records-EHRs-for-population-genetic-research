# By Aakash Nepal
# Clear workspace and memory
rm(list = ls())
gc() # Garbage collection
options(scipen = 5) # Disable scientific notation for small numbers

# Load packages ####
library_list <- c(
  "data.table", "lubridate", "tidyr", "ggplot2", "dplyr",
  "stringr", "arrow", "readxl", "Hmisc", "tidyverse", "here"
)
suppressPackageStartupMessages(
  lapply(library_list, require, character.only = TRUE)
)

# Session information
sessioninfo::session_info()

# Set working directory
setwd("/Aakash/12_GWAS_EHR_UKB")

# Root and parallel settings ####
knitr::opts_knit$set(root.dir = here()) # Use project root for knitr

# Configure data.table threads
available_threads <- parallel::detectCores()
assigned_threads <- available_threads / 2
writeLines(paste0("Threads available: ", available_threads))
writeLines(paste0("Threads assigned to data.table: ", assigned_threads))
setDTthreads(assigned_threads)

# Load required data ####
folder_path <- "/07_primary_care/data/00_spiros_paper/ukb-biomarker-phenotypes"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

EHR_data <- fread("/data/UK_biobank/primary_care/gp_clinical.txt", na.strings = c(NA_character_, ""))
EHR_data <- as_tibble(EHR_data)

# Missingness analysis
percentage_missing <- mean(is.na(EHR_data$value3)) * 100
percentage_missing_per_column <- colMeans(is.na(EHR_data)) * 100

cat("Percentage of missing data in 'value3':", percentage_missing, "%\n")
cat("Percentage of missing data per column:\n")
print(percentage_missing_per_column)
head(EHR_data)

# Define utility functions ####

# Remove outliers using median absolute deviation
identify_outliers <- function(data, column, threshold = 4) {
  median_val <- median(data[[column]], na.rm = TRUE)
  mad_val <- mad(data[[column]], na.rm = TRUE)
  
  lower_threshold <- median_val - threshold * mad_val
  upper_threshold <- median_val + threshold * mad_val
  
  data$outlier <- data[[column]] < lower_threshold | data[[column]] > upper_threshold
  return(data)
}

# Convert percentage to mmol/mol
percent_to_mmol <- function(percent) {
  10.929 * (percent - 2.152)
}


for (csv_file in csv_files){
  
  # Read the CSV file and extract metadata
  metadat <- read.csv(csv_file)
  
 # Extract the file name without extension
  file_name <- basename(csv_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)

  if(!(phenotype_name %in% c("DBP","SBP"))){
    cat("Processing:", phenotype_name, "\n")

    # Remove trailing ".00" in readcode column while preserving the "."
    metadat$readcode <- sub("(\\.\\d*[^.])?00$", "\\1", metadat$readcode)
    
    # Remove duplicate rows based on readcode and terminology
    metadat <- metadat %>% distinct(readcode, terminology, .keep_all = TRUE)

    #initialize the list:
    ttl_calc <- c()

   # Iterate through unique readcodes and data providers
    for (readcode in unique(metadat$readcode)) {
      for (val in unique(EHR_data0$data_provider)) {
        data0 <- c()
        if(any(metadat[which(metadat$readcode==readcode),]["terminology"] == "read2") & (val== 1 | val == 2 | val == 4)){
          data0<- EHR_data0 %>%
            filter(read_2==readcode)
          
          #choose according to the data_provider
          data0<- data0 %>%
            filter(data0$data_provider==c(val))
          
          #make the value_n column as numeric
          data0$value_n<- lapply(data0$value1, function(x) as.numeric(as.character(x)))
          
          #replace NA character with NA
          data0$value_n<-gsub("NA",NA,as.character(data0$value_n))
          
          #collapse the columns if NA found in another.
          data0<- data0 %>% 
            mutate(value_n = coalesce(value_n,value2))
          
          #replace the  spaces with NA
          data0$value_n <- gsub(" ",NA,as.character(data0$value_n))
          
          # bind for saving
          ttl_calc<- rbind(ttl_calc, data0)
        }
        if(any(metadat[which(metadat$readcode==readcode),]["terminology"] == "ctv3") & ( val == 3)){
          data0<- EHR_data0 %>%
            filter(read_3==readcode)
          
          #choose according to the data_provider
          data0<- data0 %>%
            filter(data0$data_provider==c(val))
          
          #make the value_n column as numeric
          data0$value_n<- lapply(data0$value1, function(x) as.numeric(as.character(x)))
          
          #replace NA character with NA
          data0$value_n<-gsub("NA",NA,as.character(data0$value_n))
          
          #collapse the columns if NA found in another.
          data0<- data0 %>% 
            mutate(value_n = coalesce(value_n,value2))
          
          #replace the  spaces with NA
          data0$value_n<-gsub(" ",NA,as.character(data0$value_n))
          
          
          
          # bind for saving 
          ttl_calc<- rbind(ttl_calc, data0)
        }
      }
    }
    
    #uniq_read2<- unique(ttl_calc$read_2)
    #uniq_read3<- unique(ttl_calc$read_3)
    
    # Filter out rows with NA values in value_n
    ttl_calc <- ttl_calc %>% mutate(value_n = as.numeric(value_n)) %>% drop_na(value_n)
    
    # Load units mapping
    measure_ments <- data.frame(
      "Raw_value" = c("MEA000", "MEA001", "MEA016", "MEA021", "MEA026", "MEA031",
                      "MEA035", "MEA037", "MEA038", "MEA047", "MEA056", "MEA057",
                      "MEA058", "MEA061", "MEA062", "MEA071", "MEA080", "MEA082",
                      "MEA083", "MEA086", "MEA090", "MEA093", "MEA095", "MEA096",
                      "MEA097", "MEA099", "MEA104", "MEA106", "MEA110", "MEA114",
                      "MEA116", "MEA124", "MEA127", "MEA142", "MEA151", "MEA153",
                      "MEA154", "MEA156", "MEA161", "MEA169", "MEA178", "MEA183",
                      "MEA184", "MEA185", "MEA187", "MEA194", "MEA205", "MEA210",
                      "MEA215", "MEA252"),
      "Cleaned" = c("", "%", "/hour", "/mL", "1", "10^12/L", "10^6/L", "10^9/L",
                    "10^9/mL", "fL", "g/dL", "g/L", "h", "IU/L", "IU/mL", "L/min",
                    "mg", "mg/dL", "mg/L", "mg/mmol", "mL/min", "mm(Hg)", "mmol/d",
                    "mmol/L", "mmol/mol", "mol/L", "ng", "ng/mL", "nmol/L", "pg",
                    "pg/dL", "s", "U/L", "umol/L", "ratio", "10^9", "L", "mmol",
                    "1/1", "mIU/L", "ug/mL", "10^-2", "10^-3", "L/L", "nmol/g",
                    "titre", "nmol/gHb/h", "g/L", "%Hb", "mmol/g")
    )
    
    ttl_calc$units <- ttl_calc$value3

    # Merge units information and fill missing units
    ttl_calc <- merge(ttl_calc, measure_ments, by.x = "value3", by.y = "Raw_value",all.x= TRUE)
    
    #collapse the columns if NA found in another.
    ttl_calc<- ttl_calc %>% 
      mutate(Cleaned = coalesce(Cleaned,units))
    
    #sort(unique(ttl_calc$data_provider))
    
    ttl_calc <- ttl_calc %>%
      select(eid,data_provider,event_dt,value_n,Cleaned)%>%
      rename(units = Cleaned,
             date_EHR = event_dt)
    
    
    #make all units to uppercase
    ttl_calc$units<-(toupper(iconv(enc2utf8(ttl_calc$units),sub="byte")))
    #ttl_calc$units<- toupper(all(stri_enc_isutf8(as.character(ttl_calc$units))))
    #unique((testt))
    
    #uniq_umol<- ttl_calc %>%
    #  filter()
    
    # Add phenotype name as a new column
    ttl_calc$concept_name<- phenotype_name
    
    # Calculate the median values for each unit (exclude NA, empty, and space units)
    median_values <- ttl_calc %>%
      filter(!is.na(units) & units != "" & units != " ") %>%
      group_by(units) %>%
      dplyr::summarize(median_value = median(value_n, na.rm = TRUE))
    
    
    # Filter out NA, empty, or spaces from median_values$units
    non_empty_units <- median_values$units[!is.na(median_values$units) & median_values$units != "" & median_values$units != " "]
    
    # Find closest unit for NA, empty, or space units in ttl_calc
    na_indices <- which(is.na(ttl_calc$units) | ttl_calc$units == "" | ttl_calc$units == " ")
    values <- ttl_calc$value_n[na_indices]
    closest_units <- non_empty_units[apply(abs(outer(values, median_values$median_value, "-")), 1, which.min)]
    ttl_calc$units[na_indices] <- closest_units
    
    uniq_units<- as.data.frame(unique(ttl_calc$units))
    
    ## date adjustments ####
    ttl_calc <- ttl_calc %>% mutate(date_EHR = ymd(as.Date(fast_strptime(date_EHR, "%d/%m/%Y"))))
    
    # Replace NA values in 'date' with a special placeholder
    ttl_calc$date_EHR[is.na(ttl_calc$date_EHR)] <- as.Date("9999-12-31")  
    
    #flag the unreal dates
    ttl_calc <- ttl_calc %>%
      mutate(unreal_date_flag = as.Date(date_EHR) < as.Date("1902-01-01") | as.Date(date_EHR) < as.Date("1901-01-01")
             | as.Date(date_EHR) < as.Date("1904-01-01")| as.Date(date_EHR) == as.Date("9999-12-31"))
    
    ttl_calc <- ttl_calc %>%
      filter(as.Date(date_EHR) >= as.Date("1990-01-01") & as.Date(date_EHR) <= as.Date("2022-01-01"))
    
    
    ### some units unification #####
    if(phenotype_name=="HbA1c"){
      #unique(data0$units)
      
      
      # Apply the conversion function only for rows where 'unit' is "%"
      ttl_calc$value_n[ttl_calc$units == "%"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "%"])
      #ttl_calc_check<- ttl_calc %>%
      #  filter(units == "PER CENT")
      ttl_calc$units[ttl_calc$units == "%"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "%TOTAL"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "%TOTAL"])
      ttl_calc$units[ttl_calc$units == "%TOTAL"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "%TOTAL HB"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "%TOTAL HB"])
      ttl_calc$units[ttl_calc$units == "%TOTAL HB"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "% TOTAL HB"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "% TOTAL HB"])
      ttl_calc$units[ttl_calc$units == "% TOTAL HB"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "%HB"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "%HB"])
      ttl_calc$units[ttl_calc$units == "%HB"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "PER CENT"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "PER CENT"])
      ttl_calc$units[ttl_calc$units == "PER CENT"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "HBA1C"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "HBA1C"])
      ttl_calc$units[ttl_calc$units == "HBA1C"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "5"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "5"])
      ttl_calc$units[ttl_calc$units == "5"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "IU/L"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "IU/L"])
      ttl_calc$units[ttl_calc$units == "IU/L"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "UNKNOWN"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "UNKNOWN"])
      ttl_calc$units[ttl_calc$units == "UNKNOWN"] <- "MMOL/MOL"
      
      ttl_calc$value_n[ttl_calc$units == "MOL/L"] <- percent_to_mmol(ttl_calc$value_n[ttl_calc$units == "MOL/L"])
      ttl_calc$units[ttl_calc$units == "MOL/L"] <- "MMOL/MOL"
      #combined_data_HbA1c$units <- "MMOL/MOL"
      
      #unique(combined_data_HbA1c$unit)
    }
    
    #remove all lower than 0
    ttl_calc$value_n[ttl_calc$value_n<0]<- NA
    
    ## outliers detection ####
    ttl_calc <- identify_outliers(ttl_calc, "value_n")
    
    # remove all outliers
    ttl_calc<- ttl_calc %>%
      filter(outlier == FALSE)
    
    #ttl_calc1<- ttl_calc %>%
    #  filter(outlier== FALSE)
    arrow::write_feather(ttl_calc, paste0("input/01_extracted_raw_data/",phenotype_name,"_raw_ext", ".feather"), compression = "zstd")
  }
  
  ### SBP ####DBP#############################################################################################
  data0<- c()
  ttl_calc<- c()
  
  if(phenotype_name == "SBP" |phenotype_name == "DBP"){
    data1<- c()
    # Condition 1: England Vision (1)
    data0<- EHR_data0 %>%
      filter(read_2 == "246.." | read_2 == "2469." | read_2 == "246A." | read_3 == "2469."| read_3 == "246A.")
    
    #choose according to the data_provider
    data1<- data0 %>%
      filter(data0$data_provider==c(1)|data0$data_provider==c(2)|data0$data_provider==c(3)| data0$data_provider==c(4))
    
    #SBP
    if(phenotype_name == "SBP"){
      cat("SBP","\n")
      
      # Condition 1: England Vision (1)
      data1$value_n <- ifelse(data1$data_provider == 1  & data1$read_2 == "246..", data1$value1, NA)
      
      # Condition 2: Scotland (2)
      data1$value_n <- ifelse(data1$data_provider == 2  & data1$read_2 == "246..", data1$value2, data1$value_n)
      
      # Condition 3: England TPP (3)
      data1$value_n <- ifelse(data1$data_provider == 3  & data1$read_3 == "2469.", data1$value1, data1$value_n)
      
      # Condition 4: Wales (4)
      data1$value_n <- ifelse(data1$data_provider == 4  &  data1$read_2 == "246..", data1$value1, data1$value_n)
      data1$value_n <- ifelse(data1$data_provider == 4  &  data1$read_2 == "2469.", data1$value1, data1$value_n)
      
      #remove " " from the data values :
      data1 <- data1 %>%
        filter(value_n != "")
      
      #add the concept name
      data1$concept_name<- phenotype_name 
      
      data1 <- data1[!is.na(as.numeric(as.character(data1$value_n))),]
      
      #change the negative values to NA
      data1$value_n[data1$value_n<0]<- NA
      
      #remove any missings from the data
      data1 <- data1 %>%
        drop_na(value_n)
      
      #make the values numeric
      data1$value_n<- as.numeric(data1$value_n)

      ## outliers detection ####
      data1 <- identify_outliers(data1, "value_n")
      
      # remove outliers
      data1<- data1 %>%
        filter(outlier == FALSE)
      
      #making all units to uppercase
      data1$units<- toupper(data1$value3)
      
      # delete a column value3
      data1$value3<- NULL
      
      
      # Calculate the median values for each unit (excluding NA, empty, and space units)
      median_values <- data1 %>%
        filter(!is.na(units) & units != "" & units != " ") %>%
        group_by(units) %>%
        dplyr::summarize(median_value = median(value_n, na.rm = TRUE))
      
      
      # Filter out NA, empty, or spaces from median_values$units
      non_empty_units <- median_values$units[!is.na(median_values$units) & median_values$units != "" & median_values$units != " "]
      
      # Find closest unit for NA, empty, or space units in data1
      na_indices <- which(is.na(data1$units) | data1$units == "" | data1$units == " ")
      values <- data1$value_n[na_indices]
      closest_units <- non_empty_units[apply(abs(outer(values, median_values$median_value, "-")), 1, which.min)]
      data1$units[na_indices] <- closest_units
      
      uniq_units<- as.data.frame(unique(data1$units))
      
      ## date adjustments ####
      data1 <- data1 %>%
        select(eid,data_provider,event_dt,value_n,units,outlier,concept_name)%>%
        rename(date_EHR = event_dt)
      
      data1 <- data1 %>% mutate(date_EHR = ymd(as.Date(fast_strptime(date_EHR, "%d/%m/%Y"))))
      
      # Replace NA values in 'date' with a special placeholder
      data1$date_EHR[is.na(data1$date_EHR)] <- as.Date("9999-12-31")  
      
      #flag the unreal dates
      data1 <- data1 %>%
        mutate(unreal_date_flag = as.Date(date_EHR) < as.Date("1902-01-01") | as.Date(date_EHR) < as.Date("1901-01-01")
               | as.Date(date_EHR) < as.Date("1904-01-01")| as.Date(date_EHR) == as.Date("9999-12-31"))
      
      data1 <- data1 %>%
        filter(as.Date(date_EHR) >= as.Date("1990-01-01") & as.Date(date_EHR) <= as.Date("2022-01-01"))
      
      
      
      arrow::write_feather(data1, paste0("input/01_extracted_raw_data/",phenotype_name,"_raw_ext", ".feather"), compression = "zstd")
    }
    
    #DBP
    if(phenotype_name == "DBP"){
      cat("DBP","\n")
      # Condition 1: England Vision (1)
      data1$value_n <- ifelse(data1$data_provider == 1 & data1$read_2 == "246..", data1$value2, NA)
      
      # Condition 2: Scotland (2)
      data1$value_n <- ifelse(data1$data_provider == 2 & data1$read_2 == "246..", data1$value1, data1$value_n)
      
      data1$value_n <- ifelse(data1$data_provider == 2 & data1$read_2 == "246A.", data1$value1, data1$value_n)
      
      # Condition 3: England TPP (3)
      data1$value_n <- ifelse(data1$data_provider == 3 & data1$read_3 == "246A.", data1$value1, data1$value_n)
      
      # Condition 4: Wales (4)
      data1$value_n <- ifelse(data1$data_provider == 4 & data1$read_2 == "246..", data1$value2, data1$value_n)
      data1$value_n <- ifelse(data1$data_provider == 4 & data1$read_2 == "246A.", data1$value1, data1$value_n)
      
      
      
      #remove " " from the data values :
      data1 <- data1 %>%
        filter(value_n != "")
      
      #add the concept name
      data1$concept_name<- phenotype_name 
      
      #make numeric:
      data1 <- data1[!is.na(as.numeric(as.character(data1$value_n))),]
      
      
      # adjust negative values to 0
      data1$value_n[data1$value_n<0]<- NA
      
      #remove any missings from the data
      data1 <- data1 %>%
        drop_na(value_n)
      
      data1$value_n<- as.numeric(data1$value_n)
      ## outliers detection ####
      data1 <- identify_outliers(data1, "value_n")

      # Remove outliers
      data1<- data1 %>%
        filter(outlier == FALSE)
      
      #make all units to uppercase
      data1$units<- toupper(data1$value3)

      #remove the values3 column
      data1$value3<- NULL

      # Calculate the median values for each unit (excluding NA, empty, and space units)
      median_values <- data1 %>%
        filter(!is.na(units) & units != "" & units != " ") %>%
        group_by(units) %>%
        dplyr::summarize(median_value = median(value_n, na.rm = TRUE))
      
      
      # Filter out NA, empty, or spaces from median_values$units
      non_empty_units <- median_values$units[!is.na(median_values$units) & median_values$units != "" & median_values$units != " "]
      
      # Find closest unit for NA, empty, or space units in data1
      na_indices <- which(is.na(data1$units) | data1$units == "" | data1$units == " ")
      values <- data1$value_n[na_indices]
      closest_units <- non_empty_units[apply(abs(outer(values, median_values$median_value, "-")), 1, which.min)]
      data1$units[na_indices] <- closest_units
      
      uniq_units<- as.data.frame(unique(data1$units))
      
      
      ## date adjustments ####
      data1 <- data1 %>%
        select(eid,data_provider,event_dt,value_n,units,outlier, concept_name)%>%
        rename(date_EHR = event_dt)
      
      #change the format of date
      data1 <- data1 %>% mutate(date_EHR = ymd(as.Date(fast_strptime(date_EHR, "%d/%m/%Y"))))
      
      # Replace NA values in 'date' with a special placeholder
      data1$date_EHR[is.na(data1$date_EHR)] <- as.Date("9999-12-31")  
      
      #flag the unreal dates
      data1 <- data1 %>%
        mutate(unreal_date_flag = as.Date(date_EHR) < as.Date("1902-01-01") | as.Date(date_EHR) < as.Date("1901-01-01")
               | as.Date(date_EHR) < as.Date("1904-01-01")| as.Date(date_EHR) == as.Date("9999-12-31"))
      
      data1 <- data1 %>%
        filter(as.Date(date_EHR) >= as.Date("1990-01-01") & as.Date(date_EHR) <= as.Date("2022-01-01"))
      
      
      #write to the compressed feather file
      arrow::write_feather(data1, paste0("input/01_extracted_raw_data/",phenotype_name,"_raw_ext", ".feather"), compression = "zstd")
    }
    
  }
  
  
}

#ttl_calc2<- ttl_calc %>%
#  filter(units == "MG/L")
#/11_GWAS_EHR_W_H/scripts/04_req_data

#dddd<- read_feather("/11_GWAS_EHR_W_H/scripts/04_req_data/ALP.feather")

# Convert the dataprovider column to a factor
#ttl_calc1$data_provider <- as.factor(ttl_calc1$data_provider)

# Create a custom color palette for the data provider groups
#my_colors <- c("#FF5733", "#33FF57", "#5733FF", "#FF33D1")

# Define custom labels
# Define custom labels for the legend
#legend_labels <- c("England Vision", "Scotland EMIS and Vision", "England TPP", "Wales")

# Create the boxplot
#p1<- ggplot(ttl_calc1, aes(x = data_provider, y = value_n, fill = data_provider)) +
#  geom_boxplot() +
#  
#  # Customize the appearance
#  labs(x = "Data Provider", y = "Values") +
#  ggtitle(ttl_calc$concept_name) +
#  scale_fill_manual(values = my_colors) +  # Use the custom color palette
#  theme_minimal() +  # Apply a minimal theme
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for readability +
#  scale_fill_discrete(name = "DataProvider", labels = legend_labels)
#p1


#uniq_read2
#uniq_read3