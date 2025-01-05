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

# load the required data ####

folder_path <- "/Aakash/07_primary_care/data/00_spiros_paper/not_in_spiros"
csv_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)


EHR_data0 <- fread("/data/UK_biobank/primary_care/gp_clinical.txt",na.strings="")
EHR_data0<- as_tibble(EHR_data0)

head(EHR_data0)
## load required functions ####
####D7 remove outliers with median absolute deviation
identify_outliers <- function(data, column, threshold = 4) {
  median_val <- median(data[[column]], na.rm = TRUE)
  mad_val <- mad(data[[column]], na.rm = TRUE)
  
  lower_threshold <- median_val - threshold * mad_val
  upper_threshold <- median_val + threshold * mad_val
  
  data$outlier <- data[[column]] < lower_threshold | data[[column]] > upper_threshold
  return(data)
}

# Function to convert percentage to mmol/mol
percent_to_mmol <- function(percent) {
  return(10.929 * (percent - 2.152))
}

#TSH
csv_file<- csv_files[1]
metadat<- read_csv(csv_file)
# Extract the file name without extension
file_name <- basename(csv_file)
phenotype_name <- tools::file_path_sans_ext(file_name)

#initialization
data0<- c()
data1<- c()
# Condition 1: England Vision (1)
data0<- EHR_data0 %>%
  filter(read_2 %in% metadat$readcode | read_3 %in% metadat$readcode )


data0_units <- unique(data0$value3)
data0_units
#check if all data providers information is there
dp<- unique(data0$data_provider)
dp

#make data ready
data1 <- data0 %>%
  mutate(value1 = as.numeric(as.character(value1)),
         value2= as.numeric(as.character(value2)))

#data1$value1<- data1$`as.numeric(data0$value1)`
#data_check1<- data0 %>%
#  filter(read_2=="442W."& value3=="MEA149")

#data1$`as.numeric(data0$value1)` <- NULL

#choose according to the data_provider
data1<- data1 %>%
  filter(data1$data_provider==c(1)|data1$data_provider==c(2)|data1$data_provider==c(3)| data1$data_provider==c(4))


# Condition 1: England Vision (1)
data1$value_n <- ifelse(data1$data_provider == 1  & data1$read_2 %in% metadat$readcode, data1$value1, NA)
#data_check<- data1%>%
#  filter(data_provider==1)

data_check <- data1 %>%
  filter(data_provider == 2)

# Condition 2: Scotland (2)
data1$value_n <- ifelse(data1$data_provider == 2  & data1$read_2 %in% metadat$readcode & is.na(data1$value1), data1$value2, data1$value_n)
data1$value_n <- ifelse(data1$data_provider == 2  & data1$read_2 %in% metadat$readcode & is.na(data1$value2), data1$value1, data1$value_n)

# Condition 3: England TPP (3)
data1$value_n <- ifelse(data1$data_provider == 3  &  data1$read_3 %in% metadat$readcode, data1$value1, data1$value_n)


data_check <- data1 %>%
  filter(data_provider == 4)

# Condition 4: Wales (4)
data1$value_n <- ifelse(data1$data_provider == 4  &  data1$read_2 %in% metadat$readcode, data1$value1, data1$value_n)

#remove " " from the data values :
data1 <- data1 %>%
  filter(value_n != "")

#add the concept name
data1$concept_name<- phenotype_name 

#filter some data
data1 <- data1[!is.na(as.numeric(as.character(data1$value_n))),]
data1$value_n[data1$value_n<=0]<- NA


#remove any missings from the data
data1 <- data1 %>%
  drop_na(value_n)

# make it numerical
data1$value_n<- as.numeric(data1$value_n)

## outliers detection ####
data1 <- identify_outliers(data1, "value_n")

#remove outliers
data1<- data1 %>%
  filter(outlier == FALSE)

#making all units to uppercase
data1$units<- toupper(data1$value3)

#delete value3 column
data1$value3<- NULL

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

## verify the extracted QC'ed EHR data #### ############################################################# ####

data1$data_provider <- as.factor(data1$data_provider)
# Create a custom color palette for the data provider groups
my_colors <- c("#FF5733", "#33FF57", "#5733FF", "#FF33D1")

# Define custom labels for the legend
legend_labels <- c("England Vision", "Scotland EMIS and Vision", "England TPP", "Wales")

# Create the boxplot
p1<- ggplot(data1, aes(x = data_provider, y = value_n, fill = data_provider)) +
  geom_boxplot() +
  
  # Customize the appearance
  labs(x = "Data Provider", y = "Values") +
  ggtitle(phenotype_name) +
  scale_fill_manual(values = my_colors) +  # Use the custom color palette
  theme_minimal() +  # Apply a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for readability +
  scale_fill_discrete(name = "DataProvider", labels = legend_labels)
p1
# the figure seems good.

#data1<- data1 %>%
#  filter(value_n > 0)

#Distribution plot
#read the file
cat("Processing:", phenotype_name, "\n")
num_bins <- 1 + log2(length(data1$value_n))
binwidth <- (max(data1$value_n) - min(data1$value_n)) / num_bins

# Create histograms for each concept_name
p<- data1 %>%
  ggplot(aes(x = value_n)) +
  geom_histogram(binwidth = binwidth , color = "black", fill = "lightblue") +
  facet_wrap(~ concept_name, scales = "free") +
  labs(title = paste0("Distribution plot for ",phenotype_name),
       x = "Value",
       y = "Frequency")
p
data10 <- c()
data10<- data1 %>%
  mutate(units = ifelse(units == "" | units == " ", NA,units),
         value_n = ifelse(value_n == "" | value_n == " ", NA,value_n))
data_unique_eids<- unique(data1$eid)
#arrow::write_feather(data1, paste0("./output/25_new_extracted_phenotype_data/",phenotype_name,"_ext", ".feather"), compression = "zstd")
plot_missing(EHR_data0)
#

