# By Aakash Nepal
# Clear the environment and perform garbage collection
rm(list = ls()) 
gc()

# Avoid scientific notation for numbers
options(scipen = 5)

# Load required packages
library(sessioninfo)
sessioninfo::session_info()

# Set the working directory
setwd("/Aakash/12_GWAS_EHR_UKB")

# Load necessary libraries
required_packages <- c(
  "data.table", "lubridate", "tidyr", "ggplot2", "dplyr", "stringr", 
  "arrow", "readxl", "Hmisc", "tidyverse", "here"
)

for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Configure root directory for `knitr` and plot settings
knitr::opts_knit$set(root.dir = here())

# Parallel settings for data.table
threads <- parallel::detectCores()
cat("Threads available:", threads, "\n")
cat("Threads assigned to data.table:", threads / 2, "\n")
setDTthreads(threads / 2)

# ====================================
# Load and process basic information
# ====================================

# Define columns to import and rename
columns_to_select <- c(
  "f.eid", "f.21022.0.0", "f.31.0.0", "f.53.0.0", 
  "f.54.0.0", "f.21001.0.0", "f.20116.0.0", "f.1558.0.0", 
  paste0("f.22009.0.", 1:10), paste0("f.20003.0.", 1:47)
)
column_names <- c(
  "f.eid", "age", "sex", "baseline_date", "centre", 
  "bmi", "smoking", "alcohol", paste0("pc", 1:10), paste0("self_med", 1:47)
)

# Load UK Biobank data
ukb_data <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = columns_to_select
)
colnames(ukb_data) <- column_names
ukb_data <- as.data.table(ukb_data)

# ====================================
# Process medication mapping
# ====================================

# Import medication mapping
med_mapping <- read_excel(
  "/data/UK_biobank/medications/UK_medication_to_ATC_code_mapping.xlsx", skip = 1
) %>% 
  select(Category, Coding.a, Medication.ATC.code, Drug.name) %>%
  as.data.frame()

# Expand ATC codes for each medication
med_mapping <- do.call(rbind, lapply(1:nrow(med_mapping), function(x) {
  atc_codes <- strsplit(med_mapping$Medication.ATC.code[x], " \\|")[[1]]
  data.frame(
    category = med_mapping$Category[x],
    ukbb_coding = med_mapping$Coding.a[x],
    atc_code = atc_codes,
    drug_name = med_mapping$Drug.name[x]
  )
}))

# Add UK Biobank coding information
ukb_coding <- read.delim("/data/UK_biobank/medications/coding4.tsv")
med_mapping <- merge(med_mapping, ukb_coding, by.x = "ukbb_coding", by.y = "coding")

# ====================================
# Create long-format dataset
# ====================================

# Reshape UK Biobank data for medication records
med_data <- melt(
  ukb_data[, c("f.eid", paste0("self_med", 1:47))],
  id.vars = "f.eid"
)[!is.na(value)]

# Map ATC codes to medications
med_data <- merge(
  med_data, unique(med_mapping[, c("ukbb_coding", "atc_code")]), 
  by.x = "value", by.y = "ukbb_coding", allow.cartesian = TRUE
)

# Reduce to one ATC code per participant
med_data <- med_data[order(f.eid, atc_code, value)]
med_data[, ind := 1:.N, by = c("f.eid", "atc_code")]
med_data <- med_data[ind == 1]

# Convert to wide format
med_data <- dcast(med_data, f.eid ~ atc_code)

# Merge back with UKB data
ukb_self <- merge(ukb_data, med_data, all.x = TRUE, by = "f.eid")

# ====================================
# Process and recode medication data
# ====================================

# Filter medication mapping for relevant ATC codes
med_mapping <- med_mapping[med_mapping$atc_code %in% colnames(ukb_self)]

# Recode ATC columns in the dataset (1 = reported, 0 = not reported)
ukb_self <- as.data.frame(ukb_self)
ukb_self[, med_mapping$atc_code] <- lapply(
  ukb_self[, med_mapping$atc_code], 
  function(x) ifelse(!is.na(x), 1, 0)
)
ukb_self <- as.data.table(ukb_self)

# ====================================
# Count reporting statistics
# ====================================

# Count participants for each ATC code
atc_summary <- data.table(
  atc_code = unique(med_mapping$atc_code),
  count_self = sapply(unique(med_mapping$atc_code), function(x) {
    nrow(ukb_self[get(x) == 1])
  }),
  n_men = sapply(unique(med_mapping$atc_code), function(x) {
    nrow(ukb_self[get(x) == 1 & sex == "Male"])
  }),
  n_women = sapply(unique(med_mapping$atc_code), function(x) {
    nrow(ukb_self[get(x) == 1 & sex == "Female"])
  })
)

# Add sex-specific prescription information
atc_summary[, sex := ifelse(
  n_men > 0 & n_women > 0, "both", 
  ifelse(n_men > 0 & n_women == 0, "men", "women")
)]



##################################
####   Augmented Phenotypes   ####
##################################

# Load the main UK Biobank data dictionary
lab_main <- fread("/data/UK_biobank/phenotypes/working_data/Data.dictionary.UKBB.main.dataset.45268.txt")

# -------------------------------
# Handle Technical Variables
# -------------------------------

# Variables: time of blood draw, fasting time, number of draws, month of visit, OLink plate ID, OLink well ID, and derived OLink PCs

# Load the list of technical variables
tech_var <- fread("/data/variables_techincal_20230414.txt")

# Generate short column names for ease of use
tech_var[, column_name := paste0("f.", id, ".0.0")] # Format field IDs as column names
tech_var[, short_name := gsub(" |\\(|\\)|,|-", "_", label)] # Create simplified, short names

# Check which variables are already available in the main dataset
tech_var[, released := column_name %in% lab_main$id.ukbb]

# Import relevant data from the main UK Biobank release
ukb_tech <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", tech_var$column_name)
)

# Rename columns to use short names
colnames(ukb_tech) <- c("f.eid", tech_var$short_name)
ukb_tech <- as.data.table(ukb_tech)

# --------------------------------
# Process Month Variable
# --------------------------------

# Extract the month from the assessment date
ukb_tech[, month := month(Date_of_attending_assessment_centre)]

# Create polynomial transformations of the month variable for non-linear modeling
ukb_tech[, month_2 := month^2] # Square of the month
ukb_tech[, month_3 := month^3] # Cube of the month

# Add the new variables to the technical variable list
tech_var <- rbind(
  tech_var,
  data.table(
    id = NA,
    label = c("Month", "Month^2", "Month^3"),
    category = "Technical",
    column_name = NA,
    short_name = c("month", "month_2", "month_3"),
    released = TRUE
  )
)

# --------------------------------
# Notes
# --------------------------------

# 1. The `month` variable represents the month of the participant's assessment visit, allowing for seasonality analysis.
# 2. Polynomial transformations (`month_2`, `month_3`) are included for potential use in non-linear models.
# 3. If additional variables related to age or OLink-specific information already exist in other datasets, they can be merged at later stages.

# Clean up the environment to free memory
gc()



# ------------------------ #
#     Body Composition     #
# ------------------------ #

# Import the list of body composition variables
body_var <- fread("/data/variables_body_composition_20230414.txt")

# Generate column names and short names for ease of use
body_var[, column_name := paste0("f.", id, ".0.0")]
body_var[, short_name := gsub(" |\\(|\\)|,|-", "_", label)]

# Load body composition data from the main release
ukb_body <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", body_var$column_name)
)

# Rename columns to short names
colnames(ukb_body) <- c("f.eid", body_var$short_name)

# Calculate Waist-to-Hip Ratio (WHR)
ukb_body <- as.data.table(ukb_body)
ukb_body[, whr := Waist_circumference / Hip_circumference]

# Add WHR to the label file
body_var <- rbind(
  body_var,
  data.table(
    id = NA,
    label = "Waist-to-Hip ratio",
    category = "Body composition",
    column_name = NA,
    short_name = "whr"
  )
)

# Calculate missing value percentages for each variable
body_var[, miss_per := sapply(short_name, function(x) {
  nrow(ukb_body[is.na(eval(as.name(x)))])
}) / nrow(ukb_body) * 100]

# ------------------------ #
#    Socioeconomic Data    #
# ------------------------ #

# Import the list of socioeconomic variables
soec_var <- fread("/data/variables_socioeconomic_20230414.txt")

# Generate column names and short names
soec_var[, column_name := paste0("f.", id, ".0.0")]
soec_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]

# Identify variables already present in the dataset
soec_var[, released := column_name %in% lab_main$id.ukbb]

# Retain only variables available in the main dataset
soec_var <- soec_var[released == TRUE]

# Load socioeconomic data from the main release
ukb_soec <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", soec_var$column_name)
)

# Rename columns to short names
colnames(ukb_soec) <- c("f.eid", soec_var$short_name)

# Analyze socioeconomic variables for data quality
table(ukb_soec$Leisure_social_activities, useNA = "always")
table(ukb_soec$Age_completed_full_time_education, useNA = "always")  # Many NAs
table(ukb_soec$Qualifications, useNA = "always")
table(ukb_soec$Average_total_household_income, useNA = "always")
table(ukb_soec$Duration_of_moderate_activity, useNA = "always")  # High missing or invalid values

# ------------------------ #
#         Biomarkers       #
# ------------------------ #

# Load QC'd biomarker levels
ukb_bio <- read_parquet("/data/UK_biobank/biomarker/working_data/Cleaned.UKB.Biomarker.20230414.parquet")

# Select relevant biomarker variables
ukb_bio <- ukb_bio %>%
  select(
    f.eid, alp, alt, crea_s, tbil, total_tg, albumin, crp, ca,
    alt_date, hdl_chol_cleaned, chol_cleaned, glucose_cleaned,
    hba1c_cleaned, urea_cleaned
  )

# Clean memory
gc(reset = TRUE)



# ------------------------ #
#          Diseases        #
# ------------------------ #

# Import Phecode assignment and convert to wide format
ukb_phe <- fread(
  "/data/UK_biobank/phecodes_light/final/UKBB.phecodes.collated.long.20210922.txt",
  sep = "\t", header = TRUE
)
ukb_phe <- dcast(ukb_phe, f.eid ~ phecode, value.var = c("date", "resource"))

# Retain only date columns
ukb_phe <- ukb_phe[, c("f.eid", grep("date", names(ukb_phe), value = TRUE)), with = FALSE]

# Import Phecode labels
lab_phe <- fread(
  "/data/UK_biobank/phecodes/data/Comparison.phecodes.cases.UKBB.UCL.Charite.20210811.txt",
  sep = "\t", header = TRUE
)
lab_phe[, id := paste0("date_", phecode)]

# Filter labels to those present in the dataset
lab_phe <- lab_phe[id %in% names(ukb_phe)]

# Add sex-specific information
lab_tmp <- fread("/data/UK_biobank/phecodes/public/phecode_definitions1.2.csv")
lab_tmp[, id := paste0("date_", phecode)]
lab_phe <- merge(lab_phe, lab_tmp[, .(id, sex)], by = "id", all.x = TRUE)
lab_phe[, sex := ifelse(is.na(sex) | sex == "", "Both", sex)]

# Add baseline date and convert to data frame
ukb_phe <- merge(ukb_phe, ukb_dat[, .(f.eid, baseline_date)], all = TRUE)
ukb_phe$baseline_date <- as.IDate(ukb_phe$baseline_date)

# Create binary case status for each Phecode
for (j in lab_phe$id) {
  bin_col <- gsub("date", "bin", j)
  ukb_phe[, (bin_col) := ifelse(
    is.na(.SD[[j]]) | .SD[[j]] > baseline_date, 0, 1
  ), .SDcols = j]
}

# Merge sex information for case counts
ukb_phe <- merge(ukb_phe, ukb_dat[, .(f.eid, sex)], by = "f.eid", all.x = TRUE)

# Compute case counts
lab_phe[, cases := apply(.SD, 1, function(x) {
  if (x[2] == "Both") {
    sum(ukb_phe[[gsub("date", "bin", x[1])]] == 1, na.rm = TRUE)
  } else {
    sum(ukb_phe[sex == x[2], gsub("date", "bin", x[1])] == 1, na.rm = TRUE)
  }
}), .SDcols = c("id", "sex")]

# Retain Phecodes with at least 100 cases
lab_phe <- lab_phe[cases >= 100]

# Retain only necessary columns in the dataset
ukb_phe <- ukb_phe[, .(
  f.eid,
  baseline_date,
  mget(paste0("bin_", lab_phe$phecode))
)]
gc(reset = TRUE)


# ------------------------ #
#           Diet           #
# ------------------------ #

# Import diet variables
diet_var <- fread("/data/variables_diet_20230428.txt")

# Create column and short names
diet_var[, column_name := paste0("f.", id, ".0.0")]
diet_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]

# Filter released variables
diet_var <- diet_var[column_name %in% lab_main$id.ukbb]

# Load diet data
ukb_diet <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", diet_var$column_name)
)

# Rename columns to short names
names(ukb_diet) <- c("f.eid", diet_var$short_name)

# Summary of diet variables
summary(ukb_diet)


# ------------------------ #
#     Blood Cell Counts    #
# ------------------------ #

# Import blood cell count variables
blood_var <- fread("/data/variables_blood_cell_20230428.txt")

# Create column and short names
blood_var[, column_name := paste0("f.", id, ".0.0")]
blood_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]

# Filter released variables
blood_var <- blood_var[column_name %in% lab_main$id.ukbb]

# Load blood cell count data
ukb_blood <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", blood_var$column_name)
)

# Rename columns to short names
names(ukb_blood) <- c("f.eid", blood_var$short_name)

# Prune outliers
ukb_blood[, (blood_var$short_name) := lapply(.SD, function(x) {
  meds <- median(x, na.rm = TRUE)
  mads <- mad(x, na.rm = TRUE)
  ifelse(x < (meds - 5 * mads) | x > (meds + 5 * mads), NA, x)
}), .SDcols = blood_var$short_name]

# Omit nucleated red blood cells from analysis
blood_var <- blood_var[!short_name %in% c("Nucleated_red_blood_cell_count", "Nucleated_red_blood_cell_percentage")]

# ------------------------ #
#         Pulmonary        #
# ------------------------ #

# Import pulmonary variables
pulm_var <- fread("/data/variables_pulmonary_20230503.txt")

# Create column and short names
pulm_var[, column_name := paste0("f.", id, ".0.0")]
pulm_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", label)]

# Filter released variables
pulm_var <- pulm_var[column_name %in% lab_main$id.ukbb]

# Load pulmonary data
ukb_pulm <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", pulm_var$column_name)
)

# Rename columns to short names
names(ukb_pulm) <- c("f.eid", pulm_var$short_name)

# Summary of pulmonary variables
summary(ukb_pulm)


# ------------------------ #
#     Cardiovascular       #
# ------------------------ #

# Import cardiovascular variables to be selected
card_var <- fread("/data/variables_cardio_20230428.txt")

# Create short names for variables by removing unwanted characters
card_var[, column_name := paste0("f.", id, ".0.0")]
card_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]

# Check if these variables are already in the dataset
card_var[, released := column_name %in% lab_main$id.ukbb]

# Filter the variables that are released and present in the dataset
card_var <- card_var[released == TRUE]

# Import cardiovascular data from the main release (only required columns)
ukb_card <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", card_var$column_name)
)

# Rename columns using the short names
names(ukb_card) <- c("f.eid", card_var$short_name)

# Show summary of the cardiovascular dataset
summary(ukb_card)


# ------------------------ #
#       Bone Health        #
# ------------------------ #

# Import bone-related variables to be selected
bone_var <- fread("/data/variables_bone_20230503.txt")

# Create short names for variables by cleaning up labels
bone_var[, column_name := paste0("f.", id, ".0.0")]
bone_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]

# Check if these variables are already in the dataset
bone_var[, released := column_name %in% lab_main$id.ukbb]

# Filter only the released variables that are present in the dataset
bone_var <- bone_var[released == TRUE]

# Import bone data from the main release (only required columns)
ukb_bone <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", bone_var$column_name)
)

# Rename columns using the short names
names(ukb_bone) <- c("f.eid", bone_var$short_name)

# Show summary of the bone health dataset
summary(ukb_bone)

# ------------------------ #
#       General Health     #
# ------------------------ #

# Import general health-related variables to be selected
health_var <- fread("/data/variables_health_20230428.txt")

# Create short names for variables by cleaning up labels
health_var[, column_name := paste0("f.", id, ".0.0")]
health_var[, short_name := gsub(" |\\(|\\)|,|-|\\/", "_", gsub(", automated reading", "", label))]

# Check if these variables are already in the dataset
health_var[, released := column_name %in% lab_main$id.ukbb]

# Filter only the released variables that are present in the dataset
health_var <- health_var[released == TRUE]

# Import health data from the main release (only required columns)
ukb_health <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", health_var$column_name)
)

# Rename columns using the short names
names(ukb_health) <- c("f.eid", health_var$short_name)

# Show summary of the general health dataset
summary(ukb_health)


# ------------------------ #
#        Pollution         #
# ------------------------ #

# Import pollution-related variables to be selected
poll_var <- fread("/data/variables_pollution_20230428.txt")

# Create short names for variables by cleaning up labels
poll_var[, column_name := paste0("f.", id, ".0.0")]
poll_var[, short_name := gsub(" |\\(|\\)|,|-|\\/|;", "_", gsub(", automated reading", "", label))]

# Check if these variables are already in the dataset
poll_var[, released := column_name %in% lab_main$id.ukbb]

# Filter only the released variables that are present in the dataset
poll_var <- poll_var[released == TRUE]

# Import pollution data from the main release (only required columns)
ukb_poll <- read_parquet(
  "/data/UK_biobank/phenotypes/working_data/parquet_files/ukb45268.parquet",
  col_select = c("f.eid", poll_var$column_name)
)

# Rename columns using the short names
names(ukb_poll) <- c("f.eid", poll_var$short_name)

# Show summary of the pollution dataset
summary(ukb_poll)


##################################
####    combined data set     ####
##################################

# ------------------------------ #
#   Combine Phenotypic Data      #
# ------------------------------ #

# Combine various datasets into one using a full join on the "f.eid" column
ukb_comb <- list(ukb.bio, ukb.blood, ukb.body, ukb.bone, ukb.card, 
                 ukb.dat[, c("f.eid", "age", "sex", "centre", "smoking", "alcohol", paste0("pc", 1:10)), with=FALSE],
                 ukb.diet, ukb.health, ukb.poll, ukb.pulm, ukb.soec, ukb.tech) %>% 
  reduce(full_join, by="f.eid")

# Add medication information by merging with self-reported data
ukb_comb <- merge(ukb_comb, ukb.self[, c("f.eid", atc.self[count.self >= 50]$atc_code), with=FALSE], by="f.eid", all.x=TRUE)

# Add disease data (phecodes) by merging with phenotypic data
ukb_comb <- merge(ukb_comb, ukb.phe, by="f.eid", all.x=TRUE)

# Convert the combined data to a data.table format for efficiency
ukb_comb <- as.data.table(ukb_comb)

# Clean up by removing intermediate datasets that are no longer needed
rm(list=c("ukb.bio", "ukb.blood", "ukb.body", "ukb.bone", "ukb.card", "ukb.diet", "ukb.health", "ukb.poll", "ukb.pulm", "ukb.soec", "ukb.tech", "ukb.phe", "ukb.self"))
gc(reset = TRUE)

# ------------------------------ #
#   Select and Rename Columns     #
# ------------------------------ #

# Select relevant columns for further analysis
ukb_data_1 <- ukb_comb %>%
  select(f.eid, alp, alt, crea_s, tbil, total_tg, albumin, crp, ca, baseline_date, 
         hdl_chol_cleaned, chol_cleaned, Eosinophill_count, Diastolic_blood_pressure, 
         Systolic_blood_pressure, Haematocrit_percentage, Haemoglobin_concentration, 
         glucose_cleaned, Basophill_count, Lymphocyte_count, hba1c_cleaned, Monocyte_count,
         Neutrophill_count, Platelet_count, urea_cleaned, White_blood_cell__leukocyte__count, 
         Red_blood_cell__erythrocyte__count, Forced_vital_capacity__FVC_, Weight, 
         Mean_corpuscular_volume, Forced_expiratory_volume_in_1_second__FEV1_, 
         Mean_corpuscular_haemoglobin, Standing_height, sex)

# Rename variables for clarity and consistency
ukb_data_1 <- ukb_data_1 %>%
  rename(
    eid = f.eid,
    Albumin = albumin,
    ALP = alp,
    ALT = alt,
    Basophills = Basophill_count,
    Calcium = ca,
    Cholesterol = chol_cleaned,
    Creatinine = crea_s,
    CRP = crp,
    DBP = Diastolic_blood_pressure,
    Eosinophills = Eosinophill_count,
    FEV1 = Forced_expiratory_volume_in_1_second__FEV1_,
    FVC = Forced_vital_capacity__FVC_,
    Glucose = glucose_cleaned,
    Haematocritperc = Haematocrit_percentage,
    Haemoglobinconc = Haemoglobin_concentration,
    HbA1c = hba1c_cleaned,
    HDL = hdl_chol_cleaned,
    Lymphocytes = Lymphocyte_count,
    MCHbconc = Mean_corpuscular_haemoglobin,
    MCV = Mean_corpuscular_volume,
    Monocytes = Monocyte_count,
    Neutrophills = Neutrophill_count,
    Platelets = Platelet_count,
    RBC = Red_blood_cell__erythrocyte__count,
    SBP = Systolic_blood_pressure,
    Totalbilirubin = tbil,
    Urea = urea_cleaned,
    WBC = White_blood_cell__leukocyte__count,
    Weight = Weight,
    Height = Standing_height,
    Triglycerides = total_tg
  )


# Convert the dataframe to long format
ukb_data <- gather(ukb_data_1, concept_name, bm_value, -eid,-baseline_date,-sex, na.rm = TRUE)

arrow::write_feather(ukb_data, paste0("input/02_ehr_ukbbiom_combined/ukb_data", ".feather"), compression = "zstd")

#####EHR ##########################################
folder_path <- "input/01_extracted_raw_data/"
feather_files <- list.files(path = folder_path, pattern = "\\.feather$", full.names = TRUE)
#feather_file<- feather_files[4]
ehr_data<- c()

## EHR files ####
for(feather_file in feather_files){
  #read the file
  ttl_calc <- arrow::read_feather(feather_file)
  # Extract the file name without extension
  file_name <- basename(feather_file)
  phenotype_name <- tools::file_path_sans_ext(file_name)
  ehr_data<-rbind(ehr_data,ttl_calc)
  
}

#save to the compressed feather files
arrow::write_feather(ehr_data, paste0("input/02_ehr_ukbbiom_combined/ehr_data", ".feather"), compression = "zstd")

