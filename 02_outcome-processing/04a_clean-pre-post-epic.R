###########***********************
#### Code Description ####
# Author: Alex - updated by Nina to ensure that we don't duplicate data for ED visits that are later admitted
# Date: 2/26/24
# Updated: 2/26/24
# Goal: Obtain Patient ZCTAs
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))

# Step 1: Load data required-------------------------------------------------------------------------

wd <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/RITM0457513"
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"

# Epic Data 

epic_dx <- read_xlsx(paste0(wd,"/","RITM0457513 - Epic Diagnoses.xlsx")) %>% 
  janitor::clean_names()

epic_enc <- read_xlsx(paste0(wd,"/","RITM0457513 - Epic Encounters.xlsx"), sheet = "Sheet2") %>% 
  janitor::clean_names()

# Pre-Epic Data

pre_epic_dx <- read_xlsx(paste0(wd,"/","RITM0457513 - Pre-Epic Diagnoses.xlsx")) %>% 
  janitor::clean_names()

pre_epic_enc <- read_xlsx(paste0(wd,"/","RITM0457513 - Pre-Epic Encounters.xlsx")) %>% 
  janitor::clean_names()



# Step 2: Clean and combine admission/visit data -------------------------------------------------------------------------


# Also adjusts for inproper import of epic dates 
# Note: NAs are noted to be created by coercion, which is likely instances with missing data 

# First, this is the cleaning for the Epic data (2020-2023)

epic <- epic_dx %>% 
  left_join(epic_enc, by = "pat_enc_csn_id") %>% 
  mutate(adt_arrival_dttm = as.numeric(adt_arrival_dttm),
         ed_departure_dttm = as.numeric(ed_departure_dttm),
         inp_adm_date = as.numeric(inp_adm_date),
         hosp_admsn_time = as.numeric(hosp_admsn_time))

# Converts to datetime format 
epic$adt_arrival_dttm <- convert_date(epic$adt_arrival_dttm, type = "excel", fraction = TRUE)
epic$ed_departure_dttm <- convert_date(epic$ed_departure_dttm, type = "excel", fraction = TRUE)
epic$inp_adm_date <- convert_date(epic$inp_adm_date, type = "excel", fraction = TRUE)
epic$hosp_admsn_time <- convert_date(epic$hosp_admsn_time, type = "excel", fraction = TRUE)
epic$age <- time_length(interval(epic$birth_date, epic$adt_arrival_dttm), "years")


# Split between inpatient and outpatient 

epic_out_purely <- epic %>% 
  filter(is.na(inp_adm_date)) %>% 
  select(pat_enc_csn_id, mrn, dx_code, respiratory_dx_yn, adt_arrival_dttm, sex, ethnicity, age, patient_address) %>% 
  mutate(class = "ED") %>% 
  rename(arrival = adt_arrival_dttm) %>% 
  mutate(from = "Epic")

epic_in_purely <- epic %>% 
  filter(!is.na(inp_adm_date) & is.na(ed_departure_dttm)) %>% 
  select(pat_enc_csn_id, mrn, dx_code, respiratory_dx_yn, hosp_admsn_time, sex, ethnicity, age, patient_address) %>% 
  mutate(class = "Inpatient") %>% 
  rename(arrival = "hosp_admsn_time") %>% 
  mutate(from = "Epic")

epic_out_to_in <- epic %>% 
  filter(!is.na(inp_adm_date) & !is.na(ed_departure_dttm)) %>% 
  select(pat_enc_csn_id, mrn, dx_code, respiratory_dx_yn, adt_arrival_dttm, sex, ethnicity, age, patient_address) %>% 
  mutate(class = "ED to Inpatient") %>% 
  rename(arrival = adt_arrival_dttm) %>% #double check that using this makes sense
  mutate(from = "Epic")


# Next, this is the cleaning for the pre-Epic data (1/2016-2020)

pre_epic <- pre_epic_dx %>% 
  left_join(pre_epic_enc, by = "account") %>% 
  rename(pat_enc_csn_id = account) %>%
  mutate(age = time_length(interval(birth_date, admit_date), "years"))


pre_epic_clean <- pre_epic %>% 
 # select(pat_enc_csn_id, diagnosis_code, respiratory_dx_yn, admit_date, sex, patient_address, ip_or_ed) %>% 
  group_by(pat_enc_csn_id, diagnosis_code) %>%
  summarize(class = case_when(
    all(ip_or_ed == "Emergency") ~ "ED",
    all(ip_or_ed == "Inpatient") ~ "Inpatient",
    any(ip_or_ed == "Emergency") & any(ip_or_ed == "Inpatient") ~ "ED to Inpatient"),
    admit_date = min(admit_date),
    sex = first(sex),
    mrn = first(mrn),
    age = min(age), 
    ethnicity = first(ethnicity), 
    patient_address = first(patient_address),
    respiratory_dx_yn = first(respiratory_dx_yn)) %>%
  rename(dx_code = diagnosis_code) %>% 
  rename(adt_arrival_dttm = admit_date) %>% 
  rename(arrival = adt_arrival_dttm) %>% 
  mutate(from = "Pre-Epic")

pre_epic_out_purely<- pre_epic_clean %>% 
  filter(class == "ED")

pre_epic_in_purely<- pre_epic_clean %>% 
  filter(class == "Inpatient")

pre_epic_out_to_in <- pre_epic_clean %>% 
  filter(class == "ED to Inpatient")


# This create a dataset of all visits 

visits_pre_post_epic <- rbind(epic_out_purely, epic_in_purely, epic_out_to_in, pre_epic_out_purely, pre_epic_in_purely, pre_epic_out_to_in) %>% 
  mutate(arrival = as_date(arrival)) 


write.fst(visits_pre_post_epic, paste0(out,"/","full_encounters.fst"))



