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


visits <- read.fst(paste0(out,"/","full_encounters.fst"))

# ZIP Crosswalk 

zip_xwalk <- read_xlsx(here("data", "census", "ZIPCodetoZCTACrosswalk2019.xlsx")) %>% 
  janitor::clean_names() %>% 
  select(zip_code, zcta) 

# NYC ZCTAs

nyc_zctas <- read_csv(here("data", "nyc", "modzcta.csv")) %>% 
  janitor::clean_names() %>%
  separate_rows(zcta, sep = ",") %>%
  select(zcta)

# ZCTAs 
# Note, these are 2010 ZCTAs

zcta_shapefile <- zctas(year = 2010) %>% 
  janitor::clean_names() %>% 
  dplyr::select(zcta5ce10, geometry) %>%
  st_transform(., crs = 3748) %>% 
  rename(zcta = zcta5ce10)


# Step 3: Extract patient zip code -------------------------------------------------------------------------

# Captures every instance with a five-digit zip code 
visits$zip <- str_replace(visits$patient_address, ".*\\b(\\d{5})\\b.*", "\\1")

# Repeats for nine-digit zip codes 
visits$zip <- str_replace(visits$zip, ".*\\b(\\d{9})\\b.*", "\\1")
 
# Restricts to a five-digit zip code 

visits <- visits %>% 
  mutate(zip = substr(zip,1,5)) %>% 
  select(-patient_address) 

# Takes out any spaces 

visits$zip <- str_replace_all(visits$zip, " ", "")

# Get a restricted zip code dataset

visits_clean_zip <- visits %>% 
  mutate(char_zip = nchar(zip)) %>% 
  filter(char_zip == 5) %>% 
  mutate(zip = as.numeric(zip)) %>%
  filter(!is.na(zip)) %>% 
  mutate(zip = as.character(zip)) %>%
  mutate(char_zip =nchar(zip)) %>%
  filter(char_zip == 5) %>%
  rename("zip_code" = "zip")

# Step 4: Combine ZCTAs to SF for mapping -------------------------------------------------------------------------

# Replace the zip_code of the addresses to the corresponding ZCTA 

zcta_visits <- left_join(visits_clean_zip, zip_xwalk, by="zip_code") %>% 
  select(-zip_code) %>% 
  mutate(zcta = as.character(zcta)) %>%
  filter(zcta %in% nyc_zctas$zcta)

write.fst(zcta_visits, paste0(out,"/","zcta_full_visits.fst"))

