###########***********************
#### Code Description ####
# Author: Alex
# Date: 7/5/24
# Updated: 7/5/24
# Goals: 
# A) Extract Medication Data from Epic and Pre-Epic Data; 
# B) Create an ICS Indicator for Asthma Medications
####**********************

# Load packages 
library(here)
library(tigris)
library(datetimeutils)
library(stringr)
library(data.table)
source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))

# Step 1: Load data required-------------------------------------------------------------------------

wd <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/RITM0457513"
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"


# Loads in Epic meds and cleans up column names

epic_meds <- read_csv(paste0(wd,"/","RITM0457513 - Epic Medication Administrations.csv"), col_names = FALSE) %>%
  janitor::clean_names() %>% 
  rename("mrn" = x2,
         "account" = x1,
         "description" = x3, 
         "arrival" = x4, # This might not be accurate; please confirm with data documentation 
         "administration_date" = x5, # This might not be accurate; please confirm with data documentation 
         "name" = x6,
         "dose" = x7,
         "route" = x8,
         "unit" = x9,
         "location" = x10) %>%
  select(mrn, account, description) %>%
  mutate(from = "Epic") #%>%
 # mutate(arrival = as_date(with_tz(administration_date, tzone = "America/New_York"))) %>%
 # select(-administration_date)

pre_meds <- read_csv(paste0(wd,"/","RITM0457513 - Pre-Epic Medication Administrations.csv")) %>% 
  separate_wider_delim("ACCOUNT|ORDER_MED_ID|MRN|DESCRIPTION|ORDERING_DATE|ADMINISTRATION_DATE|ADMINISTRATING_USER|DOSE|DOSE_UNIT|ROUTE",
                        delim = "|",
           names = c("account", "order_med_id","mrn","description",
                    "arrival","administration_date",
                    "administrating_user_dose",
                    "dose", "unit", "route"),
           too_few = 'debug') %>%
  select(mrn, account, description) %>%
  mutate(from = "Pre-Epic")

# identify level of visit -------------------------------------------------


medications <- rbind(epic_meds, pre_meds) %>%
  unique()

#rm(epic_meds)
#rm(pre_meds)

# Define medication lists
saba_list <- c("albuterol", "levalbuterol") # AJN - can remove pirbuterol and levosalbutamol

inhaled_cs_list <- c("budesonide", "fluticasone", "mometasone", "beclomethasone", "ciclesonide") 

systemic_cs_list <- c("dexamethasone", "prednisone", "methylprednisolone", "prednisolone") # AJN - added prednisolone 

terb_lama_mg_list <- c("magnesium", "ipratropium", "tiotropium", "aclidinium", "glycopyrrolate", "umeclidinium", "terbutaline")
# Create a single regex pattern
saba_pattern <- str_c(saba_list, collapse = "|") # Level 1 
i_cs_pattern <- str_c(inhaled_cs_list, collapse = "|") # Level 2 
s_cs_pattern <- str_c(systemic_cs_list, collapse = "|") # Level 3 
terb_lama_mg_pattern <- str_c(terb_lama_mg_list, collapse = "|") # Level 4
# Level 5 -> Hospitalization


# to data table
setDT(medications)

# Apply function to determine levels
medications[, medications_lower := tolower(description)]
medications[, saba := as.integer(str_detect(medications_lower, saba_pattern))]
medications[saba == 1 & str_detect(medications_lower, "if albuterol"), saba := 0] # removing it if it is actually another medication but says "if albuterol also prescribed ..."
medications[, scs := as.integer(str_detect(medications_lower, s_cs_pattern))]
medications[, ics := as.integer(str_detect(medications_lower, i_cs_pattern))]
medications[, other := as.integer(str_detect(medications_lower, terb_lama_mg_pattern))]

# use max to apply the code to all the entries for the same encounter
medications[, c("max_saba", "max_scs", "max_ics", "max_other") := .(max(saba, na.rm = TRUE), 
                                                                    max(ics, na.rm = TRUE),
                                                                    max(scs, na.rm = TRUE), 
                                                                    max(other, na.rm = TRUE)), by = c("account")]

# Apply the conditions for levels outlined by alex - level  4 to be assigned when added back to full dataset. 
medications <- medications %>%
  mutate(medication_group = case_when(
    max_saba == 1 & max_ics == 0 & max_scs == 0  & max_other == 0 ~ 1,
    max_saba == 1 & max_ics == 1 & max_scs == 0 & max_other == 0 ~ 2,
    max_saba == 1  & max_scs == 1 & max_other == 0 ~ 3,
    max_saba == 1  & max_other == 1 ~ 4,
    max_saba == 0 & max_scs == 0 & max_ics == 0 & max_other == 0 ~ 0,
    TRUE ~ 0  
  ))

# Q for Alex/Stephanie
# would there be any reason for an asthma case that CS would be prescribed without SABA? If so, may need to update code.


meds <- medications %>%
  select(account, medication_group) %>%
  unique()

##meds_pre_epic <- medications %>%
##  filter(from == "Pre-Epic") %>%
##  select(mrn, account, medication_group) %>%
##  unique()


write.fst(meds, paste0(out,"/","medi-level.fst"))
#write.fst(meds_pre_epic, paste0(out,"/","medi-level-pre-epic.fst"))


