###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/10/24
# Goals: 
# A) summarize all outcomes of interest
# B) expand grid to all zctas and dates covered
####**********************

# Load packages 
library(here)
library(tigris)
library(datetimeutils)
source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))


# Step 1: Load data required -------------------------------------------------------------------------

out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"


# List of respiratory prefixes to check
prefixes <- c("J45", "J00", "J02", "J04", "J05", "J06",
              "J21", "R05")

chd_prefixes <- c("Q20", "Q21", "Q22", "Q23", "Q24", "Q25", "Q26", "Q27", "Q28")

nm_prefixes <- c("G12", "G13", "G14", "G60", "G61", "G70", "G71", "G72", "G73")

# Create a regex pattern to match any of the prefixes
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
chd_pattern <- paste0("^(", paste(chd_prefixes, collapse = "|"), ")")
nm_pattern <- paste0("^(", paste(nm_prefixes, collapse = "|"), ")")

meds <- read.fst(paste0(out,"/","medi-level.fst")) %>%
  rename(pat_enc_csn_id = account) %>%
  mutate(pat_enc_csn_id = as.numeric(pat_enc_csn_id))

visits <- read.fst(paste0(out,"/","zcta_full_visits.fst")) %>%
  left_join(meds) %>%
  mutate(medication_group = if_else(is.na(medication_group),0,medication_group)) %>%
  mutate(date = as_date(arrival)) %>%
  mutate(cf_dx_yn = ifelse(str_detect(dx_code, "^E84") & date < as_date("2023-06-06"), 1, 0),
         scd_dx_yn = ifelse(str_detect(dx_code, "^D57") & date < as_date("2023-06-06"), 1, 0),
         chd_dx_yn = ifelse(str_detect(dx_code, chd_pattern), 1, 0),   
         nm_dx_yn = ifelse(str_detect(dx_code, nm_pattern), 1, 0)) %>%
  group_by(mrn) %>%
  mutate(cf_dx_yn = max(cf_dx_yn),
         scd_dx_yn = max(scd_dx_yn),
         chd_dx_yn = max(chd_dx_yn),
         nm_dx_yn = max(nm_dx_yn)) %>%
  ungroup() %>%
  select(-char_zip) %>%
  filter(class != "Inpatient") %>%
  mutate(respiratory_dx_yn = ifelse(str_detect(dx_code, pattern), 1, 0)) %>%
  mutate(asthma_dx_yn = ifelse(str_detect(dx_code, "^J45"), 1, 0)) %>%
  mutate(respiratory_dx_yn_ed_inp = ifelse(str_detect(dx_code, pattern) & class == "ED to Inpatient", 1, 0)) %>%
  mutate(asthma_dx_yn_ed_inp = ifelse(str_detect(dx_code, "^J45") & class == "ED to Inpatient", 1, 0)) %>%
  mutate(covid_dx_yn = ifelse(str_detect(dx_code, "^U07.1" ), 1, 0))%>%
  mutate(respiratory_dx_yn_atleast5 = ifelse(str_detect(dx_code, pattern) & age >=5, 1, 0)) %>%
  mutate(asthma_dx_yn_atleast5 = ifelse(str_detect(dx_code, "^J45")& age >=5, 1, 0)) %>%
  mutate(visits = 1) %>%
  mutate(medication_group = if_else(respiratory_dx_yn_ed_inp == 1 | asthma_dx_yn_ed_inp, 5, medication_group))

# unique individuals in the data set from ED visits
unique_ind <- visits %>%
  filter(respiratory_dx_yn == 1) %>%
  select(mrn, sex, ethnicity, age) %>%
  group_by(mrn) %>%
  dplyr::slice(1)

write.fst(unique_ind, paste0(out,"/","unique_ind.fst"))

###visits_cf <-   visits %>%
###  filter(cf_dx_yn == 1) 
###
###visits_scd <-   visits %>%
###  filter(scd_dx_yn == 1) 

# now get this to unique ED visits
visits_unique <- visits %>%
  select(-dx_code) %>%
  group_by(pat_enc_csn_id) %>%
  mutate(respiratory_dx_yn = max(respiratory_dx_yn),
         asthma_dx_yn = max(asthma_dx_yn),
         respiratory_dx_yn_ed_inp = max(respiratory_dx_yn_ed_inp),
         covid_dx_yn = max(covid_dx_yn),
         asthma_dx_yn_ed_inp = max(asthma_dx_yn_ed_inp),
         cf_dx_yn = first(cf_dx_yn),
         scd_dx_yn = first(scd_dx_yn),
         chd_dx_yn = first(chd_dx_yn),
         nm_dx_yn = first(nm_dx_yn),
         total_visits = max(visits)) %>%
  dplyr::slice(1)

write.fst(visits_unique, paste0(out,"/","total_visits.fst"))



# full data ---------------------------------------------------------------
# now zcta-day aggregation
zcta_day_visits <- visits_unique %>%
  group_by(zcta, date) %>%
  summarize(respiratory = sum(respiratory_dx_yn),
            asthma = sum(asthma_dx_yn),
            respiratory_atleast5 = sum(respiratory_dx_yn_atleast5),
            asthma_atleast5 = sum(asthma_dx_yn_atleast5),
            respiratory_ed_inp = sum(respiratory_dx_yn_ed_inp),
            asthma_ed_inp = sum(asthma_dx_yn_ed_inp),
            covid = sum(covid_dx_yn),
            total_visits = sum(total_visits),
            other_visits = total_visits - respiratory,
            level_1 = sum(medication_group == 1, na.rm = TRUE),
            level_2 = sum(medication_group == 2, na.rm = TRUE),
            level_3 = sum(medication_group == 3, na.rm = TRUE),
            level_4 = sum(medication_group == 4, na.rm = TRUE),
            level_5 = sum(medication_group == 5, na.rm = TRUE))
            
            

# expand the grid so that 0's are filled in properly

# Create a sequence of all dates within the range
all_dates <- seq.Date(as.Date("2016-01-01"), as.Date("2023-08-31"), by = "day")

# Create a data frame of all combinations of zcta and dates
full_grid <- expand.grid(zcta = unique(visits_unique$zcta), date = all_dates)

# Merge the full grid with the original dataset
df_expanded <- full_grid %>%
  left_join(zcta_day_visits, by = c("zcta", "date")) %>%
  mutate(respiratory = if_else(is.na(respiratory), 0, respiratory),
         asthma = if_else(is.na(asthma), 0, asthma),
         respiratory_atleast5 = if_else(is.na(respiratory_atleast5), 0, respiratory_atleast5),
         asthma_atleast5 = if_else(is.na(asthma_atleast5), 0, asthma_atleast5),
         respiratory_ed_inp = if_else(is.na(respiratory_ed_inp), 0, respiratory_ed_inp),
         asthma_ed_inp = if_else(is.na(asthma_ed_inp), 0, asthma_ed_inp),
         covid = if_else(is.na(covid), 0, covid),
         total_visits = if_else(is.na(total_visits), 0, total_visits),
         other_visits = if_else(is.na(other_visits), 0, other_visits),
         level_1 = if_else(is.na(level_1), 0, level_1),
         level_2 = if_else(is.na(level_2), 0, level_2),
         level_3 = if_else(is.na(level_3), 0, level_3),
         level_4 = if_else(is.na(level_4), 0, level_4),
         level_5 = if_else(is.na(level_5), 0, level_5))

write.fst(df_expanded , paste0(out,"/","zcta_day_outcomes.fst"))



# cf only -----------------------------------------------------------------
# now zcta-day aggregation
zcta_day_visits_cf <- visits_unique %>%
  filter(cf_dx_yn == 1) %>%
  group_by(zcta, date) %>%
  summarize(cf_respiratory = sum(respiratory_dx_yn),
            cf_asthma = sum(asthma_dx_yn))

# expand the grid so that 0's are filled in properly

# Create a sequence of all dates within the range
all_dates <- seq.Date(as.Date("2016-01-01"), as.Date("2023-08-31"), by = "day")

# Create a data frame of all combinations of zcta and dates
full_grid <- expand.grid(zcta = unique(visits_unique$zcta), date = all_dates)

# Merge the full grid with the original dataset
df_expanded_cf <- full_grid %>%
  left_join(zcta_day_visits_cf, by = c("zcta", "date")) %>%
  mutate(cf_respiratory = if_else(is.na(cf_respiratory), 0, cf_respiratory),
         cf_asthma = if_else(is.na(cf_asthma), 0, cf_asthma))%>%
  select(zcta, date, cf_respiratory, cf_asthma)

write.fst(df_expanded_cf , paste0(out,"/","zcta_day_outcomes_cf.fst"))


# scd only ----------------------------------------------------------------
zcta_day_visits_scd <- visits_unique %>%
  filter(scd_dx_yn == 1) %>%
  group_by(zcta, date) %>%
  summarize(scd_respiratory = sum(respiratory_dx_yn),
            scd_asthma = sum(asthma_dx_yn),
)

# expand the grid so that 0's are filled in properly

# Create a sequence of all dates within the range
all_dates <- seq.Date(as.Date("2016-01-01"), as.Date("2023-08-31"), by = "day")

# Create a data frame of all combinations of zcta and dates
full_grid <- expand.grid(zcta = unique(visits_unique$zcta), date = all_dates)

# Merge the full grid with the original dataset
df_expanded_scd <- full_grid %>%
  left_join(zcta_day_visits_scd, by = c("zcta", "date")) %>%
  mutate(scd_respiratory = if_else(is.na(scd_respiratory), 0, scd_respiratory),
         scd_asthma = if_else(is.na(scd_asthma), 0, scd_asthma)) %>%
  select(zcta, date, scd_respiratory, scd_asthma)

write.fst(df_expanded_scd , paste0(out,"/","zcta_day_outcomes_scd.fst"))


# chd only ----------------------------------------------------------------

zcta_day_visits_chd <- visits_unique %>%
  filter(chd_dx_yn == 1) %>%
  group_by(zcta, date) %>%
  summarize(chd_respiratory = sum(respiratory_dx_yn),
            chd_asthma = sum(asthma_dx_yn),
  )

# expand the grid so that 0's are filled in properly

# Create a sequence of all dates within the range
all_dates <- seq.Date(as.Date("2016-01-01"), as.Date("2023-08-31"), by = "day")

# Create a data frame of all combinations of zcta and dates
full_grid <- expand.grid(zcta = unique(visits_unique$zcta), date = all_dates)

# Merge the full grid with the original dataset
df_expanded_chd <- full_grid %>%
  left_join(zcta_day_visits_chd, by = c("zcta", "date")) %>%
  mutate(chd_respiratory = if_else(is.na(chd_respiratory), 0, chd_respiratory),
         chd_asthma = if_else(is.na(chd_asthma), 0, chd_asthma)) %>%
  select(zcta, date, chd_respiratory, chd_asthma)

write.fst(df_expanded_chd , paste0(out,"/","zcta_day_outcomes_chd.fst"))


# nm only ----------------------------------------------------------------

zcta_day_visits_nm <- visits_unique %>%
  filter(nm_dx_yn == 1) %>%
  group_by(zcta, date) %>%
  summarize(nm_respiratory = sum(respiratory_dx_yn),
            nm_asthma = sum(asthma_dx_yn),
  )

# expand the grid so that 0's are filled in properly

# Create a sequence of all dates within the range
all_dates <- seq.Date(as.Date("2016-01-01"), as.Date("2023-08-31"), by = "day")

# Create a data frame of all combinations of zcta and dates
full_grid <- expand.grid(zcta = unique(visits_unique$zcta), date = all_dates)

# Merge the full grid with the original dataset
df_expanded_nm <- full_grid %>%
  left_join(zcta_day_visits_nm, by = c("zcta", "date")) %>%
  mutate(nm_respiratory = if_else(is.na(nm_respiratory), 0, nm_respiratory),
         nm_asthma = if_else(is.na(nm_asthma), 0, nm_asthma)) %>%
  select(zcta, date, nm_respiratory, nm_asthma)

write.fst(df_expanded_nm , paste0(out,"/","zcta_day_outcomes_nm.fst"))


# join datasets together --------------------------------------------------

data <- df_expanded %>%
  full_join(df_expanded_cf) %>%
  full_join(df_expanded_scd) %>%
  full_join(df_expanded_chd) %>%
  full_join(df_expanded_nm) 


dat_smoke <- data %>%
  filter(date >= "2023-06-06" & date <= "2023-06-29")


sum(dat_smoke$scd_respiratory) # during smoke 1, 23 if we look at any time 6/6/2023 onward
# and 10 if we look 6/6-6/29

sum(dat_smoke$cf_respiratory) # during smoke 0, 2 if we look at any time 6/6/2023 onward,
# and 0 if we look 6/6-6/29

sum(dat_smoke$chd_respiratory) # during smoke 1, 85 if we look at any time 6/6/2023 onward
# and 22 if we look 6/6-6/29

sum(dat_smoke$nm_respiratory) # during smoke 1, 2 if we look at any time 6/6/2023 onward
# and 1 if we look 6/6-6/29


sum(zcta_day_visits_scd$scd_respiratory) # total 1692
sum(zcta_day_visits_cf$cf_respiratory) # total 108
sum(zcta_day_visits_chd$chd_respiratory) # total 3879
sum(zcta_day_visits_nm$nm_respiratory) # total 254

# doesnt seem like we can actually do anything with the cf and scd data

