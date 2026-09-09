###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Prepare dataset for analysis by creating month, year, dow, and scaled variables
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)

source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))

set.seed(444)

# read in data -------------------------------------------------------------------------

cov <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/cov"
pro <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates"


covid <- fread( paste0(cov ,"/","caserate-by-modzcta.csv")) %>%
  select(c("week_ending", "CASERATE_BX","CASERATE_MN"))%>%
  pivot_longer(
    cols = starts_with("CASERATE_"), # Select columns to pivot
    names_to = "county", # New column for mzcta values
    names_prefix = "CASERATE_", # Remove this prefix from column names
    values_to = "caserate" # New column for the case rate values
  )

# Get all ZCTAs as an sf object
zct <- zctas(state = "NY", year = 2010, class = "sf") %>%  # Filter for New York state
  mutate(zcta = ZCTA5CE10) %>%
  select(zcta)  %>%
  filter(grepl("^100|^104", zcta)) %>%
  mutate(county = if_else(grepl("^100", zcta), "MN", "BX"))

# join zcta information for joining ---------------------------------------

cov_zct <- covid %>% 
  full_join(zct) %>%
  mutate(date = as.Date(week_ending, format = "%m/%d/%Y")) %>%
  as.data.frame() %>%
  select(-c(week_ending, county, geometry))


write.fst(cov_zct, paste0(pro,"/","processed_covid.fst"))

