###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/10/24
# Goal: Join Covariate data to Patient ZCTAs
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# Step 1: join data -------------------------------------------------------------------------

wd <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/RITM0457513"
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
cov <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates"


visits <- read.fst(paste0(out,"/","zcta_day_outcomes.fst")) 




covariates <- read.fst(paste0(cov,"/","wf-covariates.fst")) %>%
  group_by(zcta, date) %>%
  slice(1) %>%
  ungroup()



visit_variates <- full_join(covariates, visits) %>%
  drop_na(temperature) %>%
  mutate(respiratory = if_else(is.na(respiratory), 0, respiratory),
         asthma = if_else(is.na(asthma), 0, asthma),
         respiratory_atleast5 = if_else(is.na(respiratory_atleast5), 0, respiratory_atleast5),
         asthma_atleast5 = if_else(is.na(asthma_atleast5), 0, asthma_atleast5),
         covid = if_else(is.na(asthma), 0, covid),
         respiratory_ed_inp = if_else(is.na(respiratory_ed_inp), 0, respiratory_ed_inp),
         asthma_ed_inp = if_else(is.na(asthma_ed_inp), 0, asthma_ed_inp))




write.fst(visit_variates, paste0(out,"/","analytical_data.fst"))
