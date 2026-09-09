###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Determine spatial autocorrelation
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(lme4)
library(DHARMa)
library(splines)
library(spdep)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated.fst")) 

# run the main model ------------------------------------------------------

mod = mod_main_RE_S_S = glmer.nb(asthma ~  smoke_day * time_post_intervention + time_elapsed_scaled  + holiday +
                                   year + month + dow + ns(rolling_avg_temperature_scaled, 5) + rolling_avg_precipitation_scaled +
                                   (1 | zcta),
                                 data = dat,
                                 nAGQ = 0)
summary(mod)
exp(0.360218)
