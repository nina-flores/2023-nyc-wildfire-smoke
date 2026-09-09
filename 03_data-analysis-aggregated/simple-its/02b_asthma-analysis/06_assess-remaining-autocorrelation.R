###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Determine temporal autocorrelation
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
source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst")) 

acf(dat$respiratory, na.action = na.pass, lag = 100) 
pacf(dat$respiratory, na.action = na.pass, lag = 100)


mod = glmer.nb(respiratory ~ smoke_day + 
                 time_elapsed_scaled + 
                 time_post_intervention +
                 year + 
                 month + 
                 dow + 
                 ns(temperature_scaled, df = 2) + 
                 ns(abs_hum_scaled, df = 2) +
                 total_precipitation_scaled + 
                 covid +
                 (1 |zcta),
               data = dat,
               nAGQ = 0)


# 2i.i Run ACF/PACF plots
acf(residuals(mod))
pacf(residuals(mod))

 # seems much improved