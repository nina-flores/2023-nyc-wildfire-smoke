###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Determine choice of main model 
#    (random intercepts OR random intercepts and random slopes)
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(lme4)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated.fst")) 

str(dat)

# model with random effects for zcta

mod_main_RE = glmer(respiratory ~ smoke_day + time_elapsed_scaled  + holiday +
                            year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                            (1|zcta),
                          data = dat,
                    family = "poisson", 
                    nAGQ = 0)
summary(mod_main_RE)

# model with random intercept and smoke slope
mod_main_RE_S = glmer(respiratory ~  time_elapsed_scaled  + holiday +
                        year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                         (1 + smoke_day |zcta),
                      data = dat,
                      family = "poisson", 
                      nAGQ = 0)
summary(mod_main_RE_S)

AIC(logLik(mod_main_RE))                                    
AIC(logLik(mod_main_RE_S))                                  
anova(mod_main_RE, mod_main_RE_S)



AIC(logLik(mod_main_RE))                                    # AIC = 195240.8
AIC(logLik(mod_main_RE_S))                                  # AIC = 195242.5

### AIC gets worse with the slopes, keep the simpler model.


