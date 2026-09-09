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


dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst")) 

str(dat)

# model with random effects for zcta

mod_main_RE = glmer(asthma ~ smoke_day + time_elapsed_scaled + time_post_intervention + holiday +
                            year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                            (1|zcta),
                          data = dat,
                          family = "poisson",
                          offset = log(total_population_under_18), 
                          nAGQ = 0)
summary(mod_main_RE)

# model with random intercept and smoke slope
mod_main_RE_S = glmer(asthma ~  time_elapsed_scaled + time_post_intervention + holiday +
                        year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                         (1 + smoke_day |zcta),
                      data = dat,
                      family = "poisson",
                      offset = log(total_population_under_18), 
                      nAGQ = 0)
summary(mod_main_RE_S)

AIC(logLik(mod_main_RE))                                    # AIC = 139383 
AIC(logLik(mod_main_RE_S))                                  # AIC = 139382
anova(mod_main_RE, mod_main_RE_S)


# model with random intercept, smoke slope, and time since intervention slope
mod_main_RE_S_S = glmer(asthma ~  time_elapsed_scaled  + holiday +
                        year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                        (1 + smoke_day + time_post_intervention | zcta),
                      data = dat,
                      family = "poisson",
                      offset = log(total_population_under_18), 
                      nAGQ = 0)
summary(mod_main_RE_S_S)


AIC(logLik(mod_main_RE))                                    # AIC = 81966.44 
AIC(logLik(mod_main_RE_S))                                  # AIC = 81971.36
AIC(logLik(mod_main_RE_S_S))                                  # AIC = 81980.11

### We have a worse aic with the random slopes but are interested particularly 
### in these between zcta differences and will move forward with it.


