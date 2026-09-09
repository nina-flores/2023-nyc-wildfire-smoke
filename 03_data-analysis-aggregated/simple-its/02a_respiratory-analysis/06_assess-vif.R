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
#remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
require(glmmTMB)
require(sandwich)
require(clubSandwich)
require(tsModel)

source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated.fst")) 

mod = mod_main_RE_S_S = glmer(respiratory ~  smoke_day +  time_elapsed_scaled  + holiday +
                                month + dow + ns(rolling_avg_temperature_scaled, 5) + rolling_avg_precipitation_scaled +
                                 ns(x,3) + ns(y,3) +
                                (1 | zcta),
                              data = dat,
                              family = "poisson",
                              nAGQ = 0)

summary(mod)

model_2 <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled + 
                   harmonic(doy,2,365) + dow + holiday + rolling_avg_precipitation_scaled +ns(x,3) + ns(y,3)+(1 | zcta),
                 data = dat,
                 family = "poisson",
                 nAGQ = 0)

model_3 <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled + 
                   harmonic(doy,3,365) + dow + holiday + rolling_avg_precipitation_scaled+ns(x,3) + ns(y,3) +(1 | zcta),
                 data = dat,
                 family = "poisson",
                 nAGQ = 0)

model_4 <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled + 
               harmonic(doy,4,365) + dow + holiday + rolling_avg_precipitation_scaled+ns(x,3) + ns(y,3) +(1 | zcta),
             data = dat,
             family = "poisson",
             nAGQ = 0)

model_5 <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled + 
                   harmonic(doy,5,365) + dow + holiday + rolling_avg_precipitation_scaled+ns(x,3) + ns(y,3) +(1 | zcta),
                 data = dat,
                 family = "poisson",
                 nAGQ = 0)


model_6 <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled + 
                   harmonic(doy,6,365) + dow + holiday + rolling_avg_precipitation_scaled+ns(x,3) + ns(y,3) +(1 | zcta),
                 data = dat,
                 family = "poisson",
                 nAGQ = 0)


car::vif(mod)
car::vif(model_2)
car::vif(model_3)
car::vif(model_4)
car::vif(model_5)
car::vif(model_6)



AIC(model_2)
AIC(model_3)
AIC(model_4) # this one
AIC(model_5)
AIC(model_6)

# vif is really high when using the month and year indicators but not when using 
# harmonics for season instead


summary(model_4)



model_mult <- glmer(respiratory ~ smoke_day1 + smoke_day2 + smoke_day3 +  time_elapsed_scaled +ns(x,3) + ns(y,3) +
                   harmonic(doy,4,365) + dow + holiday + rolling_avg_precipitation_scaled+ caserate_scaled +(1 | zcta),
                 data = dat,
                 family = "poisson",
                 nAGQ = 0)

summary(model_mult)



model <- glmer(respiratory ~ smoke_day +  time_elapsed_scaled +ns(x,3) + ns(y,3) +
                        harmonic(doy,4,365) + dow + holiday + rolling_avg_precipitation_scaled+ caserate_scaled +(1 | zcta),
                      data = dat,
                      family = "poisson",
                      nAGQ = 0)

summary(model)

exp(0.306169)





# model_ns_mult <- glmer(respiratory ~ smoke_day1 + smoke_day2 + smoke_day3 +  ns(time_elapsed_scaled,7.5*4) + 
#                          dow + holiday + rolling_avg_precipitation_scaled + caserate_scaled+ (1 | zcta),
#                       data = dat,
#                       family = "poisson",
#                       nAGQ = 0)
# 
# summary(model_ns_mult)
# 
# 
# attr(ns(dat$time_elapsed_scaled, 12), "knots")
# 
# 
#  # recover the scaling parameters used to create time_elapsed_scaled
#    # (assuming you did something like scale(time_elapsed) originally)
#    center <- attr(scale(dat$time_elapsed), "scaled:center")
#  scale_ <- attr(scale(dat$time_elapsed), "scaled:scale")
#  
#    knots_scaled <- attr(ns(dat$time_elapsed_scaled, 12), "knots")
#  knots_unscaled <- knots_scaled * scale_ + center
#  
#    # if time_elapsed is "days since study start" (or similar), convert to dates:
#    study_start <- min(dat$date)
#  knot_dates <- study_start + knots_unscaled
# 
#  knot_dates
#  
#  attr(ns(dat$time_elapsed_scaled, 12), "Boundary.knots")
# 
#  
#  
#  boundary_scaled <- attr(ns(dat$time_elapsed_scaled, 12), "Boundary.knots")
#  boundary_unscaled <- boundary_scaled * scale_ + center
#  boundary_dates <- study_start + boundary_unscaled
#  boundary_dates
#  
#  
#  
#  