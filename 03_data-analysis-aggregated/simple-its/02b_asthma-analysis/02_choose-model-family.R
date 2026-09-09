###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Determine choice of model family
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(lme4)
library(DHARMa)
library(splines)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst")) 

# model with random effects for zcta
mod_main_RE = glmer(asthma ~ smoke_day + time_elapsed_scaled + time_post_intervention + holiday +
                      year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation +
                      (1|zcta),
                    data = dat,
                    family = "poisson",
                    offset = log(total_population_under_18), 
                    nAGQ = 0)

summary(mod_main_RE)
simulationOutput <- simulateResiduals(fittedModel = mod_main_RE, plot = F)
testDispersion(simulationOutput) # There is no overdispersion, 0.22913, p-value = 0.52
plot(simulationOutput)
testZeroInflation(simulationOutput) # no zero inflation either 1.06, p-value = 0.544

# test negative binomial
mod_main_RE_nb = glmer.nb(asthma ~ smoke_day + time_elapsed_scaled + time_post_intervention + holiday +
                            year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation +
                            (1|zcta),
                          data = dat,
                          offset = log(total_population_under_18), 
                          nAGQ = 0)


summary(mod_main_RE_nb)
simulationOutput_nb <- simulateResiduals(fittedModel = mod_main_RE_nb, plot = F)
testDispersion(mod_main_RE_nb) # dispersion = 0.20095, p-value = 0.368
plot(simulationOutput_nb)
testZeroInflation(simulationOutput_nb) # ratioObsSim = 0.81303, p-value = 0.376
testUniformity(mod_main_RE_nb) 
testOutliers(mod_main_RE_nb) 

# Will move forward with poisson. 

