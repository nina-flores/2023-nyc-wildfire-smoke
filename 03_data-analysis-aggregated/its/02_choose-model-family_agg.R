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
library(MASS)
library(tsModel)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"

dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated-agg.fst")) 

mod_main <- glm(respiratory~ smoke_day + time_elapsed_scaled +
                         as.factor(holiday) + harmonic(doy,4,365) +
                         rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled + covid_post  ,
                       data = dat,
                       family = poisson)

summary(mod_main)
testDispersion(mod_main) # bad in terms of dispersion, dispersion = 5.4968, p-value < 0.00000000000000022
simulationOutput <- simulateResiduals(fittedModel = mod_main, plot = F)

mod_main_nb<- glm.nb(respiratory ~ smoke_day + time_elapsed_scaled +
                  as.factor(holiday) + harmonic(doy,4,365) + dow +
                  rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled+covid_post  ,
                data = dat)

summary(mod_main_nb)

testDispersion(mod_main_nb) # better! dispersion = 0.88274, p-value < 0.00000000000000022
simulationOutput <- simulateResiduals(fittedModel = mod_main_nb, plot = F) # plot looks much better here too
plot(simulationOutput)



mod_main_nb1 <- glmmTMB(respiratory ~ smoke_day + time_elapsed_scaled + 
                       as.factor(holiday) + harmonic(doy,8,365) + dow +covid_acute+covid_post +
                       ns(rolling_avg_temperature_scaled, 5) + 
                         rolling_avg_precipitation_scaled  ,
                      family = nbinom1(),
                     data = dat)

summary(mod_main_nb1)

testDispersion(mod_main_nb1) # better! dispersion = 0.88274, p-value < 0.00000000000000022
simulationOutput <- simulateResiduals(fittedModel = mod_main_nb1, plot = F) # plot looks much better here too
plot(simulationOutput)

# will move forward with negative binomial model



# same for asthma ---------------------------------------------------------

mod_main <- glm(asthma ~ smoke_day + time_elapsed_scaled +
                  as.factor(holiday) + harmonic(doy,4,365) +
                  rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +covid_post  ,
                data = dat,
                family = poisson)

summary(mod_main)
testDispersion(mod_main) # 2.5347, p-value < 0.00000000000000022
simulationOutput <- simulateResiduals(fittedModel = mod_main, plot = F)
plot(simulationOutput)


mod_asthma <- glmmTMB(asthma ~ smoke_day + time_elapsed_scaled +
                  as.factor(holiday) + harmonic(doy,8,365)  +
                  ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                    covid_post + covid_acute ,
                    family = nbinom1(),
                data = dat)
summary(mod_asthma )
testDispersion(mod_asthma ) # 0.85623, p-value < 0.00000000000000022
simulationOutput <- simulateResiduals(fittedModel = mod_asthma,  plot = F)
plot(simulationOutput)

mod_asthma1 <- glmmTMB(asthma ~ smoke_day + time_elapsed_scaled +
                        as.factor(holiday) + harmonic(doy,8,365)  +
                        ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                        covid_post + covid_acute ,
                      family = nbinom1(),
                      data = dat)

mod_asthma2 <- glmmTMB(asthma ~ smoke_day + time_elapsed_scaled +
                         as.factor(holiday) + harmonic(doy,8,365)  +
                         ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                         covid_post + covid_acute ,
                       family = nbinom2(),
                       data = dat)


AIC(mod_asthma2, mod_asthma1)



# Conclusion from these analyses is that negative binomial distribution is the 
# way to go and in particular the regression diagnostics look best with nbinom1()

