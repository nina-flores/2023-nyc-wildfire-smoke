###########***********************
#### Code Description ####
# Author: Nina
# Goal: Test different spline df combinations - NEGATIVE BINOMIAL
####**********************
####*

# Load packages
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(MASS)      # glm.nb
library(DHARMa)
library(splines)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# Read in data -------------------------------------------------------------------------
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated-agg.fst"))


mod<- glmmTMB(respiratory ~ smoke_day + time_elapsed_scaled +
                        as.factor(holiday) + harmonic(doy,8,365)  +
                        ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                        covid_post + covid_acute ,
                      family = nbinom1(),
                      data = dat)


# Pearson residuals from the glm.nb fit
resid_pearson <- residuals(mod, type = "pearson")

# Visual check
acf(resid_pearson, main = "ACF of Pearson residuals - NB model")
pacf(resid_pearson, main = "PACF of Pearson residuals - NB model")

# Formal tests at a few relevant lags
Box.test(resid_pearson, lag = 7,   type = "Ljung-Box")   # weekly
Box.test(resid_pearson, lag = 30,  type = "Ljung-Box")   # monthly
Box.test(resid_pearson, lag = 365, type = "Ljung-Box")   # annual



# same for asthma ---------------------------------------------------------

mod<- glmmTMB(asthma ~ smoke_day + time_elapsed_scaled +
                as.factor(holiday) + harmonic(doy,10,365)  +
                ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                covid_post + covid_acute ,
              family = nbinom1(),
              data = dat)



# Pearson residuals from the glm.nb fit
resid_pearson <- residuals(mod, type = "pearson")

# Visual check
acf(resid_pearson, main = "ACF of Pearson residuals - NB model")
pacf(resid_pearson, main = "PACF of Pearson residuals - NB model")

# Formal tests at a few relevant lags
Box.test(resid_pearson, lag = 7,   type = "Ljung-Box")   # weekly
Box.test(resid_pearson, lag = 30,  type = "Ljung-Box")   # monthly
Box.test(resid_pearson, lag = 365, type = "Ljung-Box")   # annual

