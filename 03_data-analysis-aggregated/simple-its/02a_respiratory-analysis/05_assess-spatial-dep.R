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

# read in zctas to get distance matrix ------------------------------------


# Get all ZCTAs as an sf object
zct <- zctas(state = "NY", year = 2010, class = "sf") %>%  # Filter for New York state
  mutate(zcta = ZCTA5CE10) %>%
  dplyr::select(zcta)


# add to full data --------------------------------------------------------

dat_spat <- dat %>% left_join(zct) %>%
  st_as_sf() 


centroids <- st_centroid(dat_spat$geometry)
centroid_coords <- st_coordinates(centroids)

dat_spat <- dat_spat %>%
  mutate(x = centroid_coords[,1],
         y = centroid_coords[,2])



# run the main model ------------------------------------------------------

mod = mod_main_RE_S_S = glmer(respiratory ~  smoke_day +  time_elapsed_scaled  + holiday +
                                year + month + dow + ns(rolling_avg_temperature_scaled, 5) + rolling_avg_precipitation_scaled +
                                factor(borough) + ns(x,3) + ns(y,3) +
                                (1 | zcta),
                              data = dat_spat,
                              family = "poisson",
                              nAGQ = 0)
summary(mod)

# test spatial autocorrelation, conditional on RE's -----------------------

# double check with recalculated residuals 
res <- simulateResiduals(mod)
res2 <- recalculateResiduals(res, group = dat_spat$zcta)

testSpatialAutocorrelation(res2,
                           x = aggregate(dat_spat $x, list(dat$zcta), mean)$x,
                           y = aggregate(dat_spat $y, list(dat$zcta), mean)$x)

# There were spatial autocorrelation issues that go away with adding lat/long
# splines to the model


