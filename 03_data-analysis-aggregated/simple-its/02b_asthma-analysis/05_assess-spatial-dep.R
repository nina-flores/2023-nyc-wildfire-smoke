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


dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst")) 

# read in zctas to get distance matrix ------------------------------------


# Get all ZCTAs as an sf object
zct <- zctas(state = "NY", year = 2010, class = "sf") %>%  # Filter for New York state
  mutate(zcta = ZCTA5CE10) %>%
  dplyr::select(zcta)

# make new grouping based on what we're using for the data:
zct_100.4 <- read.fst( paste0(an ,"/","zct_1004.fst"))
zct_104.4 <- read.fst( paste0(an ,"/","zct_1044.fst"))
zct_100.1 <- read.fst( paste0(an ,"/","zct_1001.fst"))
zct_104.1 <- read.fst( paste0(an ,"/","zct_1041.fst"))


zct_mod <- zct %>%
  mutate(zcta = if_else(zcta %in% zct_100.4$zcta, "10000.4", zcta),
         zcta = if_else(zcta %in% zct_104.4$zcta, "10400.4", zcta),
         zcta = if_else(zcta %in% zct_100.1$zcta, "10000.1", zcta),
         zcta = if_else(zcta %in% zct_104.1$zcta, "10400.1", zcta)) %>%
  group_by(zcta) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup() 

  
# add to full data --------------------------------------------------------

dat_spat <- dat %>% left_join(zct_mod) %>%
  st_as_sf() 

#plot(dat_spat["geometry"])

#nb <- poly2nb(dat_spat, queen = TRUE)
#weights_list <- nb2listw(nb, style = "B")

centroids <- st_centroid(dat_spat$geometry)
centroid_coords <- st_coordinates(centroids)

dat_spat <- dat_spat %>%
  mutate(x = centroid_coords[,1],
         y = centroid_coords[,2])


# run the main model ------------------------------------------------------

mod = mod_main_RE_S_S = glmer(asthma ~ smoke_day + time_post_intervention + time_elapsed_scaled  + holiday +
                                year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                                (1 | zcta),
                              data = dat,
                              family = "poisson",
                              offset = log(total_population_under_18), 
                              nAGQ = 0)
summary(mod)
ranef(mod)
# test spatial autocorrelation, conditional on RE's -----------------------

# double check with recalculated residuals 
res <- simulateResiduals(mod)
res2 <- recalculateResiduals(res, group = dat_spat$zcta)

testSpatialAutocorrelation(res2,
                           x = aggregate(dat_spat $x, list(dat$zcta), mean)$x,
                           y = aggregate(dat_spat $y, list(dat$zcta), mean)$x)




# run the main model ------------------------------------------------------

mod = glmer(asthma ~ smoke_day + time_post_intervention + time_elapsed_scaled  + holiday +
                                year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                                borough +
                                (1 | zcta),
                              data = dat,
                              family = "poisson",
                              offset = log(total_population_under_18), 
                              nAGQ = 0)
summary(mod)
ranef(mod)
# test spatial autocorrelation, conditional on RE's -----------------------

# double check with recalculated residuals 
res <- simulateResiduals(mod)
res2 <- recalculateResiduals(res, group = dat_spat$zcta)

testSpatialAutocorrelation(res2,
                           x = aggregate(dat_spat $x, list(dat$zcta), mean)$x,
                           y = aggregate(dat_spat $y, list(dat$zcta), mean)$x)


# including borough terms largely reduces spatial autocorrelation. We could 
# do better probably by including a distance from nyp term in the future. 


