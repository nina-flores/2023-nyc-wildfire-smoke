###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Prepare dataset for analysis by creating month, year, dow, and scaled variables

# we decided to remove 2020 since it is such an anomaly in terms of ED usage - for 
# example there are 23 times where the numbers of respiratory ED visits are 0
# in this data and they are all in 2020. 
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(tis)
library(lubridate)
library(zoo)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat_restrict <- read.fst(paste0(out,"/","analytical_data.fst")) %>%
  group_by(zcta) %>%
  summarize(
    sum_visits = sum(respiratory),
    sum_asthm = sum(asthma)) %>%
  ungroup() %>%
  mutate(sum_visits_all = sum(sum_visits)) %>%
  group_by(zcta) %>%
  mutate(percent_resp_visits = 100*(sum_visits/sum_visits_all)) #%>%
  #filter(percent_resp_visits >= .5) 
  

sum(dat_restrict$sum_visits)
sum(dat_restrict$sum_visits)



dat <- read.fst(paste0(out,"/","analytical_data.fst")) %>%
  filter(zcta %in% dat_restrict$zcta) %>%
  group_by(date) %>%
  summarize(
    respiratory = sum(respiratory),
    asthma = sum(asthma),
    total_population_under_18 = sum(total_population_under_18),
    #covid = sum(covid),
    caserate = mean(caserate), 
    respiratory_ed_inp = sum(respiratory_ed_inp),
    asthma_ed_inp = sum(asthma_ed_inp),
    temperature = mean(temperature),
    total_precipitation= mean(total_precipitation),
    ah = mean(ah),
    nevi_quartile = mean(nevi_quartile),
    nevi = mean(nevi),
    level_1 = sum(level_1),
    level_2 = sum(level_2),
    level_3 = sum(level_3),
    level_4 = sum(level_4),
    level_5 = sum(level_5)
    ) %>%
  mutate(respiratory_rate = (respiratory/total_population_under_18)*10000,
         asthma_rate = (asthma/total_population_under_18)*10000,
         rolling_avg_temperature = rollapply(temperature, width = 7, FUN = mean, partial = TRUE, align = "right"),
         rolling_avg_precipitation = rollapply(total_precipitation, width = 7, FUN = mean, partial = TRUE, align = "right"),
         level_1_rate = (level_1/total_population_under_18)*10000,
         level_2_rate = (level_2/total_population_under_18)*10000,
         level_3_rate = (level_3/total_population_under_18)*10000,
         level_4_rate = (level_4/total_population_under_18)*10000,
         level_5_rate = (level_5/total_population_under_18)*10000) %>%
  ungroup() %>%
  mutate(smoke_day = as.factor(if_else(date %in% c("2023-06-06","2023-06-07","2023-06-08") ,1,0)),
         smoke_day1 = as.factor(if_else(date %in% c("2023-06-06") ,1,0)),
         smoke_day2 = as.factor(if_else(date %in% c("2023-06-07") ,1,0)),
         smoke_day3 = as.factor(if_else(date %in% c("2023-06-08") ,1,0)),
         holiday = isHoliday(date),
         holiday = if_else(holiday == TRUE, 1,0),
         year = as.factor(year(date)),
         month = as.factor(month(date)),
         dow = as.factor(wday(date)),
         doy = yday(date),
         temperature_scaled = scale(temperature),
         total_precipitation_scaled = scale(total_precipitation),
         rolling_avg_temperature_scaled = scale(rolling_avg_temperature),
         rolling_avg_precipitation_scaled = scale(rolling_avg_precipitation),
         caserate_scaled = scale(caserate)) %>%
  mutate(post_treatment = if_else(date > as.Date("2023-06-08"),1,0),
         covid_acute = ifelse(date >= "2020-03-01" & date <= "2020-12-31", 1, 0),
         covid_post = ifelse(date >= "2022-01-01", 1, 0)) %>%
  arrange(date) %>%
  mutate(time_elapsed = row_number(),
         time_post_intervention = cumsum(post_treatment))  %>%
  ungroup() %>%
  mutate(time_elapsed_scaled = scale(time_elapsed))

sum(dat$respiratory)
# 91661/98494 # this method retains 93% of the visits

92831/99771

write.fst(dat, paste0(an ,"/","analytical_data_prepped-updated-agg.fst"))

