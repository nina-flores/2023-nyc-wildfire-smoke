# The purpose of this script is to use tidycensus to pull in the census
# zcta level varaible of total population under 18 for each year

require(tidyverse)
require(tidycensus)
require(fst)
census_api_key("b30c02786d1c9fddff73822ffe0dd6c8dd31b82e")

pro <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates"

# pull in data ------------------------------------------------------------

v22 <- load_variables(2022, "acs5", cache = TRUE)

v20_dec <- load_variables(2020, "pl", cache = TRUE)

zt_data <- get_acs(geography = "zcta", 
                       variables = c(
                         male_under_5 = "B01001_003",
                         male_5_to_9 = "B01001_004",
                         male_10_to_14 = "B01001_005",
                         male_15_to_17 = "B01001_006",
                         female_under_5 = "B01001_027",
                         female_5_to_9 = "B01001_028",
                         female_10_to_14 = "B01001_029",
                         female_15_to_17 = "B01001_030"
                       ), 
                       year =2022,
                       survey = "acs5") # Specify the ACS 5-year survey

# Clean the data and sum the relevant age groups
zt_data_cleaned <- zt_data %>%
  select(-moe) %>%
  pivot_wider(names_from = variable, values_from = estimate) %>%
  mutate(total_population_under_18 = male_under_5 + male_5_to_9 + male_10_to_14 + male_15_to_17 +
           female_under_5 + female_5_to_9 + female_10_to_14 + female_15_to_17) %>% 
  select(GEOID, total_population_under_18) 

write.fst(zt_data_cleaned, paste0(pro,"/","processed_denominator.fst"))

