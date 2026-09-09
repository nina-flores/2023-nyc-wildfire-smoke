require(data.table)
require(fst)
require(dplyr)
require(lubridate)


# prep nldas --------------------------------------------------------------

setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/cov")
nldas <- read.fst("nldas_daily.fst") %>%
  mutate(zcta = as.character(ZCTA5CE10)) %>%
  dplyr::select(-ZCTA5CE10) %>%
  rename("date" = "date_eastern")


# prep_pm -----------------------------------------------------------------

setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates")
pm <- read.fst("IDW_pm_zcta.fst") %>%
  mutate(zcta = as.character(ZCTA5CE10),
         date = as.Date(date, format = "%m/%d/%Y")) %>%
  dplyr::select(-ZCTA5CE10) 


# prep covid --------------------------------------------------------------
setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates")
covid <- read.fst("processed_covid.fst")


# prep denom --------------------------------------------------------------
setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates")
denom <- read.fst("processed_denominator.fst") %>%
  rename(zcta = GEOID)



# prep nevi ---------------------------------------------------------------
setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/nevi")
nevi <- read.csv("nevi_zip_final.csv") %>%
  mutate(zip_code = as.character(zip))

# transform to zcta 
# ZIP Crosswalk 
zip_xwalk <- read_xlsx(here("data", "census", "ZIPCodetoZCTACrosswalk2019.xlsx")) %>% 
  janitor::clean_names() %>% 
  dplyr::select(zip_code, zcta) 


zcta_nevi <- left_join(nevi, zip_xwalk, by="zip_code") %>% 
  dplyr::select(-zip_code) %>% 
  mutate(zcta = as.character(zcta)) %>%
  mutate(nevi_quartile = ntile(nevi,4))


# combine data ------------------------------------------------------------

covariates <- nldas %>%
  full_join(pm) %>%
  full_join(covid) %>%
  mutate(caserate = if_else(is.na(caserate), 0, caserate)) %>%
  full_join(zcta_nevi) %>%
  left_join(denom) %>%
  na.omit() %>%
  filter(date < "2023-09-01") %>%
  dplyr::select(date,
         zcta,
         total_population_under_18,
         caserate,
         temperature,
         ah,
         rh,
         total_precipitation,
         pm,
         nevi,
         nevi_quartile)

setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/processed-covariates")
write.fst(covariates, "wf-covariates.fst")


# plot pm over time -------------------------------------------------------
# Create the time series plot
ggplot(covariates, aes(x = date, y = pm)) +
  geom_line() +
  labs(title = expression("Daily PM"[2.5]), x = "Date", y = expression("PM"[2.5])) +
  theme_minimal()

covariates_cut <- covariates %>%
  filter(date > "2023-06-01" & date < "2023-06-29")

ggplot(covariates_cut, aes(x = date, y = pm)) +
  geom_line() +
  labs(title = expression("Daily PM"[2.5]), x = "Date", y = expression("PM"[2.5])) +
  theme_minimal() +
  geom_vline(xintercept = as.Date("2023-06-06"), color = "orange", linetype = "longdash")+
  geom_vline(xintercept = as.Date("2023-06-08"), color = "orange", linetype = "longdash")


