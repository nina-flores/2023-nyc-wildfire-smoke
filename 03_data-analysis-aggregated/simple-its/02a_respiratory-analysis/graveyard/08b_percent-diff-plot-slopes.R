# Load necessary libraries
library(here)
library(tigris)
library(datetimeutils)
library(lme4)
library(DHARMa)
library(splines)
library(dplyr)
library(ggplot2)
library(lubridate)

# Source custom libraries (if needed)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# Load data
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an, "/", "analytical_data_prepped.fst"))

# Fit main model
mod_main = glmer.nb(respiratory ~ smoke_day + time_post_intervention + ns(time_elapsed_scaled, 15) + holiday +
                      ns(rolling_avg_temperature_scaled, 5) + rolling_avg_precipitation_scaled +
                      (1 + smoke_day + time_post_intervention | zcta),
                    data = dat,
                    nAGQ = 0)

summary(mod_main)
car::vif(mod_main)
# Extract random effects for ZCTA
random_effects <- ranef(mod_main)$zcta
random_effects_df <- as.data.frame(random_effects)

# Get fixed effects for smoke_day and time_post_intervention
fixed_effects <- fixef(mod_main)
fixed_smoke_day <- fixed_effects["smoke_day1"]  # assuming smoke_day is a factor
fixed_time_post_intervention <- fixed_effects["time_post_intervention"]

# Combine fixed effects with random effects to get ZCTA-specific slopes
random_effects_df$zcta_slope_smoke_day <- fixed_smoke_day + random_effects_df$smoke_day1
random_effects_df$zcta_slope_time_post_intervention <- fixed_time_post_intervention + random_effects_df$time_post_intervention

random_effects_df <- random_effects_df %>%
  rownames_to_column(var = "zcta") %>% # Move ZCTA from row names to a column
  select(zcta, zcta_slope_smoke_day, zcta_slope_time_post_intervention)

# Create trajectory data for plotting
dat_trajectory <- dat %>%
  select(zcta, time_elapsed_scaled, smoke_day, time_post_intervention, nevi_quartile) %>%
  distinct() %>%
  left_join(random_effects_df, by = "zcta") %>%
  group_by(zcta, time_post_intervention) %>%
  mutate(predicted_trajectory_smoke = zcta_slope_smoke_day * as.numeric(smoke_day),
         predicted_trajectory_time_post_intervention = zcta_slope_time_post_intervention * time_post_intervention,
         predicted_its = exp(predicted_trajectory_smoke + predicted_trajectory_time_post_intervention)) %>%
  filter(time_post_intervention > 0)

# Create counterfactual predictions (smoke_day = 0)
dat_trajectory_counterfactual <- dat %>%
  select(zcta, time_elapsed_scaled, smoke_day, time_post_intervention, nevi_quartile) %>%
  distinct() %>%
  left_join(random_effects_df, by = "zcta") %>%
  filter(time_post_intervention > 0) %>%
  mutate(smoke_day = 0,  # Set smoke_day to 0 for the counterfactual scenario
         time_post_intervention2 = 0,
         predicted_trajectory_smoke_counterf = zcta_slope_smoke_day * as.numeric(smoke_day),
         predicted_trajectory_time_post_intervention_counterf = zcta_slope_time_post_intervention * time_post_intervention2,
         predicted_its_counterfactual = exp(predicted_trajectory_smoke_counterf + predicted_trajectory_time_post_intervention_counterf))

# Combine actual and counterfactual predictions
dat_trajectory <- dat_trajectory %>%
  left_join(dat_trajectory_counterfactual %>% select(zcta, time_post_intervention, predicted_its_counterfactual), 
            by = c("zcta", "time_post_intervention"))

# Calculate percent difference
dat_trajectory <- dat_trajectory %>%
  mutate(percent_diff = ((predicted_its - predicted_its_counterfactual) / predicted_its_counterfactual) * 100)

# Plot percent difference for each ZCTA over time_post_intervention
ggplot(dat_trajectory, aes(x = time_post_intervention)) +
  geom_line(aes(y = percent_diff, color = zcta)) +
  facet_wrap(~nevi_quartile) +
  labs(x = "Time Post Intervention", 
       y = "Percent Difference (%)", 
       title = "ZCTA-specific Percent Difference Between Actual and Counterfactual Values") +
  theme_minimal() 


