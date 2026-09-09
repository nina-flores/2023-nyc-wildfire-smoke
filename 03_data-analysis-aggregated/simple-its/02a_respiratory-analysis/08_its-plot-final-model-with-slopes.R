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
mod_main = glmer.nb(respiratory ~ time_elapsed_scaled + holiday +
                   year + month + dow + ns(rolling_avg_temperature_scaled, 5) + rolling_avg_precipitation_scaled +
                   (1 + smoke_day + time_post_intervention | zcta),
                 data = dat,
                 nAGQ = 0)


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
  select(zcta, time_elapsed_scaled, smoke_day, time_post_intervention) %>%
  distinct() %>%
  left_join(random_effects_df, by = "zcta") %>%
  mutate(predicted_trajectory_smoke = zcta_slope_smoke_day * as.numeric(smoke_day),
         predicted_trajectory_time_post_intervention = zcta_slope_time_post_intervention * time_post_intervention,
         predicted_its = exp(predicted_trajectory_smoke + predicted_trajectory_time_post_intervention)) %>%
  filter(time_post_intervention > 0)

# Plot ZCTA-specific trajectories for both smoke_day and time_post_intervention
ggplot(dat_trajectory, aes(x = time_post_intervention)) +
  geom_line(aes(y = predicted_its, color = zcta)) +
  #facet_wrap(~ zcta) +
  labs(x = "Time", y = "Predicted Values", title = "ZCTA-specific Trajectories for Smoke Day and Time Post Intervention") +
 # scale_color_manual(values = c("ITS slope" = "#a6bddb")) +
  theme_minimal()
