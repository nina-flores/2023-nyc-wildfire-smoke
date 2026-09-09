###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: create its plot for final model
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

# main model 

mod_main = glmer.nb(respiratory ~ smoke_day + time_post_intervention + time_elapsed_scaled  + holiday +
                      year + month + dow + rolling_avg_temperature_scaled + rolling_avg_precipitation_scaled +
                      (1 | zcta),
                    data = dat,
                    nAGQ = 0)

summary(mod_main)

# 3a Make data for plot
#    Notes: Set intervention to FALSE to predict the counterfactual ED visits
data_for_its_plot <- dat %>% mutate(smoke_day = as.factor("0")) %>%
  mutate(time_post_intervention = 0)

# 3b Predict ED visits for actual intervention scenario
data_for_its_plot$preds_actual = predict(mod_main, type = "link")
data_for_its_plot$preds_actual_se = predict(mod_main, type = "link", se.fit = TRUE)$se.fit


# 3c Predict for counterfactual
data_for_its_plot$preds_counterf = predict(mod_main, data_for_its_plot, type = "link")
data_for_its_plot$preds_counterf_se = predict(mod_main, data_for_its_plot,type = "link", se.fit = TRUE)$se.fit


# 3d Average predictions for all monitors and calculate percent difference
data_for_its_plot <- data_for_its_plot %>%
  group_by(time_elapsed) %>% 
  summarise(preds_actual = mean(preds_actual),
            preds_counterf = mean(preds_counterf),
            preds_counterf_se = sqrt(sum(preds_counterf_se^2) / n()),
            preds_actual_se = sqrt(sum(preds_actual_se^2) / n()),
            preds_counterf_lower = preds_counterf - 1.96 * preds_counterf_se,
            preds_counterf_upper = preds_counterf + 1.96 * preds_counterf_se,
            preds_actual_lower = preds_actual - 1.96 * preds_actual_se,
            preds_actual_upper = preds_actual + 1.96 * preds_actual_se) %>%
  mutate(preds_actual = exp(preds_actual),
         preds_counterf = exp(preds_counterf),
         preds_actual_lower = exp(preds_actual_lower),
         preds_counterf_lower = exp(preds_counterf_lower),
         preds_actual_upper = exp(preds_actual_upper),
         preds_counterf_upper = exp(preds_counterf_upper),
         # Calculate the percent difference
         percent_diff = ((preds_actual - preds_counterf) / preds_counterf) * 100,
         percent_diff_lower = ((preds_actual_lower - preds_counterf_lower) / preds_counterf_lower) * 100,
         percent_diff_upper = ((preds_actual_upper - preds_counterf_upper) / preds_counterf_upper) * 100)

# 3e Add date 
dt <- dat %>% group_by(zcta) %>% 
  dplyr::select(date, time_elapsed) %>% distinct()

# 3f Average by day and create rolling daily averages
data_for_its_plot2 <- data_for_its_plot %>%
  left_join(dt, by = "time_elapsed", copy = TRUE) %>%
  filter(date >= lubridate::as_date("2023-01-01"))

# Filtered data for confidence intervals
intervention_date <- lubridate::as_date("2023-06-06")
data_for_ci <- data_for_its_plot2 %>% filter(date >= intervention_date)


# Update the plot to show percent difference
its_results_plot <- data_for_its_plot2 %>%
  ggplot(aes(x = date)) +
  geom_rect(aes(xmin = lubridate::as_date("2023-06-06"), xmax = lubridate::as_date("2023-06-08"), ymin = -Inf, ymax = Inf), fill = "grey90", alpha = 0.1) +
  geom_point(aes(y = percent_diff), alpha = 0.8, shape = "circle open", color = "blue") +
  geom_ribbon(data = data_for_ci, aes(ymin = percent_diff_lower, ymax = percent_diff_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = percent_diff), size = 1, color = "blue") +
  xlab("Day (June 2023)") +
  ylab("Percent Difference (%) between Actual and Counterfactual") +
  theme_bw(base_size = 16) +
  scale_x_date(limits = c(lubridate::as_date("2023-06-01"), lubridate::as_date("2023-06-30"))) +  # Set x-axis limits to June
  theme(
    legend.position = c(0.15, 0.15),  # Position of the legend
    legend.box = "rect",              # Add box around the legend
    legend.background = element_rect(color = "black", fill = "white", size = 0.5),  # Customize the legend box
    legend.title = element_blank()    # Remove the legend title
  ) 

# Print the plot
print(its_results_plot)
