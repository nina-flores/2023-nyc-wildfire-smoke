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
source(paste0(here("scripts","outcome-processing", "0_read_libraries.R")))

# read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"


dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst")) %>%
  mutate(smoke_day = as.factor(smoke_day))


# main model 




###mod_main = glmer.nb(respiratory ~ 
###                            smoke_day + 
###                            time_elapsed + 
###                            time_post_intervention +
###                            year + 
###                            month + 
###                            dow + 
###                            holiday + 
###                            ns(temperature_scaled, df = 2) + 
###                            ns(abs_hum_scaled, df = 2) +
###                            total_precipitation_scaled + 
###                            covid +
###                            (1|zcta),
###                          data = dat,
###                          nAGQ = 0)
###

mod_main = glmer.nb(respiratory ~ 
                      smoke_day * I(time_elapsed - 2714) + # first day of smoke way was on 2714 day
                      year + 
                      month + 
                      dow + 
                      holiday + 
                      ns(temperature_scaled, df = 2) + 
                      ns(abs_hum_scaled, df = 2) +
                      total_precipitation_scaled + 
                      covid +
                      (1|zcta),
                    data = dat,
                    nAGQ = 0)

summary(mod_main)

# 3a Make data for plot
#    Notes: Set intervention to FALSE to predict the counterfactual ED visits
data_for_its_plot <- dat %>% mutate(smoke_day = "0") 

# 3b Predict ED visits for actual intervention scenario
data_for_its_plot$preds_actual = predict(mod_main)
#data_for_its_plot$preds_actual_se = predict(mod_main, se.fit = TRUE)$se.fit


# 3c Predict for counterfactual
data_for_its_plot$preds_counterf = predict(mod_main, data_for_its_plot)
#data_for_its_plot$preds_counterf_se = predict(mod_main, data_for_its_plot, se.fit = TRUE)$se.fit


# 3d Average predictions for all monitors
data_for_its_plot_summarized <- data_for_its_plot %>% group_by(time_elapsed) %>% 
  summarise(preds_actual = mean(preds_actual),
            preds_counterf = mean(preds_counterf)) %>%
  mutate(preds_actual = exp(preds_actual),
         preds_counterf = exp(preds_counterf))

# 3e Add date 
dt <- dat %>% group_by(zcta) %>% 
  dplyr::select(date, time_elapsed) %>% distinct()

# 3f Average by day and create rolling daily averages
data_for_its_plot2 <- data_for_its_plot %>% 
  left_join(dt, by = "time_elapsed", copy = TRUE) %>% 
  filter(date >= lubridate::as_date("2023-01-01"))

# 3g Create colors vector
colors <- c("Actual" = "cornflowerblue", "Counterfactual" = "orange")

# 3h Create plot
its_results_plot <- data_for_its_plot2 %>% ggplot(aes(x = date)) +

  geom_point(aes(y = preds_counterf, color = "Counterfactual"), 
             alpha = 0.8, shape = "circle open") +
  geom_point(aes(y = preds_actual, color = "Actual"),  
             alpha = 0.8) + 

  geom_line(aes(y = preds_counterf, color = "Counterfactual"), size = 1) +
  geom_line(aes(y = preds_actual, color = "Actual"), size = 1)  +
  xlab("Date") + ylab(expression("Asthma and asthma-assocatied respiratory ED visits")) + labs(color = "")  +  # Intervention point
  scale_color_manual(values = colors) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.1, 0.1))
its_results_plot


