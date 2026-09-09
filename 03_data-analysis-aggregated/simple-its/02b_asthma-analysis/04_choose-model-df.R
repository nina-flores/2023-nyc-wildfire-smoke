###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/12/24
# Goal: Test different spline df combinations on model with random effects
####**********************
####*
####*For some reason, using the rolling average of temperature and precip
####* rather than daily, leads to singular fits when assessing non-linearity.
####* Will leave these linear. 

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(lme4)
library(DHARMa)
library(splines)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# Read in data -------------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an ,"/","analytical_data_prepped.fst"))

# List of df values for testing 
df_values <- c(1, 2, 3, 4, 5)

# Empty list to store models and AIC values
models <- list()
aic_values <- data.frame(Temperature_df = integer(), Precipitation_df = integer(), AIC = numeric(), stringsAsFactors = FALSE)

# Loop through combinations of df for temperature and precipitation
for (temp_df in df_values) {
  for (precip_df in df_values) {
    
    # Model with ns() applied for both temperature and precipitation
    mod <- glmer(asthma ~ smoke_day + time_post_intervention + time_elapsed_scaled + holiday + 
                   year + month + dow + ns(rolling_avg_temperature_scaled, df = temp_df) + 
                   ns(rolling_avg_precipitation_scaled, df = precip_df) +
                   (1 | zcta),
                 data = dat,
                 family = "poisson",
                 offset = log(total_population_under_18), 
                 nAGQ = 0)
    
    # Store the model and its AIC value
    models[[paste0("mod_temp_", temp_df, "_precip_", precip_df)]] <- mod
    aic_values <- rbind(aic_values, data.frame(Temperature_df = temp_df, Precipitation_df = precip_df, AIC = AIC(mod)))
  }
}

# Display AIC values for all models
print(aic_values)

# Get the model with the lowest AIC
best_model <- aic_values[which.min(aic_values$AIC), ]
print(paste("Best model is with temperature df =", best_model$Temperature_df, 
            "and precipitation df =", best_model$Precipitation_df))

summary(best_model_object)

# Predictions for the best model --------------------------------------------------------

# Create a new data frame for prediction
new_data <- data.frame(
  rolling_avg_temperature_scaled = seq(min(dat$rolling_avg_temperature_scaled), max(dat$rolling_avg_temperature_scaled), length.out = 100),
  year = factor(names(which.min(table(dat$year))), levels = levels(dat$year)),
  zcta = factor(names(which.min(table(dat$zcta))), levels = levels(dat$zcta)),
  month = factor(names(which.min(table(dat$month))), levels = levels(dat$month)),
  dow = factor(names(which.max(table(dat$dow))), levels = levels(dat$dow)),
  holiday = 0,
  rolling_avg_precipitation_scaled = mean(dat$rolling_avg_precipitation_scaled),
  time_elapsed_scaled = mean(dat$time_elapsed_scaled),
  smoke_day = as.character(1),
  time_post_intervention = mean(dat$time_post_intervention)
)



# Extract best model object from the list
best_model_object <- models[[paste0("mod_temp_", best_model$Temperature_df, "_precip_", best_model$Precipitation_df)]]

# Predict using the best model
new_data$pred <- predict(best_model_object, newdata = new_data, type = "response")

# Generate predicted values using the best model

# Plot predictions
ggplot(new_data, aes(rolling_avg_temperature_scaled)) + 
  geom_line(aes(y = pred), color = "blue", size = 1) +
  xlab("Rolling Avg Temperature (Scaled)") +
  ylab("Predicted asthma Visits")






# Create a new data frame for prediction
new_data <- data.frame(
  rolling_avg_precipitation_scaled = seq(min(dat$rolling_avg_precipitation_scaled), max(dat$rolling_avg_precipitation_scaled), length.out = 100),
  year = factor(names(which.min(table(dat$year))), levels = levels(dat$year)),
  zcta = factor(names(which.min(table(dat$zcta))), levels = levels(dat$zcta)),
  month = factor(names(which.min(table(dat$month))), levels = levels(dat$month)),
  dow = factor(names(which.max(table(dat$dow))), levels = levels(dat$dow)),
  holiday = 0,
  rolling_avg_temperature_scaled = mean(dat$rolling_avg_temperature_scaled),
  time_elapsed_scaled = mean(dat$time_elapsed_scaled),
  smoke_day = as.character(1),
  time_post_intervention = mean(dat$time_post_intervention)
)

# Predict using the best model
new_data$pred <- predict(best_model_object, newdata = new_data, type = "response")

# Generate predicted values using the best model

# Plot predictions
ggplot(new_data, aes(rolling_avg_precipitation_scaled)) + 
  geom_line(aes(y = pred), color = "blue", size = 1) +
  xlab("Rolling Avg Precip (Scaled)") +
  ylab("Predicted asthma Visits")



