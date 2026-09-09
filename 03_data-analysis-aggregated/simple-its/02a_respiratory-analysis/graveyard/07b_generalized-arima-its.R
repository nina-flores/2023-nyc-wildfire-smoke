# Load necessary libraries
library(glarma)
library(ggplot2)
library(dplyr)
library(fst)

# read in data -------------------------------------------------------------------------
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
resu <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/results"


# Read in the data (modify the file path accordingly)
dat <- read.fst(paste0(an, "/", "analytical_data_prepped_nyc.fst")) %>%
  mutate(smoke_day = as.factor(smoke_day),
         rolling_avg_temperature = scale(rolling_avg_temperature),
         rolling_avg_precipitation = scale(rolling_avg_precipitation))

# Step 1: Fit the model on data where smoke_day == 0 --------------------------------------
# Filter the data to include only observations where smoke_day == 0
dat_no_smoke <- dat %>% filter(smoke_day == 0)

# Create the model matrix for the no-smoke data
X_no_smoke <- model.matrix(~ ns(date, 32)  + holiday +
                             dow + rolling_avg_temperature +
                             rolling_avg_precipitation + caserate, data = dat_no_smoke)

#initial(y = dat_no_smoke$respiratory,X = X_no_smoke,
 #       type = "NegBin")

# Fit the GARIMA model using the no-smoke data
fit_garima_no_smoke <- glarma(
  y = dat_no_smoke$respiratory,  # The count response variable (dependent variable)
  X = X_no_smoke,  # Covariates
  type = "NegBin",  # Poisson family (you can also use "NegBin" for Negative Binomial)
 # phiLags = c(1:8), 
  thetaLags = c(14, 30),
 # alphaInit = initial()
  
)


# Summary of the fitted model
summary(fit_garima_no_smoke)



# Filter the next 24 observations where smoke_day == 1
dat_smoke_24 <- dat %>% filter(smoke_day == 1) 

# Create the model matrix for the smoke-day observations (next 24 observations)
X_smoke_24 <- model.matrix(~ ns(date, 32) + holiday +
                             dow + rolling_avg_temperature +
                             rolling_avg_precipitation + caserate, data = dat_smoke_24)


# Forecast the next 24 days where smoke_day == 1 using the fitted model
forecast_results <- glarma::forecast(fit_garima_no_smoke, n.ahead = 24, 
                                     newdata = X_smoke_24)


# Extract the predicted counts (forecasted values)
predicted_forecast <- forecast_results$Forecast

# Extract the observed counts for the 24 days where smoke_day == 1
observed_counts_smoke_24 <- dat_smoke_24$respiratory

# Extract the corresponding dates (or time variable) for the 24 smoke-day observations
dates_smoke_24 <- dat_smoke_24$time_post_intervention

# Create a data frame with observed and predicted values for the smoke-day observations
plot_data <- data.frame(
  date = dates_smoke_24,
  observed = observed_counts_smoke_24,
  predicted = predicted_forecast
)

# Calculate confidence intervals from the forecast results
upper_bound <- forecast_results$Upper
lower_bound <- forecast_results$Lower

# Add confidence intervals to the plot_data data frame
plot_data$upper_bound <- upper_bound
plot_data$lower_bound <- lower_bound

# Plot observed vs. forecasted counts with confidence intervals for the smoke-day observations
ggplot(plot_data, aes(x = date)) +
  geom_line(aes(y = observed, color = "Observed"), size = 1) +  # Plot observed counts
  geom_line(aes(y = predicted, color = "Predicted"), size = 1, linetype = "dashed") +  # Plot predicted counts
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "grey") +  # Add confidence intervals
  labs(title = "Observed vs. Forecasted Counts (Smoke Days)",
       x = "Days",
       y = "Respiratory Count") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title



print(summary(fit_garima_no_smoke))



