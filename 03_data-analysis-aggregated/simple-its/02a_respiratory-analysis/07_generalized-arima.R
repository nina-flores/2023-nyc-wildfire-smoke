library(glarma)


# read in data -------------------------------------------------------------------------
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
resu <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/results"

dat <- read.fst(paste0(an, "/", "analytical_data_prepped_nyc.fst"))  %>%
  mutate(smoke_day = as.factor(smoke_day),
         rolling_avg_temperature = scale(rolling_avg_temperature),
         rolling_avg_precipitation = scale( rolling_avg_precipitation)) %>%
  filter(smoke_day == 0)

X = model.matrix(~  ns(date, 32)  + holiday +
                    dow + rolling_avg_temperature +
                   + rolling_avg_precipitation + caserate, data = dat)


# Fit the GARIMA model using Poisson distribution with ARMA(1,1) structure
fit_garima <- glarma(
  y = dat$respiratory,  # The count response variable (dependent variable)
  X = X,  # Covariates
  type = "NegBin",  # Poisson family (you can also use "NegBin" for Negative Binomial)
  phiLags = c(1:8), 
  thetaLags = c(14, 30),
  alphaInit = initial()

)


# Summary of the fitted model
summary(fit_garima)

plot(fit_garima)

# Extract fitted values
fitted_values <- fitted(fit_garima)

# Extract residuals
residuals_values <- residuals(fit_garima)

# Plot the autocorrelation function of residuals
#acf(residuals_values, main = "ACF of Residuals")



# Number of days to predict
n_predict_days <- 365

# Extract the predicted counts from the fitted model
predicted_counts <- fitted(fit_garima)

# Subset the predicted counts for the final 24 days
predicted_counts_last_24 <- tail(predicted_counts, n_predict_days)

# Subset the observed counts (respiratory) for the final 24 days
observed_counts_last_24 <- tail(dat$respiratory, n_predict_days)

# Subset the corresponding dates (or time variable) for the final 24 days
dates_last_24 <- tail(dat$time_elapsed, n_predict_days)

# Create a data frame with observed and predicted values for the last 24 days
plot_data <- data.frame(
  date = dates_last_24,
  observed = observed_counts_last_24,
  predicted = predicted_counts_last_24
)

# Calculate confidence intervals for the predicted counts (based on residuals)
upper_bound <- predicted_counts_last_24 + 1.96 * sqrt(predicted_counts_last_24)
lower_bound <- predicted_counts_last_24 - 1.96 * sqrt(predicted_counts_last_24)

# Add confidence intervals to the plot_data data frame
plot_data$upper_bound <- upper_bound
plot_data$lower_bound <- lower_bound

# Plot observed vs. predicted counts with confidence intervals for the final 24 days
ggplot(plot_data, aes(x = date)) +
  geom_line(aes(y = observed, color = "Observed"), size = 1) +  # Plot observed counts
  geom_line(aes(y = predicted, color = "Predicted"), size = 1, linetype = "dashed") +  # Plot predicted counts
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "grey") +  # Add confidence intervals
  labs(title = "Observed vs. Predicted Counts",
       x = "Days following ",
       y = "Respiratory Count") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title







