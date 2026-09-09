# Load necessary libraries
library(tscount)
library(ggplot2)
library(dplyr)
library(fst)
library(splines)

# Read in the data (modify the file path accordingly)
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
resu <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/results"

dat <- read.fst(paste0(an, "/", "analytical_data_prepped_nyc.fst")) %>%
  mutate(smoke_day = as.factor(smoke_day),
         rolling_avg_temperature = scale(rolling_avg_temperature),
         rolling_avg_precipitation = scale(rolling_avg_precipitation)) %>%
  mutate(holiday = if_else(month == 6, 0, holiday))

# Step 1: Fit the model on data where smoke_day == 0 --------------------------------------
# Filter the data to include only observations where smoke_day == 0
dat_no_smoke <- dat %>% filter(smoke_day == 0)

# Prepare the response variable (counts of respiratory events)
y_no_smoke <- dat_no_smoke$respiratory

# Prepare the covariates matrix (without smoke_day)
xreg_no_smoke <- model.matrix(~ holiday + dow + ns(date, 40) +  
                                ns(rolling_avg_temperature, 5) + 
                                rolling_avg_precipitation + caserate, 
                              data = dat_no_smoke)


#vif_values <- car::vif(lm(respiratory ~ time_elapsed_scaled + holiday + year + month + dow + 
#                           rolling_avg_temperature + rolling_avg_precipitation + caserate , data = dat_no_smoke))
#print(vif_values)

# Fit the tsglm model with external regressors and Negative Binomial distribution
fit_tsglm_no_smoke <- tsglm(ts = y_no_smoke,
                            xreg = xreg_no_smoke,
                            model = list(past_obs = c(1:7), past_mean = c(14, 28)),  
                            link = "log",  # Log-link for count data
                            distr = "nbinom")  # Negative Binomial distribution


# Summary of the fitted model
summary(fit_tsglm_no_smoke)
#plot(fit_tsglm_no_smoke)

# see what this looks like: -----------------------------------------------
# Generate the predicted values using the fitted tsglm model
predicted_values <- fit_tsglm_no_smoke$fitted.values
# Create a data frame containing the actual observed values and the predicted values
plot_data <- data.frame(
  date = dat_no_smoke$time_elapsed_scaled,  # Use your time variable here
  observed = y_no_smoke,  # Actual respiratory counts
  predicted = predicted_values  # Predicted counts
) #%>% filter(date >= 1.5) 
# Plot the observed and predicted values
ggplot(plot_data, aes(x = date)) +
  geom_line(aes(y = observed, color = "Observed"), size = 1) +  # Plot observed counts
  geom_line(aes(y = predicted, color = "Predicted"), size = 1, linetype = "dashed") +  # Plot predicted counts
  labs(title = "Observed vs. Predicted respiratory Counts",
       x = "Time",
       y = "respiratory Counts") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove legend title


# Step 2: Prepare data for the final 24 days (where smoke_day == 1) -----------------------
# Filter the data for the last 24 days where smoke_day == 1
dat_smoke_24 <- dat %>%
  filter(smoke_day == 1) #%>%
  #filter(time_post_intervention <=15)

# Prepare the covariates matrix for the forecasted period (smoke days)
xreg_smoke_24 <- model.matrix(~ holiday + dow + ns(date, 40) +  
                                ns(rolling_avg_temperature, 5) + rolling_avg_precipitation + caserate, 
                              data = dat_smoke_24)



# Step 3: Forecast the final 24 days ------------------------------------------------------
# Forecast the final 24 days using the fitted model
#forecast_results <- predict(fit_tsglm_no_smoke, 
 #                           n.ahead = 24, 
  #                          newxreg = xreg_smoke_24,
   #                         method = "conddistr")

n_ahead <- 24  # Forecast horizon
forecasts <- numeric(n_ahead)  # Store forecasted values
lower_bounds <- numeric(n_ahead)  # Initialize lower bound array
upper_bounds <- numeric(n_ahead)  # Initialize upper bound array

for (i in 1:n_ahead) {
  # Predict 1 step ahead using the 'conddistr' method
  forecast_result <- predict(
    fit_tsglm_no_smoke,
    n.ahead = 1,
    newobs = if (i > 1) forecasts[i - 1] else NULL,
    newxreg = xreg_smoke_24[i, , drop = FALSE],
    level = 0.95,  # 95% prediction interval
    method = "conddistr"
  )
  
  # Store forecast and confidence interval values
  forecasts[i] <- forecast_result$pred
  lower_bounds[i] <- forecast_result$interval[1, "lower"]
  upper_bounds[i] <- forecast_result$interval[1, "upper"]
}

# Create a data frame for plotting
plot_data <- data.frame(
  date = dat_smoke_24$time_post_intervention,
  observed = dat_smoke_24$respiratory,
  predicted = forecasts,
  lower_bound = lower_bounds,
  upper_bound = upper_bounds
)

# Plot the results
ggplot(plot_data, aes(x = date)) +
  geom_line(aes(y = observed, color = "Observed")) +
  geom_line(aes(y = predicted, color = "Predicted"), linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "grey") +
  labs(title = "Observed vs. Forecasted Counts (Smoke Days)",
       x = "Date", y = "respiratory Counts") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())

plot(fit_tsglm_no_smoke, ask = FALSE)

