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
         rolling_avg_precipitation = scale(rolling_avg_precipitation)) 

# Step 1: Fit the model on data where smoke_day == 0 --------------------------------------
# Filter the data to include only observations where smoke_day == 0
dat_no_smoke <- dat %>% filter(smoke_day == 0)

# Prepare the response variable (counts of asthma events)
y_no_smoke <- dat_no_smoke$asthma

# Prepare the covariates matrix (without smoke_day)
xreg_no_smoke <- model.matrix(~ holiday + dow + ns(date, 40) +  
                                ns(rolling_avg_temperature, 5) + 
                                rolling_avg_precipitation + caserate, 
                              data = dat_no_smoke)


# first get baseline to compare:
mean_baseline <- mean(y_no_smoke)
baseline_predictions <- rep(mean_baseline, length(y_no_smoke))
mse_baseline <- mean((y_no_smoke - baseline_predictions)^2)
mae_baseline <- mean(abs(y_no_smoke - baseline_predictions))
smape_baseline <- mean(2 * abs(y_no_smoke - baseline_predictions) / (abs(y_no_smoke) + abs(baseline_predictions))) * 100

print(paste("Baseline MSE:", mse_baseline))
print(paste("Baseline MAE:", mae_baseline))
print(paste("Baseline SMAPE:", smape_baseline))



# Define hyperparameter grids for past_obs and past_mean lags
obs_lags <- list(c(1), c(1:3), c(1:5), c(1:7), c(1:10))  # Different lag combinations
mean_lags <- list(c(14, 28), c(14, 28, 90), c(14, 28, 90, 180), c(14, 28, 90, 365))

# Store results
results <- list()
best_smape <- Inf  # Track the best SMAPE score
best_mse <- Inf  # Track the best MSE score
best_mae <- Inf  # Track the best MAE score
best_smape_model <- NULL  # Store the best model based on SMAPE
best_mse_model <- NULL  # Store the best model based on MSE
best_mae_model <- NULL  # Store the best model based on MAE

for (obs in obs_lags) {
  for (mean in mean_lags) {
    # Fit the model
    model <- tsglm(
      ts = y_no_smoke,
      xreg = xreg_no_smoke,
      model = list(past_obs = obs, past_mean = mean),
      link = "log",
      distr = "nbinom"
    )
    
    
    # Get predicted values (fitted.values)
    predicted_values <- model$fitted.values
    
    # Calculate MSE
    mse <- mean((y_no_smoke - predicted_values)^2)
    
    # Calculate MAE
    mae <- mean(abs(y_no_smoke - predicted_values))
    
    # Calculate SMAPE
    smape <- mean(2 * abs(y_no_smoke - predicted_values) / (abs(y_no_smoke) + abs(predicted_values))) * 100
    
    # Store the results with a unique key
    key <- paste(paste(obs, collapse = ","), paste(mean, collapse = ","), sep = "_")
    results[[key]] <- list(model = model, MSE = mse, MAE = mae, SMAPE = smape)
    
    # Track the best MSE model
    if (mse < best_mse) {
      best_mse <- mse
      best_mse_model <- model
    }
    
    # Track the best MAE model
    if (mae < best_mae) {
      best_mae <- mae
      best_mae_model <- model
    }
    
    # Track the best SMAPE model
    if (smape < best_smape) {
      best_smape <- smape
      best_smape_model <- model
    }
    
    cat("Obs Lags:", obs, "Mean Lags:", mean, "MSE:", mse, "MAE:", mae, "SMAPE:", smape, "\n")
  }
}

# Summary of the best models
cat("Best Model MSE:", best_mse, "\n")
summary(best_mse_model)

cat("Best Model MAE:", best_mae, "\n")
summary(best_mae_model)

cat("Best Model SMAPE:", best_smape, "\n")
summary(best_smape_model)



# resp
# Obs Lags: 1 2 3 4 5 6 7 8 9 10 Mean Lags: 14 28 MSE: 64.26519 MAE: 5.930153 SMAPE: 22.51717 
 
# asthma:
# Obs Lags: 1 2 3 4 5 6 7 8 9 10 Mean Lags: 14 28 90 180 MSE: 13.54535 MAE: 2.731552 SMAPE: 41.36822 


