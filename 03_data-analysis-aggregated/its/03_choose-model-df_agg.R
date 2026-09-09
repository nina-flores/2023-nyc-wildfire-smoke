###########***********************
#### Code Description ####
# Author: Nina
# Goal: Test different spline df / harmonic combinations using glmmTMB(family = nbinom1())
#       models include covid_post + covid_acute; each outcome is tuned in three
#       sequential steps: precipitation df -> temperature df -> harmonic pairs,
#       each step holding the other two at their current best/default value
####**********************
####*

# Load packages
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(DHARMa)
library(splines)
library(glmmTMB)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# Read in data -------------------------------------------------------------------------
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated-agg.fst"))

# Helper: fit a glmmTMB(nbinom1()) model, catching warnings/errors, and return
# diagnostics analogous to the glm.nb $converged / $theta checks used previously -----
fit_nb1_safe <- function(formula, data, label) {
  
  mod <- tryCatch({
    glmmTMB(formula, family = nbinom1(), data = data)
  }, warning = function(w) {
    message(paste("Warning for", label, ":", conditionMessage(w)))
    suppressWarnings(glmmTMB(formula, family = nbinom1(), data = data))
  }, error = function(e) {
    message(paste("Error for", label, ":", conditionMessage(e)))
    NULL
  })
  
  if (is.null(mod)) return(NULL)
  
  conv_code <- mod$fit$convergence
  pdHess    <- isTRUE(mod$sdr$pdHess)
  converged_ok <- (!is.null(conv_code) && conv_code == 0 && pdHess)
  
  disp_val <- tryCatch(sigma(mod), error = function(e) NA_real_)
  disp_unstable <- is.na(disp_val) || disp_val > 1e4
  
  if (!converged_ok) {
    warning(paste("Model", label, "did not converge (convergence code =",
                  conv_code, ", pdHess =", pdHess, ")"))
    return(NULL)
  } else if (disp_unstable) {
    warning(paste("Model", label, "has an unstable dispersion estimate (dispersion =",
                  round(disp_val, 3), ") - inspect before trusting"))
    return(NULL)
  }
  
  list(model = mod, dispersion = disp_val, converged = converged_ok, AIC = AIC(mod))
}


# =========================================================================================
# RESPIRATORY -------------------------------------------------------------------------
# =========================================================================================

# Step 1: precipitation df, holding temperature linear and harmonics = 4 -------------------
precip_values <- c(2, 3, 4, 5)

resp_precip_models <- list()
resp_precip_aic <- data.frame(Precip_df = integer(), AIC = numeric(),
                              dispersion = numeric(), stringsAsFactors = FALSE)

for (precip_df in precip_values) {
  
  f <- as.formula(paste0(
    "respiratory ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,4,365) + dow  + ",
    "rolling_avg_temperature_scaled + ns(rolling_avg_precipitation_scaled, df = ", precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("respiratory precip df =", precip_df))
  if (is.null(fit)) next
  
  resp_precip_models[[paste0("precip_", precip_df)]] <- fit$model
  resp_precip_aic <- rbind(resp_precip_aic, data.frame(
    Precip_df = precip_df, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(resp_precip_aic)
best_resp_precip_df <- resp_precip_aic$Precip_df[which.min(resp_precip_aic$AIC)]
print(paste("Best respiratory precip df =", best_resp_precip_df))


# Step 2: temperature df, holding precip at the best df from Step 1, harmonics = 4 ---------
temp_values <- c(2, 3, 4, 5)

resp_temp_models <- list()
resp_temp_aic <- data.frame(Temperature_df = integer(), AIC = numeric(),
                            dispersion = numeric(), stringsAsFactors = FALSE)

for (temp_df in temp_values) {
  
  f <- as.formula(paste0(
    "respiratory ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,4,365) + dow + ",
    "ns(rolling_avg_temperature_scaled, df = ", temp_df, ") + ",
    "ns(rolling_avg_precipitation_scaled, df = ", best_resp_precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("respiratory temp df =", temp_df))
  if (is.null(fit)) next
  
  resp_temp_models[[paste0("temp_", temp_df)]] <- fit$model
  resp_temp_aic <- rbind(resp_temp_aic, data.frame(
    Temperature_df = temp_df, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(resp_temp_aic)
best_resp_temp_df <- resp_temp_aic$Temperature_df[which.min(resp_temp_aic$AIC)]
print(paste("Best respiratory temperature df =", best_resp_temp_df))


# Step 3: harmonic pairs, holding temp/precip df at their Step 1-2 best values -------------
harmonic_values <- 1:5

resp_harm_models <- list()
resp_harm_aic <- data.frame(Harmonics = integer(), AIC = numeric(),
                            dispersion = numeric(), stringsAsFactors = FALSE)

for (k in harmonic_values) {
  
  f <- as.formula(paste0(
    "respiratory ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,", k, ",365) + dow + ",
    "ns(rolling_avg_temperature_scaled, df = ", best_resp_temp_df, ") + ",
    "ns(rolling_avg_precipitation_scaled, df = ", best_resp_precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("respiratory harmonics =", k))
  if (is.null(fit)) next
  
  resp_harm_models[[paste0("harm_", k)]] <- fit$model
  resp_harm_aic <- rbind(resp_harm_aic, data.frame(
    Harmonics = k, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(resp_harm_aic)
best_resp_harmonics <- resp_harm_aic$Harmonics[which.min(resp_harm_aic$AIC)]
print(paste("Best respiratory harmonic pairs =", best_resp_harmonics))

best_resp_model <- resp_harm_models[[paste0("harm_", best_resp_harmonics)]]
summary(best_resp_model)


# =========================================================================================
# ASTHMA --------------------------------------------------------------------------------
# =========================================================================================

# Step 1: precipitation df, holding temperature linear and harmonics = 4 -------------------
asthma_precip_models <- list()
asthma_precip_aic <- data.frame(Precip_df = integer(), AIC = numeric(),
                                dispersion = numeric(), stringsAsFactors = FALSE)

for (precip_df in precip_values) {
  
  f <- as.formula(paste0(
    "asthma ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,4,365) + dow  + ",
    "rolling_avg_temperature_scaled + ns(rolling_avg_precipitation_scaled, df = ", precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("asthma precip df =", precip_df))
  if (is.null(fit)) next
  
  asthma_precip_models[[paste0("precip_", precip_df)]] <- fit$model
  asthma_precip_aic <- rbind(asthma_precip_aic, data.frame(
    Precip_df = precip_df, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(asthma_precip_aic)
best_asthma_precip_df <- asthma_precip_aic$Precip_df[which.min(asthma_precip_aic$AIC)]
print(paste("Best asthma precip df =", best_asthma_precip_df))


# Step 2: temperature df, holding precip at the best df from Step 1, harmonics = 4 ---------
asthma_temp_models <- list()
asthma_temp_aic <- data.frame(Temperature_df = integer(), AIC = numeric(),
                              dispersion = numeric(), stringsAsFactors = FALSE)

for (temp_df in temp_values) {
  
  f <- as.formula(paste0(
    "asthma ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,4,365) + dow + ",
    "ns(rolling_avg_temperature_scaled, df = ", temp_df, ") + ",
    "ns(rolling_avg_precipitation_scaled, df = ", best_asthma_precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("asthma temp df =", temp_df))
  if (is.null(fit)) next
  
  asthma_temp_models[[paste0("temp_", temp_df)]] <- fit$model
  asthma_temp_aic <- rbind(asthma_temp_aic, data.frame(
    Temperature_df = temp_df, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(asthma_temp_aic)
best_asthma_temp_df <- asthma_temp_aic$Temperature_df[which.min(asthma_temp_aic$AIC)]
print(paste("Best asthma temperature df =", best_asthma_temp_df))


# Step 3: harmonic pairs, holding temp/precip df at their Step 1-2 best values -------------
asthma_harm_models <- list()
asthma_harm_aic <- data.frame(Harmonics = integer(), AIC = numeric(),
                              dispersion = numeric(), stringsAsFactors = FALSE)

for (k in harmonic_values) {
  
  f <- as.formula(paste0(
    "asthma ~ smoke_day + time_elapsed_scaled + as.factor(holiday) + ",
    "harmonic(doy,", k, ",365) + dow + ",
    "ns(rolling_avg_temperature_scaled, df = ", best_asthma_temp_df, ") + ",
    "ns(rolling_avg_precipitation_scaled, df = ", best_asthma_precip_df, ") + ",
    "covid_post + covid_acute"
  ))
  
  fit <- fit_nb1_safe(f, dat, paste("asthma harmonics =", k))
  if (is.null(fit)) next
  
  asthma_harm_models[[paste0("harm_", k)]] <- fit$model
  asthma_harm_aic <- rbind(asthma_harm_aic, data.frame(
    Harmonics = k, AIC = fit$AIC, dispersion = fit$dispersion
  ))
}

print(asthma_harm_aic)
best_asthma_harmonics <- asthma_harm_aic$Harmonics[which.min(asthma_harm_aic$AIC)]
print(paste("Best asthma harmonic pairs =", best_asthma_harmonics))

best_asthma_model <- asthma_harm_models[[paste0("harm_", best_asthma_harmonics)]]
summary(best_asthma_model)

# Resp harmonics: 2 or 4
# Asthma harmonics: 2 

#Resp precip: linear
#Asthma precip: linear

#Resp temp: 5
#Asthma temp: 5
