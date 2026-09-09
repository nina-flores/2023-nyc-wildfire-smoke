#### Wildfire smoke ITS analysis - citywide respiratory/asthma visits
#### Three approaches, run in the order:
####   1. Before/arter IRR 
####      the commonly used/published comparison - like Chen paper 
####   2. NB - excess-visit calculation, secondary
####      comparison
####   3. Two-stage ARIMAX ITS (auto.arima, fit on PRE-EVENT data only,
####      forecast the event window as the counterfactual) - excess-visit
####      calculation, primary approach; skips refitting if cached output
####      already exists
#### Plotting has been moved to a separate script.
####***********************

rm(list = ls())


# Packages ----------------------------------------------------------------

library(here)
library(dplyr)
library(purrr)
library(splines)
library(fastDummies)
library(glmmTMB)
library(forecast)
library(fst)
library(forcats)
library(tsModel)
library(DHARMa)


source("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R")


# Paths -------------------------------------------------------------------

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"

out_dir <- file.path(an, "models")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}



# Event dates -------------------------------------------------------------

event_start <- as.Date("2023-06-06")
event_end   <- as.Date("2023-06-08")


event_dates <- data.frame(
  day = c("6/6", "6/7", "6/8"),
  date = seq(event_start, event_end, by = "day")
)


# Reference period for the published-design IRR ---------------------------
# Matches the design used in the published JAMA research letter: two
# adjacent, day-of-week-matched nonsmoke windows flanking the event.

ref_before_start <- as.Date("2023-05-30")
ref_before_end   <- as.Date("2023-06-01")

ref_after_start  <- as.Date("2023-06-13")
ref_after_end    <- as.Date("2023-06-15")



# ============================================================================
# 1. Read and prep data
# ============================================================================

dat <- read.fst(
  file.path(an, "analytical_data_prepped-updated-agg.fst")
) %>%
  arrange(date) %>%
  mutate(
    date_num = as.integer(date - min(date)) + 1
  )




# ============================================================================
# 2. Published-design IRR (JAMA-style) - RR result only, run first
# ============================================================================
# Event-window rate vs pooled before+after reference-period rate.
# CI uses the standard two-sample Poisson/rate-ratio log-scale SE:
#   SE(ln IRR) = sqrt(1/visits_event + 1/visits_ref)
# This matches fmsb::riskratio() when population-time denominators are
# large relative to event counts (true here), and reproduces the
# published NYC wildfire-smoke asthma IRR (1.44, 95% CI 1.31-1.58) when
# run on the same event/reference windows.

# Generalized IRR function - takes explicit event/reference date VECTORS
# rather than start/end ranges, so the same function covers both the
# pooled "Overall" comparison (event window vs pooled before+after window)
# and each day-specific comparison (single event day vs its two
# day-of-week-matched reference days, one from the before window and one
# from the after window).

compute_published_irr <- function(data,
                                  outcome_col,
                                  event_dates_vec,
                                  ref_dates_vec) {
  
  n_event_days <- length(event_dates_vec)
  n_ref_days   <- length(ref_dates_vec)
  
  visits_event <- data %>%
    filter(date %in% event_dates_vec) %>%
    summarise(total = sum(.data[[outcome_col]])) %>%
    pull(total)
  
  visits_ref <- data %>%
    filter(date %in% ref_dates_vec) %>%
    summarise(total = sum(.data[[outcome_col]])) %>%
    pull(total)
  
  rate_event <- visits_event / n_event_days
  rate_ref   <- visits_ref / n_ref_days
  
  IRR <- rate_event / rate_ref
  
  se_log  <- sqrt(1 / visits_event + 1 / visits_ref)
  log_irr <- log(IRR)
  
  conf.low  <- exp(log_irr - 1.96 * se_log)
  conf.high <- exp(log_irr + 1.96 * se_log)
  
  data.frame(
    term = "published_design_irr",
    estimate = IRR,
    conf.low = conf.low,
    conf.high = conf.high
  )
  
}


# Overall event/reference windows (pooled, as before) ----------------------

event_dates_overall <- seq(event_start, event_end, by = "day")

ref_dates_overall <- c(
  seq(ref_before_start, ref_before_end, by = "day"),
  seq(ref_after_start,  ref_after_end,  by = "day")
)


# Day-specific matches -------------------------------------------------------
# Each event day is matched to the before-window and after-window days
# exactly 7 days apart (same day of week), consistent with the
# day-of-week-matched reference windows used for the overall comparison.

day_match <- data.frame(
  day             = event_dates$day,
  event_date      = event_dates$date,
  ref_before_date = seq(ref_before_start, ref_before_end, by = "day"),
  ref_after_date  = seq(ref_after_start,  ref_after_end,  by = "day")
)


get_irr_for_outcome <- function(outcome_col, outcome_label) {
  
  overall <- compute_published_irr(
    dat, outcome_col,
    event_dates_vec = event_dates_overall,
    ref_dates_vec   = ref_dates_overall
  ) %>%
    mutate(day = "Overall")
  
  by_day <- map_dfr(seq_len(nrow(day_match)), function(i) {
    compute_published_irr(
      dat, outcome_col,
      event_dates_vec = day_match$event_date[i],
      ref_dates_vec   = c(day_match$ref_before_date[i], day_match$ref_after_date[i])
    ) %>%
      mutate(day = day_match$day[i])
  })
  
  bind_rows(overall, by_day) %>%
    mutate(outcome = outcome_label, model = "Published design")
  
}


rr_published <- bind_rows(
  get_irr_for_outcome("respiratory", "Respiratory"),
  get_irr_for_outcome("asthma",      "Asthma")
)


print(rr_published)

saveRDS(rr_published, file.path(out_dir, "rr_published_design.rds"))
write.csv(rr_published, file.path(out_dir, "rr_published_design.csv"), row.names = FALSE)




# ============================================================================
# 3. NB models (nbinom1, no AR1) - secondary comparison, excess-visit
#    calculation only (no RR extraction for this method)
# ============================================================================

fit_nb <- function(formula, data) {

  glmmTMB(
    formula,
    family = nbinom1(),
    data = data
  )
  
}



# Single smoke-day models -----------------------------------------------


mod_respiratory <- fit_nb(
  
  respiratory ~ smoke_day +
    covid_post + 
    covid_acute +
    time_elapsed_scaled +
    as.factor(holiday) +
    harmonic(doy, 2, 365) +
    dow +
    ns(rolling_avg_temperature_scaled,5) +
    rolling_avg_precipitation_scaled ,
  
  dat
  
)


mod_asthma <- fit_nb(
  
  asthma ~ smoke_day +
    covid_post +
    covid_acute +
    time_elapsed_scaled +
    as.factor(holiday) +
    harmonic(doy,2,365) +
    dow +
    ns(rolling_avg_temperature_scaled,5) +
    rolling_avg_precipitation_scaled,
  
  dat
  
)




# Lagged smoke-day models -----------------------------------------------


mod_respiratory_multi <- fit_nb(
  
  respiratory ~ smoke_day1 +
    smoke_day2 +
    smoke_day3 +
    covid_post +
    covid_acute +
    time_elapsed_scaled +
    as.factor(holiday) +
    harmonic(doy,2,365) +
    dow +
    ns(rolling_avg_temperature_scaled,5) +
    rolling_avg_precipitation_scaled ,
  
  dat
  
)



mod_asthma_multi <- fit_nb(
  
  asthma ~ smoke_day1 +
    smoke_day2 +
    smoke_day3 +
    covid_post +
    covid_acute +
    time_elapsed_scaled +
    as.factor(holiday) +
    harmonic(doy,2,365) +
    dow +
    ns(rolling_avg_temperature_scaled,5) +
    rolling_avg_precipitation_scaled ,
  
  dat
  
)


summary(mod_respiratory)
summary(mod_respiratory_multi)
summary(mod_asthma)
summary(mod_asthma_multi)




# ============================================================================
# 4. Check model coefficients
# ============================================================================

names(fixef(mod_respiratory)$cond)
names(fixef(mod_asthma)$cond)
names(fixef(mod_respiratory_multi)$cond)
names(fixef(mod_asthma_multi)$cond)


check_smoke_effect <- function(model, smoke_var) {
  
  
  coefs <- fixef(model)$cond
  
  if (!smoke_var %in% names(coefs)) {
    
    warning(
      paste(
        "Variable",
        smoke_var,
        "not found in model coefficients"
      )
    )
    
    return(
      data.frame(
        variable = smoke_var,
        beta = NA_real_,
        RR = NA_real_,
        percent_change = NA_real_
      )
    )
    
  }
  
  
  
  beta <- unname(coefs[smoke_var])
  
  data.frame(
    
    variable = smoke_var,
    beta = beta,
    RR = exp(beta),
    percent_change = (exp(beta)-1)*100
    
  )
  
}




coef_checks <- bind_rows(
  
  check_smoke_effect(
    mod_respiratory,
    "smoke_day"
  ),
  
  check_smoke_effect(
    mod_asthma,
    "smoke_day"
  ),
  
  check_smoke_effect(
    mod_respiratory_multi,
    "smoke_day1"
  ),
  
  check_smoke_effect(
    mod_asthma_multi,
    "smoke_day1"
  )
  
)



print(coef_checks)




# ============================================================================
# 5. Calculate excess visits - NB models
# ============================================================================


get_excess_single_nb <- function(model, data, outcome_col) {
  
  # Counterfactual = no smoke exposure
  data_cf <- data %>%
    mutate(
      smoke_day = 0
    )
  
  
  
  pred.cf <- predict(
    
    model,
    newdata = data_cf,
    type = "link",
    se.fit = TRUE,
    re.form = NA
    
  )
  
  
  data_cf$pred.cf <- exp(pred.cf$fit)
  
  data_cf$pred.cf_lower <-
    exp(pred.cf$fit - 1.96 * pred.cf$se.fit)
  
  data_cf$pred.cf_upper <-
    exp(pred.cf$fit + 1.96 * pred.cf$se.fit)
  
  
  data.event <- data_cf %>%
    
    filter(
      date >= event_start &
        date <= event_end
    )
  
  
  observed <- sum(data.event[[outcome_col]])
  
  counterfactual <- sum(data.event$pred.cf)
  
  
  data.frame(
    
    day = "Overall",
    observed = observed,
    counterfactual = counterfactual,
    estimate = observed - counterfactual,
    conf.low = observed - sum(data.event$pred.cf_upper),
    conf.high = observed - sum(data.event$pred.cf_lower)
    
  )
  
}




get_excess_multi_nb <- function(model, data, outcome_col) {
  
  
  # Counterfactual = all smoke exposure removed
  
  data_cf <- data %>%
    
    mutate(
      smoke_day1 = 0,
      smoke_day2 = 0,
      smoke_day3 = 0
      
    )
  
  
  
  pred.cf <- predict(
    
    model,
    newdata = data_cf,
    type = "link",
    se.fit = TRUE,
    re.form = NA
    
  )
  
  
  data_cf$pred.cf <- exp(pred.cf$fit)
  
  
  data_cf$pred.cf_lower <- exp(pred.cf$fit - 1.96 * pred.cf$se.fit)
  
  
  data_cf$pred.cf_upper <- exp(pred.cf$fit + 1.96 * pred.cf$se.fit)
  
  
  
  
  map_dfr(seq_len(nrow(event_dates)), function(i){
    
    
    
    d <- data_cf %>%
      
      filter(
        date == event_dates$date[i]
      )
    
    
    
    observed <- d[[outcome_col]]
    
    counterfactual <- d$pred.cf
    
    
    
    data.frame(
      
      day = event_dates$day[i],
      observed = observed,
      counterfactual = counterfactual,
      estimate = observed - counterfactual,
      conf.low = observed - d$pred.cf_upper,
      conf.high = observed - d$pred.cf_lower
      
    )
    
    
  })
  
}




# ============================================================================
# 6. Generate NB excess estimates
# ============================================================================


excess_respiratory_nb <-
  
  get_excess_single_nb(
    mod_respiratory,
    dat,
    "respiratory"
  ) %>%
  mutate(
    outcome = "Respiratory",
    model = "Single smoke day",
    source = "NB (nbinom1)"
  )



excess_respiratory_multi_nb <-
  
  get_excess_multi_nb(
    mod_respiratory_multi,
    dat,
    "respiratory"
  ) %>%
  mutate(
    outcome = "Respiratory",
    model = "Lagged smoke days",
    source = "NB (nbinom1)"
  )



excess_asthma_nb <-
  
  get_excess_single_nb(
    mod_asthma,
    dat,
    "asthma"
  ) %>%
  mutate(
    outcome = "Asthma",
    model = "Single smoke day",
    source = "NB (nbinom1)"
  )



excess_asthma_multi_nb <-
  
  get_excess_multi_nb(
    mod_asthma_multi,
    dat,
    "asthma"
  ) %>%
  mutate(
    outcome = "Asthma",
    model = "Lagged smoke days",
    source = "NB (nbinom1)"
  )




# ============================================================================
# 7. Final NB excess table
# ============================================================================


excess_results_nb <- bind_rows(
  excess_respiratory_nb,
  excess_respiratory_multi_nb,
  excess_asthma_nb,
  excess_asthma_multi_nb
)



print(excess_results_nb)



excess_results_nb %>%
  select(
    outcome,
    model,
    day,
    observed,
    counterfactual,
    estimate
  ) %>%  print()




# ============================================================================
# 8. NB diagnostic panels (DHARMa residuals - fit quality check)
# ============================================================================
# NOTE: diagnostic plots use ggplot/patchwork, which are otherwise no
# longer loaded above since final plotting was moved out of this script.
# Load them locally here so diagnostics still work standalone.

library(ggplot2)
library(patchwork)

make_diagnostic_panel <- function(model, outcome_label) {
  
  sim_res <- simulateResiduals(model, n = 1000)
  
  resid_df <- data.frame(
    fitted = predict(model, type = "response"),
    resid  = residuals(sim_res, quantileFunction = qnorm)
  )
  
  p1 <- ggplot(resid_df, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.4, size = 1, color = "grey30") +
    geom_hline(yintercept = 0, color = "firebrick", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "steelblue", linewidth = 0.7) +
    labs(x = "Fitted values", y = "Residuals",
         title = "A. Residuals vs. Fitted") +
    theme_bw(base_size = 12)
  
  acf_obj <- acf(resid_df$resid, plot = FALSE, lag.max = 30)
  acf_df <- data.frame(lag = acf_obj$lag[,1,1], acf = acf_obj$acf[,1,1])
  ci_line <- qnorm(0.975) / sqrt(acf_obj$n.used)
  
  p2 <- ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-ci_line, ci_line), color = "firebrick", linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0), color = "grey30") +
    labs(x = "Lag (days)", y = "ACF",
         title = "B. Autocorrelation of Residuals") +
    theme_bw(base_size = 12)
  
  p3 <- ggplot(resid_df, aes(sample = resid)) +
    stat_qq(alpha = 0.4, size = 1, color = "grey30") +
    stat_qq_line(color = "firebrick", linetype = "dashed") +
    labs(x = "Theoretical quantiles", y = "Sample quantiles",
         title = "C. Normal Q-Q Plot") +
    theme_bw(base_size = 12)
  
  (p1 | p2 | p3) +
    plot_annotation(title = outcome_label,
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
}

diag_respiratory <- make_diagnostic_panel(mod_respiratory, "Respiratory Visits (NB, nbinom1, no AR1)")
diag_respiratory

diag_asthma <- make_diagnostic_panel(mod_asthma, "Asthma Visits (NB, nbinom1, no AR1)")
diag_asthma




# ============================================================================
# 9. Two-stage ARIMAX ITS - fit on PRE-EVENT data, forecast event window
#    NOTE: smoke_day is NOT in xreg at all - forecasting into the event
#    period using pre-event dynamics IS the counterfactual. No zeroing-out
#    needed, matching the classical SF county ITS script exactly.
#    Skips refitting (auto.arima) if cached .rds output already exists for
#    that outcome (slow to fit) - but the xreg/design-matrix construction
#    below always runs regardless of caching, since xreg_event_* is needed
#    for the section 10 bootstrap simulation even when the model itself is
#    loaded from cache.
#
#    Harmonic count matches the NB models per outcome: 8 harmonic pairs
#    for respiratory, 10 for asthma. Two separate xreg matrices are built
#    (respiratory / asthma) instead of one shared xreg, since the two
#    outcomes use different numbers of Fourier terms.
# ============================================================================

arimax_respiratory_path <- file.path(out_dir, "arimax_respiratory.rds")
arimax_asthma_path      <- file.path(out_dir, "arimax_asthma.rds")

# Number of bootstrap paths used to build the Overall excess CI (section 10) -
# simulates future paths from the fitted ARIMA model's residual distribution
# rather than summing the daily forecast interval bounds, so day-to-day
# forecast-error correlation is properly propagated into the multi-day total.
ARIMAX_SIM_NPATHS <- 5000
ARIMAX_SIM_SEED   <- 1

dat_city <- dat %>%
  # drop leap days so frequency = 365 stays consistent across leap years
  # (2016, 2020) in the training window - ts() with frequency = 365
  # assumes exactly 365 obs/year, and leaving Feb 29 in shifts doy-based
  # harmonics out of phase with the ts index past that point
  filter(!(format(date, "%m") == "02" & format(date, "%d") == "29")) %>%
  mutate(
    cos1 = cos(2*pi*1*doy/365), sin1 = sin(2*pi*1*doy/365),
    cos2 = cos(2*pi*2*doy/365), sin2 = sin(2*pi*2*doy/365),
    cos3 = cos(2*pi*3*doy/365), sin3 = sin(2*pi*3*doy/365),
    cos4 = cos(2*pi*4*doy/365), sin4 = sin(2*pi*4*doy/365),
    cos5 = cos(2*pi*5*doy/365), sin5 = sin(2*pi*5*doy/365),
    cos6 = cos(2*pi*6*doy/365), sin6 = sin(2*pi*6*doy/365),
    cos7 = cos(2*pi*7*doy/365), sin7 = sin(2*pi*7*doy/365),
    cos8 = cos(2*pi*8*doy/365), sin8 = sin(2*pi*8*doy/365),
    cos9 = cos(2*pi*9*doy/365), sin9 = sin(2*pi*9*doy/365),
    cos10 = cos(2*pi*10*doy/365), sin10 = sin(2*pi*10*doy/365)
    
    
  )
temp_ns <- ns(dat_city$rolling_avg_temperature_scaled, df = 5)
colnames(temp_ns) <- paste0("temp_ns", 1:ncol(temp_ns))
dat_city <- bind_cols(dat_city, as.data.frame(temp_ns))

dat_city <- dummy_cols(dat_city, select_columns = "dow")

print(sapply(dat_city, class))   # inspect - fix any factor/character columns below if needed

dat_city <- dat_city %>%
  mutate(holiday = as.numeric(as.character(holiday)))

dow_cols  <- grep("^dow_", names(dat_city), value = TRUE)[-1]  # drop reference level
temp_cols <- grep("^temp_ns", names(dat_city), value = TRUE)

# Shared covariates - no harmonic terms here, those are added per-outcome
# below since respiratory uses 8 harmonic pairs and asthma uses 10,
# matching each outcome's NB model spec.
xreg_cols_base <- c(
  "time_elapsed_scaled",
  "covid_acute",
  "covid_post",
  dow_cols,
  "holiday",
  "rolling_avg_precipitation_scaled",
  temp_cols
)

harmonics_2  <- c("cos1","sin1","cos2","sin2")
#harmonics_10 <- c(harmonics_8, "cos9","sin9","cos10","sin10")

xreg_cols_respiratory <- c(xreg_cols_base, harmonics_2)
xreg_cols_asthma      <- c(xreg_cols_base, harmonics_2)

# Split into training (pre-event) and event window
train_dat <- dat_city %>% filter(date < event_start)
event_dat <- dat_city %>% filter(date >= event_start & date <= event_end)

start_year <- as.numeric(format(min(train_dat$date), "%Y"))
start_doy  <- as.numeric(format(min(train_dat$date), "%j"))

build_xreg <- function(data, cols) {
  data[, cols] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
}

xreg_train_respiratory <- build_xreg(train_dat, xreg_cols_respiratory)
xreg_event_respiratory <- build_xreg(event_dat, xreg_cols_respiratory)

xreg_train_asthma <- build_xreg(train_dat, xreg_cols_asthma)
xreg_event_asthma <- build_xreg(event_dat, xreg_cols_asthma)

stopifnot(
  nrow(xreg_train_respiratory) == nrow(train_dat),
  nrow(xreg_train_asthma) == nrow(train_dat),
  !anyNA(xreg_train_respiratory), !anyNA(xreg_event_respiratory),
  !anyNA(xreg_train_asthma), !anyNA(xreg_event_asthma)
)

fit_arimax_its <- function(outcome_col, xreg_train, xreg_event) {
  
  y_train <- ts(train_dat[[outcome_col]], start = c(start_year, start_doy), frequency = 365)
  
  model <- auto.arima(
    y_train,
    xreg = xreg_train,
    seasonal = TRUE,
    stepwise = FALSE,
    trace = TRUE
  )
  
  pred <- forecast(model, xreg = xreg_event, h = nrow(xreg_event))
  # bootstrap alternative, if residuals are non-normal:
  # pred <- forecast(model, xreg = xreg_event, h = nrow(xreg_event), bootstrap = TRUE, npaths = 1000)
  
  list(model = model, pred = pred)
}

if (file.exists(arimax_respiratory_path)) {
  message("Cached ARIMAX (respiratory) found - loading instead of refitting.")
  arimax_respiratory <- readRDS(arimax_respiratory_path)
} else {
  arimax_respiratory <- fit_arimax_its("respiratory", xreg_train_respiratory, xreg_event_respiratory)
  summary(arimax_respiratory$model)
  checkresiduals(arimax_respiratory$model)
  saveRDS(arimax_respiratory, arimax_respiratory_path)
}

if (file.exists(arimax_asthma_path)) {
  message("Cached ARIMAX (asthma) found - loading instead of refitting.")
  arimax_asthma <- readRDS(arimax_asthma_path)
} else {
  arimax_asthma <- fit_arimax_its("asthma", xreg_train_asthma, xreg_event_asthma)
  summary(arimax_asthma$model)
  checkresiduals(arimax_asthma$model)
  saveRDS(arimax_asthma, arimax_asthma_path)
}


# Safety net only - event_dat is now always built above regardless of
# ARIMAX caching, kept in case the section
# above is ever restructured again.
if (!exists("event_dat")) {
  event_dat <- dat %>% filter(date >= event_start & date <= event_end)
}


# ============================================================================
# 10. ARIMAX excess-count extraction - day-specific AND overall
# ============================================================================
# Day-specific estimates/CIs still come from the forecast object's own
# marginal prediction intervals (pred$lower/pred$upper) - those are fine as
# single-day, marginal intervals.
#
# The Overall CI is NOT the sum of the daily bounds. Summing conf.low across
# days (and conf.high across days) implicitly assumes the low (or high) end
# of the interval occurs on every day simultaneously, which is almost never
# true and produces an artificially wide interval that ignores the actual
# day-to-day correlation in multi-step ARIMA forecast errors. Instead,
# Overall's CI comes from simulating many bootstrapped future paths from the
# fitted model, summing each simulated path's excess across the event
# window, and taking the 2.5th/97.5th percentiles of that distribution of
# sums - this automatically incorporates whatever correlation actually
# exists between days' forecast errors, without assuming normality or
# assuming the days are independent.

simulate_arimax_overall_ci <- function(model, xreg_event, observed_vec,
                                       npaths = ARIMAX_SIM_NPATHS,
                                       seed   = ARIMAX_SIM_SEED) {
  
  n_days <- nrow(xreg_event)
  
  set.seed(seed)
  
  # n_days x npaths matrix - each column is one bootstrap-simulated future
  # path drawn from the fitted model's residual distribution, conditioned
  # on xreg_event (so seasonal/weather/day-of-week structure is preserved
  # in every simulated path, only the residual shocks vary).
  sim_paths <- replicate(
    npaths,
    as.numeric(
      simulate(
        model,
        nsim      = n_days,
        future    = TRUE,
        xreg      = xreg_event,
        bootstrap = TRUE
      )
    )
  )
  
  # Total excess (observed - simulated) for each simulated path, summed
  # across the event window - this is the actual joint distribution of the
  # 3-day total, not an assumption about how the days' errors relate.
  total_excess_sim <- sum(observed_vec) - colSums(sim_paths)
  
  ci <- quantile(total_excess_sim, probs = c(0.025, 0.975))
  
  data.frame(
    conf.low  = unname(ci[1]),
    conf.high = unname(ci[2])
  )
  
}


get_excess_arimax <- function(fit_obj, data_event, outcome_col, xreg_event) {
  
  pred <- fit_obj$pred
  
  day_level <- data.frame(
    day       = event_dates$day,
    obs       = data_event[[outcome_col]],
    pred      = as.numeric(pred$mean),
    pred_lower = as.numeric(pred$lower[, 2]),  # 95% CI
    pred_upper = as.numeric(pred$upper[, 2])
  ) %>%
    mutate(
      estimate  = obs - pred,
      conf.low  = obs - pred_upper,
      conf.high = obs - pred_lower
    ) %>%
    select(day, estimate, conf.low, conf.high)
  
  overall_ci <- simulate_arimax_overall_ci(
    model        = fit_obj$model,
    xreg_event   = xreg_event,
    observed_vec = data_event[[outcome_col]]
  )
  
  print(overall_ci)
  
  overall <- data.frame(
    day       = "Overall",
    estimate  = sum(day_level$estimate),
    conf.low  = overall_ci$conf.low,
    conf.high = overall_ci$conf.high
  )
  
  bind_rows(overall, day_level)
}

excess_respiratory_arima <- get_excess_arimax(
  arimax_respiratory, event_dat, "respiratory", xreg_event_respiratory
) %>%
  mutate(outcome = "Respiratory", model = "Day-specific forecast", source = "ARIMA errors")

excess_asthma_arima <- get_excess_arimax(
  arimax_asthma, event_dat, "asthma", xreg_event_asthma
) %>%
  mutate(outcome = "Asthma", model = "Day-specific forecast", source = "ARIMA errors")


print(excess_respiratory_arima)
print(excess_asthma_arima)




# ============================================================================
# 11. Combine NB and ARIMAX excess estimates for downstream plotting
#     (plotting itself lives in a separate script). The published-design
#     result stays separate as its own RR-only output (rr_published_design).
# ============================================================================

excess_combined <- bind_rows(
  
  excess_results_nb %>%
    select(outcome, model, day, estimate, conf.low, conf.high, source),
  
  excess_respiratory_arima %>%
    select(outcome, model, day, estimate, conf.low, conf.high, source),
  
  excess_asthma_arima %>%
    select(outcome, model, day, estimate, conf.low, conf.high, source)
  
)

print(excess_combined)

saveRDS(excess_combined, file.path(out_dir, "excess_combined.rds"))
write.csv(excess_combined, file.path(out_dir, "excess_combined.csv"), row.names = FALSE)