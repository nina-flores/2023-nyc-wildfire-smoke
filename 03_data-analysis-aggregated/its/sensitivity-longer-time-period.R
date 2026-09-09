#### ARIMAX 2-week projection - excess visits overall + daily
#### Reads in the FINAL fitted ARIMAX models (cached .rds, fit on pre-event
#### data only) and forecasts a 14-day window instead of the original 3-day
#### event window, to get daily excess counts + an overall 2-week excess
#### total across that longer horizon.
####***********************

rm(list = ls())

# Packages ----------------------------------------------------------------

library(here)
library(dplyr)
library(purrr)
library(splines)
library(fastDummies)
library(forecast)
library(fst)
library(forcats)
library(tsModel)

source("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R")

# Paths ---------------------------------------------------------------------

an      <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
out_dir <- file.path(an, "models")

arimax_respiratory_path <- file.path(out_dir, "arimax_respiratory.rds")
arimax_asthma_path      <- file.path(out_dir, "arimax_asthma.rds")

stopifnot(file.exists(arimax_respiratory_path), file.exists(arimax_asthma_path))

arimax_respiratory <- readRDS(arimax_respiratory_path)
arimax_asthma      <- readRDS(arimax_asthma_path)

# 2-week projection window ---------------------------------------------------
# ASSUMPTION: starts at the original event start and runs 14 days. Change
# proj_start/proj_end here if you meant a different window.

proj_start <- as.Date("2023-06-06")
proj_end   <- proj_start + 9   # 10 days total

proj_dates <- data.frame(
  day  = format(seq(proj_start, proj_end, by = "day"), "%m/%d"),
  date = seq(proj_start, proj_end, by = "day")
)

ARIMAX_SIM_NPATHS <- 5000
ARIMAX_SIM_SEED   <- 1

# ============================================================================
# 1. Read and rebuild the same xreg design used at training time
# ============================================================================
# Must exactly match the harmonic count, dow dummies, temp spline, and
# column set used when arimax_respiratory / arimax_asthma were originally
# fit, or forecast() will silently misalign coefficients.

dat <- read.fst(
  file.path(an, "analytical_data_prepped-updated-agg.fst")
) %>%
  arrange(date) %>%
  mutate(date_num = as.integer(date - min(date)) + 1)

dat_city <- dat %>%
  filter(!(format(date, "%m") == "02" & format(date, "%d") == "29")) %>%
  mutate(
    cos1 = cos(2*pi*1*doy/365), sin1 = sin(2*pi*1*doy/365),
    cos2 = cos(2*pi*2*doy/365), sin2 = sin(2*pi*2*doy/365)
  )

temp_ns <- ns(dat_city$rolling_avg_temperature_scaled, df = 5)
colnames(temp_ns) <- paste0("temp_ns", 1:ncol(temp_ns))
dat_city <- bind_cols(dat_city, as.data.frame(temp_ns))

dat_city <- dummy_cols(dat_city, select_columns = "dow") %>%
  mutate(holiday = as.numeric(as.character(holiday)))

dow_cols  <- grep("^dow_", names(dat_city), value = TRUE)[-1]
temp_cols <- grep("^temp_ns", names(dat_city), value = TRUE)

xreg_cols_base <- c(
  "time_elapsed_scaled",
  "covid_acute",
  "covid_post",
  dow_cols,
  "holiday",
  "rolling_avg_precipitation_scaled",
  temp_cols
)

harmonics_2 <- c("cos1", "sin1", "cos2", "sin2")

xreg_cols_respiratory <- c(xreg_cols_base, harmonics_2)
xreg_cols_asthma      <- c(xreg_cols_base, harmonics_2)

build_xreg <- function(data, cols) {
  data[, cols] %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
}

proj_dat <- dat_city %>% filter(date >= proj_start & date <= proj_end)

stopifnot(nrow(proj_dat) == nrow(proj_dates))

xreg_proj_respiratory <- build_xreg(proj_dat, xreg_cols_respiratory)
xreg_proj_asthma      <- build_xreg(proj_dat, xreg_cols_asthma)

stopifnot(
  !anyNA(xreg_proj_respiratory),
  !anyNA(xreg_proj_asthma)
)

# ============================================================================
# 2. Forecast the 14-day window from each fitted model
# ============================================================================

forecast_2wk <- function(fit_obj, xreg_proj) {
  forecast(fit_obj$model, xreg = xreg_proj, h = nrow(xreg_proj))
}

pred_respiratory <- forecast_2wk(arimax_respiratory, xreg_proj_respiratory)
pred_asthma      <- forecast_2wk(arimax_asthma, xreg_proj_asthma)

# ============================================================================
# 3. Overall (14-day total) excess CI via bootstrap simulation
# ============================================================================
# Same logic as the original 3-day Overall CI: simulate many bootstrapped
# future paths from the fitted model's residual distribution, sum each
# path's excess across all 14 days, and take the 2.5th/97.5th percentiles
# of that distribution - properly propagates day-to-day forecast-error
# correlation over the longer horizon instead of summing daily bounds.

simulate_arimax_overall_ci <- function(model, xreg_proj, observed_vec,
                                       npaths = ARIMAX_SIM_NPATHS,
                                       seed   = ARIMAX_SIM_SEED) {
  
  n_days <- nrow(xreg_proj)
  
  set.seed(seed)
  
  sim_paths <- replicate(
    npaths,
    as.numeric(
      simulate(
        model,
        nsim      = n_days,
        future    = TRUE,
        xreg      = xreg_proj,
        bootstrap = TRUE
      )
    )
  )
  
  total_excess_sim <- sum(observed_vec) - colSums(sim_paths)
  
  ci <- quantile(total_excess_sim, probs = c(0.025, 0.975))
  
  data.frame(
    conf.low  = unname(ci[1]),
    conf.high = unname(ci[2])
  )
  
}

# ============================================================================
# 4. Daily + overall excess-visit table
# ============================================================================

get_excess_2wk <- function(fit_obj, data_proj, outcome_col, xreg_proj, dates_df) {
  
  pred <- forecast_2wk(fit_obj, xreg_proj)
  
  day_level <- data.frame(
    day        = dates_df$day,
    date       = dates_df$date,
    obs        = data_proj[[outcome_col]],
    pred       = as.numeric(pred$mean),
    pred_lower = as.numeric(pred$lower[, 2]),  # 95% CI
    pred_upper = as.numeric(pred$upper[, 2])
  ) %>%
    mutate(
      estimate  = obs - pred,
      conf.low  = obs - pred_upper,
      conf.high = obs - pred_lower
    ) %>%
    select(day, date, estimate, conf.low, conf.high)
  
  overall_ci <- simulate_arimax_overall_ci(
    model        = fit_obj$model,
    xreg_proj    = xreg_proj,
    observed_vec = data_proj[[outcome_col]]
  )
  
  overall <- data.frame(
    day       = "Overall",
    date      = NA,
    estimate  = sum(day_level$estimate),
    conf.low  = overall_ci$conf.low,
    conf.high = overall_ci$conf.high
  )
  
  bind_rows(overall, day_level)
}

excess_respiratory_2wk <- get_excess_2wk(
  arimax_respiratory, proj_dat, "respiratory", xreg_proj_respiratory, proj_dates
) %>%
  mutate(outcome = "Respiratory", window = "2-week projection", source = "ARIMA errors")

excess_asthma_2wk <- get_excess_2wk(
  arimax_asthma, proj_dat, "asthma", xreg_proj_asthma, proj_dates
) %>%
  mutate(outcome = "Asthma", window = "2-week projection", source = "ARIMA errors")

print(excess_respiratory_2wk)
print(excess_asthma_2wk)

excess_combined_2wk <- bind_rows(excess_respiratory_2wk, excess_asthma_2wk)

print(excess_combined_2wk)

saveRDS(excess_combined_2wk, file.path(out_dir, "excess_combined_2wk.rds"))
write.csv(excess_combined_2wk, file.path(out_dir, "excess_combined_2wk.csv"), row.names = FALSE)


excess_plot_2wk <- excess_combined_2wk %>%
  mutate(
    outcome = factor(outcome, levels = c("Respiratory", "Asthma")),
    day     = factor(day, levels = c("Overall", proj_dates$day))
  ) %>%
  ggplot(
    aes(
      x = day,
      y = estimate,
      ymin = conf.low,
      ymax = conf.high
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_pointrange(
    size = 0.5,
    color = "black"
  ) +
  
  facet_wrap(
    ~ outcome,
    scales = "free_y"
  ) +
  
  labs(
    x = NULL,
    y = "Excess visits (95% CI)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

excess_plot_2wk


#### 

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an ,"/","analytical_data_prepped-updated-agg.fst")) %>%
  mutate(smoke_day = as.factor(if_else(date %in% c("2023-06-06",
                                                   "2023-06-07",
                                                   "2023-06-08",
                                                   "2023-06-09",
                                                   "2023-06-10",
                                                   "2023-06-11",
                                                   "2023-06-12",
                                                   "2023-06-13",
                                                   "2023-06-14",
                                                   "2023-06-15",
                                                   "2023-06-16",
                                                   "2023-06-17",
                                                   "2023-06-18",
                                                   "2023-06-19") ,1,0))) %>%
  mutate(smoke_day1 = as.factor(if_else(date %in% c("2023-06-06") ,1,0)),
         smoke_day2 = as.factor(if_else(date %in% c("2023-06-07") ,1,0)),
         smoke_day3 = as.factor(if_else(date %in% c("2023-06-09") ,1,0)),
         smoke_day4 = as.factor(if_else(date %in% c("2023-06-10") ,1,0)),
         smoke_day5 = as.factor(if_else(date %in% c("2023-06-11") ,1,0)),
         smoke_day6 = as.factor(if_else(date %in% c("2023-06-12") ,1,0)),
         smoke_day7 = as.factor(if_else(date %in% c("2023-06-13") ,1,0)),
         smoke_day8 = as.factor(if_else(date %in% c("2023-06-14") ,1,0)),
         smoke_day9 = as.factor(if_else(date %in% c("2023-06-15") ,1,0)),
         smoke_day10 = as.factor(if_else(date %in% c("2023-06-16") ,1,0)),
         smoke_day11 = as.factor(if_else(date %in% c("2023-06-17") ,1,0)),
         smoke_day12 = as.factor(if_else(date %in% c("2023-06-18") ,1,0)),
         smoke_day13 = as.factor(if_else(date %in% c("2023-06-19") ,1,0)))
  
  


mod<- glmmTMB(respiratory ~ smoke_day + time_elapsed_scaled +
                as.factor(holiday) + harmonic(doy,2,365)  +
                ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                covid_post + covid_acute ,
              family = nbinom1(),
              data = dat)
summary(mod)


mod2 <- glmmTMB(respiratory ~ smoke_day1 +
                  smoke_day2 + 
                  smoke_day3 +
                  smoke_day4 +
                  smoke_day5 +
                  smoke_day6 +
                  smoke_day7 +
                  smoke_day8 +
                  smoke_day9 +
                  smoke_day10 +
                  smoke_day11 +
                  smoke_day12 +
                  smoke_day13 +
                time_elapsed_scaled +
                as.factor(holiday) + harmonic(doy,2,365)  +
                ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                covid_post + covid_acute ,
              family = nbinom1(),
              data = dat)
summary(mod2)


mod_asthma <- glmmTMB(asthma ~ smoke_day + time_elapsed_scaled +
                as.factor(holiday) + harmonic(doy,2,365)  +
                ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                covid_post + covid_acute ,
              family = nbinom1(),
              data = dat)
summary(mod_asthma)


mod_asthma2 <- glmmTMB(asthma ~ smoke_day1 +
                  smoke_day2 + 
                  smoke_day3 +
                  smoke_day4 +
                  smoke_day5 +
                  smoke_day6 +
                  smoke_day7 +
                  smoke_day8 +
                  smoke_day9 +
                  smoke_day10 +
                  smoke_day11 +
                  smoke_day12 +
                  smoke_day13 +
                  time_elapsed_scaled +
                  as.factor(holiday) + harmonic(doy,2,365)  +
                  ns(rolling_avg_temperature_scaled,5) + rolling_avg_precipitation_scaled +
                  covid_post + covid_acute ,
                family = nbinom1(),
                data = dat)
summary(mod_asthma2)

