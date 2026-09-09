###########***********************
#### Code Description ####
# Goal: Recreate model_4_mult, aggregated to a single citywide daily series,
#       using auto.arima() to identify the temporal correlation structure.
#
#   model_4_mult <- glmer(respiratory ~ smoke_day + time_elapsed_scaled +
#                            harmonic(doy,4,365) + dow + holiday +
#                            rolling_avg_precipitation_scaled + caserate_scaled +
#                            ns(x,3) + ns(y,3) + (1 | zcta),
#                          data = dat, family = "poisson", nAGQ = 0)
#
# WHAT CHANGES WHEN AGGREGATING TO A SINGLE CITYWIDE SERIES:
#   - (1 | zcta)          -> dropped. auto.arima fits one series; there is
#                            no panel/random-effect structure once collapsed.
#   - ns(x, 3), ns(y, 3)  -> dropped. These are zcta-level spatial coordinates
#                            and have no meaningful citywide equivalent once
#                            ummed/aggregated across zcta.
#   - Everything else (smoke_day, time_elapsed_scaled, harmonic(doy,4,365),
#     dow, holiday, rolling_avg_precipitation_scaled, caserate_scaled) is
#     date-level (identical across zcta on a given day) or aggregatable, and
#     is recreated below.
#
# NOTE ON HARMONIC(doy, 4, 365): this specifies 4 harmonic pairs (8 columns:
# cos1/sin1 ... cos4/sin4), matching harmonic() from the tsModel package.
####**********************

# 0. Setup ------------------------------------------------------------------
library(forecast)   # auto.arima, checkresiduals
library(dplyr)
library(fastDummies) # dummy_cols
library(fst)

source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
dat <- read.fst(paste0(an, "/", "analytical_data_prepped-updated.fst"))


# 1. Collapse the panel to a single citywide daily series --------------------
# Aggregate RAW variables only here - build harmonics AFTER aggregating
# (see section 2), not before, so they reflect the citywide series itself.
#
# date-level variables (doy, dow, holiday, time_elapsed_scaled) are identical
# across zcta on a given day, so first() is appropriate and introduces no
# aggregation-order concern.

dat_city <- dat %>%
  group_by(date) %>%
  summarise(
    respiratory = sum(respiratory, na.rm = TRUE),
    smoke_day    = first(smoke_day, na.rm = TRUE),        
    caserate_scaled = mean(caserate_scaled, na.rm = TRUE),
    rolling_avg_precipitation_scaled = mean(rolling_avg_precipitation_scaled, na.rm = TRUE),
    time_elapsed_scaled = first(time_elapsed_scaled),
    holiday = first(holiday),
    doy     = first(doy),
    dow     = first(dow),
    .groups = "drop"
  ) %>%
  arrange(date)

# 2. Build harmonic(doy, 4, 365) terms on the aggregated series -------------
dat_city <- dat_city %>%
  mutate(
    cos1 = cos(2 * pi * 1 * doy / 365), sin1 = sin(2 * pi * 1 * doy / 365),
    cos2 = cos(2 * pi * 2 * doy / 365), sin2 = sin(2 * pi * 2 * doy / 365),
    cos3 = cos(2 * pi * 3 * doy / 365), sin3 = sin(2 * pi * 3 * doy / 365),
    cos4 = cos(2 * pi * 4 * doy / 365), sin4 = sin(2 * pi * 4 * doy / 365)
  )

dat_city <- dummy_cols(dat_city, select_columns = "dow")


# 3. Force everything numeric before building xreg ---------------------------
# (glmer's as.factor(holiday)/dow suggests these may be stored as
# factor/character - as.matrix() on a data frame coerces the WHOLE matrix
# to one type if any single column isn't numeric, which silently breaks
# auto.arima. Check and fix explicitly.)

print(sapply(dat_city, class))   # inspect - fix any factor/character columns below if needed

dat_city <- dat_city %>%
  mutate(
    smoke_day = as.numeric(smoke_day)-1,
    holiday   = as.numeric(as.character(holiday))   # handles factor-coded 0/1 - adjust manually if "Yes"/"No"
  )


# 4. Build the ts object and xreg matrix -------------------------------------
start_year <- as.numeric(format(min(dat_city$date), "%Y"))
y <- ts(dat_city$respiratory, start = c(start_year, 1), frequency = 365)

dow_cols <- grep("^dow_", names(dat_city), value = TRUE)

xreg_cols <- c(
  "smoke_day",
  "time_elapsed_scaled",
  "cos1", "sin1", "cos2", "sin2", "cos3", "sin3", "cos4", "sin4",
  dow_cols[-1],     # drop one dow dummy as reference level
  "holiday",
  "rolling_avg_precipitation_scaled",
  "caserate_scaled"
)

xreg <- dat_city[, xreg_cols] %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

stopifnot(is.numeric(xreg), nrow(xreg) == length(y))
anyNA(xreg)  # confirm no NAs were introduced


# 5. auto.arima ---------------------------------------------------------------
model <- auto.arima(
  y,
  xreg = xreg,
  seasonal = TRUE,
  stepwise = FALSE,
  trace = TRUE
)

summary(model)
checkresiduals(model)
confint(model)

# The selected (p,d,q) order tells you the temporal correlation structure
# present after controlling for smoke_day, trend, 4-harmonic seasonality,
# dow, holiday, precipitation, and case rate. E.g.:
#   - ARMA(1,0,0) / AR1  -> supports using ar1() in glmmTMB for the full
#                           zcta-level panel model
#   - ARMA(2,0,0) or higher -> a single ar1() term likely under-corrects;
#                           consider higher-order structure or ou()
#   - Little/no remaining AR/MA structure -> the harmonic(doy,4,365) +
#     time_elapsed_scaled specification may already adequately control
#     for serial correlation
