require(data.table)
require(fst)
require(dplyr)
require(lubridate)

setwd("~/Desktop/projects/casey cohort/nyc-wildfires/data/cov")
nldas <- read.fst("nldas_cov_2017_2023.fst")

# Convert data frame to data.table
dt <- as.data.table(nldas)

# Combine date and hour columns into a single string and convert to datetime
dt[, date_time_utc := ymd_hm(paste(date, hour), tz = "UTC")]

# Convert the datetime from UTC to Eastern Time
dt[, date_time_eastern := with_tz(date_time_utc, tzone = "America/New_York")]

dt[, date_eastern := date(date_time_eastern)]


vars_to_avg <- c("temperature", 
                 "specific_humidity",
                 "rh",
                # "ah", 
                 "pressure",
                 "total_precipitation")

# Group by date and zcta, and calculate the average for each variable
agg_data <- dt[, lapply(.SD, mean), by = .(date_eastern, ZCTA5CE10), .SDcols = vars_to_avg]

write.fst(agg_data, "nldas_daily.fst")


