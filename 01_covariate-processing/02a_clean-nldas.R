require(data.table)
require(fst)
require(dplyr)
require(lubridate)

setwd("~/Desktop/projects/casey cohort/nyc wildfires/data/cov")
a <- fread("nldas-zcta-2016.csv")
b <- fread("nldas-zcta-2017.csv")
c <- fread("nldas-zcta-2018.csv")
d <- fread("nldas-zcta-2019.csv")
e <- fread("nldas-zcta-2020.csv")
f <- fread("nldas-zcta-2021.csv")
g <- fread("nldas-zcta-2022.csv")
h <- fread("nldas-zcta-2023.csv")


nldas <- rbind(a,b,c,d,e, f, g, h)

rm(a,b,c,d,e,f,g,h)


setDT(nldas)
# Split "system:index" column and drop the leading 'A' from date
nldas <- nldas[, c("date", "hour") := tstrsplit(`system:index`, "_", fixed = TRUE)[1:2]]
# Drop the first letter of the date string
nldas <- nldas[, date := substr(date, 2, nchar(date))]
# Convert date to Date class
nldas <-nldas[, date := as.Date(as.character(date), format = "%Y%m%d")]

##' Convert specific humidity to relative humidity
##'
##' converting specific humidity into relative humidity
##' NCEP surface flux data does not have RH
##' from Bolton 1980 The computation of Equivalent Potential Temperature
##' \url{http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html}
##' @title qair2rh
##' @param qair specific humidity, dimensionless (e.g. kg/kg) ratio of water mass / total air mass
##' @param temp degrees C
##' @param press pressure in mb
##' @return rh relative humidity, ratio of actual water mixing ratio to saturation mixing ratio
##' @export
##' @author David LeBauer
qair2rh <- function(qair, temp, press){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  e <- qair * press / (0.378 * qair + 0.622)
  rh <- e / es
  rh[rh > 1] <- 1
  rh[rh < 0] <- 0
  return(rh)
}



absolute_humidity <- function(air_temp, rh) {
  # https://carnotcycle.wordpress.com/2012/08/04/how-to-convert-relative-humidity-to-absolute-humidity/
  (6.112 * exp((17.67 * air_temp) / (air_temp + 243.5)) * rh * 2.1674) / 
    (273.15 + air_temp)
}


nldas <-nldas %>%
  mutate(pres_mb = .01 * pressure,
         rh = qair2rh(specific_humidity, temperature, pres_mb),
         rh_percentage = rh *100,
         ah = absolute_humidity(temperature, rh_percentage)) %>%
  group_by(date, hour, ZCTA5CE10) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(ZCTA5CE10, date, hour, temperature, specific_humidity,rh, ah, pressure, total_precipitation)

write.fst(nldas, "nldas_cov_2017_2023.fst")

head(nldas)


