# Load required packages
library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(tidygeocoder)
library(viridis)
library(purrr)
library(patchwork)

my_wd <- "/Users/ninaflores/Library/CloudStorage/OneDrive-SharedLibraries-UW/casey_cohort - Documents/"
wf_loc <- "data/environmental/wildfires/raw_data/us_2006-2023_zcta-wf-pm/smokePM2pt5_predictions_daily_zcta_20060101-20231231.rds"
fig <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/figures"


# --------------------------------------------------------------------------
# Read in data
# --------------------------------------------------------------------------

# read in wf pm data
wf_dta <- readRDS(paste0(my_wd, wf_loc)) %>%
  filter(date >= as.Date("2023-01-01"))

# read in nyc 2019 zctas for nyc

# Set tigris options
options(tigris_use_cache = TRUE)  # Cache the shapefile for future use

# Pull in the ZCTA shapefiles for 2019
zctas_2019 <- zctas(year = 2019, cb = TRUE, class = "sf")

# Filter for NYC ZCTAs (FIPS codes for NYC counties: Bronx (36005), Kings (36047), New York (36061), Queens (36081), Richmond (36085))
nyc_county_fips <- c("36005", "36047", "36061", "36081", "36085")

# Load county boundaries to intersect with ZCTAs
nyc_counties <- counties(state = "NY", year = 2019, class = "sf", cb = TRUE) %>%
  filter(GEOID %in% nyc_county_fips)

# --------------------------------------------------------------------------
# Intersect ZCTAs with NYC boundaries
# --------------------------------------------------------------------------

nyc_zctas_2019 <- st_intersection(zctas_2019, nyc_counties) %>%
  select(GEOID10, COUNTYFP) %>%
  rename("GEOID" = "GEOID10")

nyc_wf <- nyc_zctas_2019 %>%
  left_join(wf_dta) %>%
  st_transform(2263)

catchment_counties <- nyc_zctas_2019 %>%
  filter(COUNTYFP == "005" | COUNTYFP == "061") %>%
  group_by(COUNTYFP) %>%
  st_union()

# --------------------------------------------------------------------------
# Geocode the address of CHONY
# --------------------------------------------------------------------------

address_df <- data.frame(address = c("3959 Broadway, New York, NY 10032"))

location <- address_df %>%
  geocode(address = address, method = "osm") %>% # Geocode using OpenStreetMap
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% # WGS84 (EPSG:4326)
  st_transform(2263)

# --------------------------------------------------------------------------
# Daily wildfire smoke PM2.5 (event window)
# --------------------------------------------------------------------------

fig1 <- ggplot() +
  geom_sf(data = nyc_wf %>%
            filter(date %in% as.Date(c("2023-06-06",
                                       "2023-06-07",
                                       "2023-06-08"))),
          aes(fill = smokePM_pred, geometry = geometry),
          color = NA) +
  scale_fill_viridis(name = expression("PM"[2.5] ~ "(µg/m³)"), option = "inferno", direction = -1) +
  facet_wrap(~ date, labeller = labeller(date = function(x) gsub(" 0", " ", format(as.Date(x), "%B %d, %Y")))) +  geom_sf(data = location, aes(color = "MS-CHONY"), size = 1.75, shape = 8) +
  scale_color_manual(
    values = c("MS-CHONY" = "#009EFF"),
    name = ""
  ) +
  labs(title = expression("b. Daily wildfire smoke PM"[2.5])) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm"),
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2)
  )

# add in background pm average

nyccas_dir <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/nyccas"

# --------------------------------------------------------------------------
#  Clean NYC ZCTA
# --------------------------------------------------------------------------

nyc_zctas <- nyc_zctas_2019 %>%
  group_by(GEOID) %>%
  summarise(
    geometry = st_union(geometry),
    .groups = "drop"
  )

nyc_zctas_vect <- terra::vect(nyc_zctas)

# --------------------------------------------------------------------------
# Read each annual NYCCAS raster and calculate ZCTA mean
#
# aa8  = 2016
# aa9  = 2017
# ...
# aa15 = 2023
# --------------------------------------------------------------------------

nyccas_years <- tibble(
  year = 2016:2023,
  folder = paste0("aa", 8:15, "_pm300m")
)

nyccas_annual <- map_dfr(seq_len(nrow(nyccas_years)), function(i) {
  
  yr <- nyccas_years$year[i]
  folder <- nyccas_years$folder[i]
  
  message("Processing ", yr)
  
  # Read ArcInfo Grid
  r <- terra::rast(
    file.path(nyccas_dir, folder)
  )
  
  # Put ZCTAs into the raster's CRS
  zctas_r <- terra::project(
    nyc_zctas_vect,
    terra::crs(r)
  )
  
  # Mean PM2.5 within each ZCTA
  ext <- terra::extract(
    r,
    zctas_r,
    fun = mean,
    na.rm = TRUE
  )
  
  tibble(
    GEOID = nyc_zctas$GEOID,
    year = yr,
    nyccas_pm25 = ext[[2]]
  )
})

# --------------------------------------------------------------------------
# Average the annual ZCTA estimates across 2016-2023
# --------------------------------------------------------------------------

nyccas_2016_2023 <- nyccas_annual %>%
  group_by(GEOID) %>%
  summarise(
    nyccas_pm25_2016_2023 = mean(
      nyccas_pm25,
      na.rm = TRUE
    ),
    n_years = sum(!is.na(nyccas_pm25)),
    .groups = "drop"
  )

# Check that all ZCTAs have data for all 8 years
table(nyccas_2016_2023$n_years)

summary(nyccas_2016_2023$nyccas_pm25_2016_2023)

# --------------------------------------------------------------------------
# Join average PM2.5 back to ZCTA geometry
# --------------------------------------------------------------------------

nyccas_avg_sf <- nyc_zctas %>%
  left_join(
    nyccas_2016_2023,
    by = "GEOID"
  )

# --------------------------------------------------------------------------
# Annual average daily PM2.5 (2016-2023)
# --------------------------------------------------------------------------

fig_nyccas <- ggplot() +
  geom_sf(
    data = nyccas_avg_sf,
    aes(fill = nyccas_pm25_2016_2023, geometry = geometry),
    color = NA
  ) +
  scale_fill_viridis(
    name = expression("PM"[2.5] ~ "(µg/m³)"),
    option = "inferno",
    direction = -1
  ) +
  geom_sf(data = location, aes(color = "MS-CHONY"), size = 1.75, shape = 8) +
  scale_color_manual(
    values = c("MS-CHONY" = "#009EFF"),
    name = ""
  ) +
  labs(title = expression("a. Average daily PM"[2.5]*" (2016-2023)")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(1, "cm")
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2)
  )

fig_nyccas

# --------------------------------------------------------------------------
# Combine
# --------------------------------------------------------------------------

fig_combined <- fig_nyccas + fig1 +
  plot_layout(
    ncol = 2,
    widths = c(1.1, 3)
  )

fig_combined


ggsave(
  filename = file.path(fig, "pm_maps.pdf"),
  plot = fig_combined,
  width = 12.5, height = 5, dpi = 300
)
