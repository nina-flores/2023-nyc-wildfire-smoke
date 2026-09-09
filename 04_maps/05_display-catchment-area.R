###########***********************
#### Code Description ####
# Author: Nina
# Date: 7/13/26
# Goal: Show the catchment area of the hospital
####**********************

# Load packages 
rm(list=ls())
library(here)
library(tigris)
library(datetimeutils)
library(tis)
library(lubridate)
library(patchwork)
source(paste0("~/Documents/repos/nyc-wildfire-smoke/02_outcome-processing/0_read_libraries.R"))

# read in data -------------------------------------------------------------------------

out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
fig <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/figures"


dat <- read.fst(paste0(out,"/","analytical_data.fst")) %>%
  group_by(zcta) %>%
  summarize(
    sum_visits = sum(respiratory),
    sum_asthm = sum(asthma)) %>%
  ungroup() %>%
  mutate(sum_visits_all = sum(sum_visits)) %>%
  group_by(zcta) %>%
  mutate(percent_resp_visits = 100*(sum_visits/sum_visits_all))


# Get all ZCTAs as an sf object
zct <- zctas(state = "NY", year = 2010, class = "sf") %>%  # Filter for New York state
  mutate(zcta = ZCTA5CE10) %>%
  select(zcta)

visits_zct <- dat %>%
  left_join(zct) %>%
  st_as_sf() %>%
  mutate(
    percent_resp_visits_bin = cut(
      percent_resp_visits,
      breaks = c(-Inf, .5, 4, 8, Inf),
      labels = c("< 0.5%", "0.5-4%", "4-8%", "8%+"),
      right = FALSE
    )
  )

# Create a point for the CHONY
hospital_point <- st_sf(
  name = "CUIMC",
  geometry = st_sfc(st_point(c(-73.9435, 40.8404)), crs = 4326)
)



p1 <- ggplot(visits_zct) +
  geom_sf(aes(fill = percent_resp_visits_bin), linewidth = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "inferno", direction = -1,begin = 0.35, end = 1, na.value = "grey90") +
  geom_sf(data = hospital_point, aes(color = name), fill = "#009EFF", shape = 8, size = 1.75) +
  scale_color_manual(values = c("CUIMC" = "#009EFF"), name = NULL, labels = "MS-CHONY") +
  labs(fill = "% visits") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  ggtitle("a. % respiratory visits, by ZCTA")


dat_nevi <- read.fst(paste0(out, "/", "analytical_data.fst")) %>%
  group_by(zcta) %>%
  select(nevi_quartile) %>%
  slice(1)

nevi_zct <- dat_nevi %>%
  left_join(zct) %>%
  st_as_sf() %>%
  mutate(nevi_quartile = factor(nevi_quartile, levels = c("1", "2", "3", "4")))

p2 <- ggplot(nevi_zct) +
  geom_sf(aes(fill = nevi_quartile), linewidth = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "inferno", direction = -1,begin = 0.35, end = 1, na.value = "grey90") +
  geom_sf(data = hospital_point, aes(color = name), fill = "#009EFF", shape = 8, size = 1.75) +
  scale_color_manual(values = c("CUIMC" = "#009EFF"), name = NULL, labels = "MS-CHONY") +
  labs(fill = "NEVI quartile") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  ggtitle("c. NEVI quartile of ZCTA")

# Combine side by side
p1 + p2


visits_zct_filtered <- visits_zct %>%
  mutate(catchment = if_else(percent_resp_visits >= .5, "Yes", "No"))

inferno_2 <- viridis(2, option = "inferno", direction = -1, begin =1, end = 0.45)

p3 <- ggplot(visits_zct_filtered) +
  geom_sf(aes(fill = catchment), linewidth = 0.5) +
  scale_fill_manual(values = c("Yes" = inferno_2[1], "No" = inferno_2[2]))+
  geom_sf(data = hospital_point, aes(color = name), fill = "#009EFF", shape = 8, size = 1.75) +
  scale_color_manual(values = c("CUIMC" = "#009EFF"), name = NULL, labels = "MS-CHONY") +
  theme_minimal(base_size = 12) +
  labs(fill = "") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  ggtitle("b. Included ZCTAs")

# Save combined figure -------------------------------------------------------------------------
ggsave(
  filename = file.path(fig, "catchment_area_maps.pdf"),
  plot = p1 / p3 / p2,
  width = 5, height = 8.5, dpi = 300
)

