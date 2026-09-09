###########***********************
#### Replot wildfire smoke ITS results
#### ARIMAX only
####***********************
rm(list = ls())

# Packages ----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(forcats)

# Paths -------------------------------------------------------------------
an <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
out_dir <- file.path(an, "models")

# ============================================================================
# Read saved results
# ============================================================================
# Columns: outcome, model, day, estimate, conf.low, conf.high, source.
plot_dat_combined <- readRDS(
  file.path(out_dir, "excess_combined.rds")
) %>%
  mutate(day = factor(day, levels = c("Overall", "6/6", "6/7", "6/8")),
         outcome = factor(outcome, levels = c("Respiratory", "Asthma")))

# ============================================================================
# Update labels + filter to ARIMAX only
# ============================================================================
plot_dat_arimax <- plot_dat_combined %>%
  mutate(
    source = recode(
      source,
      "NB (nbinom1)" = "GLM",
      "ARIMA errors" = "ARIMAX"
    )
  ) %>%
  filter(source == "ARIMAX")

# Check updated data
print(plot_dat_arimax)

# Save filtered version
saveRDS(
  plot_dat_arimax,
  file.path(out_dir, "plot_dat_arimax_only.rds")
)

# ============================================================================
# Plot
# ============================================================================
figure <- ggplot(
  plot_dat_arimax,
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
    size = 0.7,
    color = "black"
  ) +
  facet_wrap(
    ~ outcome
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 5, only.loose = FALSE)
  ) +
  labs(
    x = NULL,
    y = "Excess visits (95% CI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

figure

fig <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/figures"
ggsave(
  filename = file.path(fig, "reg_results_arimax_only.pdf"),
  plot = figure,
  width = 12, height = 5, dpi = 300
)