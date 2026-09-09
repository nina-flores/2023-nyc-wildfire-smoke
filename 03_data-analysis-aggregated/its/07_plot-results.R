###########***********************
#### Replot wildfire smoke ITS results
#### Rename model labels only
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
# The analysis script saves "excess_combined.rds", not
# "plot_dat_combined.rds" 
# Columns: outcome, model, day, estimate, conf.low, conf.high, source.

plot_dat_combined <- readRDS(
  file.path(out_dir, "excess_combined.rds")
) %>%
  mutate(day = factor(day, levels = c("Overall","6/6", "6/7", "6/8")),
         outcome = factor(outcome, levels = c("Respiratory", "Asthma")))




# ============================================================================
# Update labels
# ============================================================================


plot_dat_combined <- plot_dat_combined %>%
  mutate(
    
    source = recode(
      source,
      "NB (nbinom1)" = "GLM",
      "ARIMA errors" = "ARIMAX"
    ),
    
    
    source = factor(
      source,
      levels = c(
        "GLM",
        "ARIMAX"
      )
    )
    
  )


# Check updated data
print(plot_dat_combined)


# Save updated version
saveRDS(
  plot_dat_combined,
  file.path(out_dir, "plot_dat_combined_relabelled.rds")
) 



# ============================================================================
# Plot
# ============================================================================

figure <- ggplot(
  plot_dat_combined,
  aes(
    x = day,
    y = estimate,
    ymin = conf.low,
    ymax = conf.high,
    color = source
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_pointrange(
    size = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(
    ~ outcome,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "GLM" = "#88419d",
      "ARIMAX" = "black"
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 5, only.loose = FALSE)
  )+
  labs(
    x = NULL,
    y = "Excess visits (95% CI)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

figure

fig <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/figures"

ggsave(
  filename = file.path(fig, "reg_results.pdf"),
  plot = figure,
  width = 12, height = 5, dpi = 300
)