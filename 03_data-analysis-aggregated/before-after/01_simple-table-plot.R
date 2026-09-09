library(dplyr)
library(tidyr)
library(ggplot2)
library(fst)

# read in data -------------------------------------------------------------------------
out <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/outcome/processed-outcome"
an  <- "~/Desktop/projects/casey cohort/nyc-wildfires/data/analytical-dat"
fig <- "/Users/ninaflores/Desktop/projects/casey cohort/nyc-wildfires/figures"

strip_leading_zero <- function(x) gsub("(?<=/|^)0", "", x, perl = TRUE)

tab1 <- read.fst(paste0(an, "/", "analytical_data_prepped-updated-agg.fst")) %>%
  mutate(date = as.Date(date)) %>%
  mutate(group = case_when(
    date >= as.Date("2023-05-30") & date <= as.Date("2023-06-01") ~ "Before",
    date >= as.Date("2023-06-06") & date <= as.Date("2023-06-08") ~ "During",
    date >= as.Date("2023-06-13") & date <= as.Date("2023-06-15") ~ "After",
    TRUE ~ "Other"
  )) %>%
  filter(group != "Other") %>%
  group_by(group, date) %>%
  summarize(asthma_total = sum(asthma),
            respiratory_total = sum(respiratory),
            .groups = "drop") %>%
  mutate(group = factor(group, levels = c("Before", "During", "After"))) %>%
  arrange(group, date) %>%
  mutate(day_num = row_number(), .by = group,
         bar_label = strip_leading_zero(format(date, "%m/%d")))

tab1_overall <- tab1 %>%
  group_by(group) %>%
  summarize(asthma_total = sum(asthma_total),
            respiratory_total = sum(respiratory_total),
            min_date = min(date),
            max_date = max(date),
            .groups = "drop") %>%
  mutate(day_num = 0,
         bar_label = paste0(strip_leading_zero(format(min_date, "%m/%d")), "-",
                            strip_leading_zero(format(max_date, "%m/%d")))) %>%
  select(-min_date, -max_date)

tab1_combined <- bind_rows(tab1 %>% select(-date), tab1_overall) %>%
  mutate(group_label = recode(group,
                              "Before" = "Week prior",
                              "During" = "During",
                              "After"  = "Week after"),
         group_label = factor(group_label, levels = c("Week prior", "During", "Week after")))

facet_labels <- c("0" = "Overall", "1" = "Day 1", "2" = "Day 2", "3" = "Day 3")

tab1_long <- tab1_combined %>%
  pivot_longer(cols = c(asthma_total, respiratory_total),
               names_to = "outcome",
               values_to = "count") %>%
  mutate(outcome = recode(outcome,
                          "asthma_total" = "Asthma",
                          "respiratory_total" = "Respiratory"),
         outcome = factor(outcome, levels = c("Respiratory", "Asthma")))

# read in published-design IRR results -----------------------------------------------------
# read in published-design IRR results -----------------------------------------------------
rr <- readRDS(file.path(an, "models", "rr_published_design.rds")) %>%
  mutate(
    day_num = case_when(
      day == "Overall" ~ 0,
      day == "6/6"      ~ 1,
      day == "6/7"      ~ 2,
      day == "6/8"      ~ 3
    ),
    outcome = factor(outcome, levels = c("Respiratory", "Asthma")),
    irr_label = sprintf("IRR = %.2f (%.2f-%.2f)", estimate, conf.low, conf.high)  )




inferno_pal <- viridis(10, option = "inferno", direction = -1, begin = 0.15, end = 1)

# preview to see them
scales::show_col(inferno_pal)

dark   <- inferno_pal[8]   # darker purple, near the "begin" side
light1 <- inferno_pal[1]
light2 <- inferno_pal[4]  # lighter end


p <- ggplot(tab1_long, aes(x = group_label, y = count, fill = group_label)) +
  geom_col() +
  geom_text(aes(label = bar_label), vjust = -0.4, size = 3.2) +
  geom_label(
    data = rr,
    aes(x = "During", y = Inf, label = irr_label),
    inherit.aes = FALSE,
    vjust = 1.1, size = 3.2,
    label.size = 0.3, label.r = unit(0.2, "lines"),
    fill = "white", alpha = 0.85, lineheight = 0.9
  ) +
  facet_grid(outcome ~ day_num, scales = "free_y",
             labeller = labeller(day_num = facet_labels)) +
  scale_fill_manual(values = c("Week prior" = light1,
                               "During"     = dark,
                               "Week after" = light2)) +
  labs(x = NULL, y = "Number of ED visits", fill = "Period") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28)))

p

# Save combined figure -------------------------------------------------------------------------
ggsave(
  filename = file.path(fig, "simple_before_after_plot.pdf"),
  plot = p,
  width = 11, height = 8, dpi = 300
)
