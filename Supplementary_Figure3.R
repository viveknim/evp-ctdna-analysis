# ============================================================================
# Supplementary Figure 3 — ctDNA Collection Characteristics
# EVP ctDNA Project
#
# Cohort: First-line EVP patients with non-zero baseline ctDNA AND
#         ≥1 post-EVP ctDNA draw within 3 months of EVP initiation
#         (same cohort as Figures 3 and 4)
#
# Panel A: Distribution of baseline ctDNA values (log10 scale)
# Panel B: Number of ctDNA draws per patient (all draws)
# Panel C: Time intervals between consecutive ctDNA draws
# Panel D: Time from EVP initiation to first post-EVP ctDNA draw
#
# Layout: 2×2 grid, 10 × 8 inches
# ============================================================================

library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(lubridate)
library(cowplot)

# ---- USER: Update these paths -----------------------------------------------
data_dir      <- "data"
outcomes_path <- file.path(data_dir, "ev_ctdna_survival_tab.csv")
ctdna_path    <- file.path(data_dir, "ctdna_data.csv")
output_dir    <- "output"
# -----------------------------------------------------------------------------


# =============================================================================
# Shared theme and external title helper
# =============================================================================

theme_pub <- function() {
  theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 10, hjust = 0),
      plot.subtitle    = element_text(size = 8, color = "gray40", hjust = 0),
      axis.title       = element_text(face = "bold", size = 9),
      axis.text        = element_text(size = 8, color = "black"),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt")
    )
}

add_panel_title <- function(panel_grid, title_text, subtitle_text = NULL) {
  t_strip <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 10,
               x = 0.02, hjust = 0, vjust = 0.5)
  if (!is.null(subtitle_text) && nchar(subtitle_text) > 0) {
    s_strip <- ggdraw() +
      draw_label(subtitle_text, color = "gray40", size = 7,
                 x = 0.02, hjust = 0, vjust = 0.5)
    plot_grid(t_strip, s_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 0.75, 9))
  } else {
    plot_grid(t_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 10))
  }
}


# =============================================================================
# 1. Load outcomes and first-line EVP cohort
# =============================================================================

outcomes <- read_csv(outcomes_path) %>%
  mutate(evp_c1d1 = as.Date(evp_c1d1, format = "%m/%d/%y"))

fl_evp <- outcomes %>% filter(drug == 1, tx_line == 1)

cat("First-line EVP patients:", nrow(fl_evp), "\n")


# =============================================================================
# 2. Load ctDNA data
# =============================================================================

ctdna_raw <- read_csv(ctdna_path) %>%
  left_join(fl_evp %>% select(record_id, evp_c1d1), by = "record_id") %>%
  filter(record_id %in% fl_evp$record_id) %>%
  mutate(
    ctdna_date    = as.Date(ctdna_date, format = "%m/%d/%y"),
    evp_c1d1      = as.Date(evp_c1d1,  format = "%m/%d/%y"),
    time_from_evp = as.numeric(ctdna_date - evp_c1d1) / 30.44
  ) %>%
  filter(!is.na(redcap_repeat_instrument))


# =============================================================================
# 3. Define analytic cohort (Figures 3 & 4)
#    - ≥1 post-EVP draw within 3 months
#    - Non-zero baseline ctDNA (last pre-EVP draw)
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

baseline_per_pt <- ctdna_raw %>%
  filter(time_from_evp <= 0) %>%
  group_by(record_id) %>%
  filter(time_from_evp == max(time_from_evp)) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, baseline_ctdna = ctdna_value, baseline_time = time_from_evp)

analytic_ids <- baseline_per_pt %>%
  filter(!is.na(baseline_ctdna), baseline_ctdna != 0) %>%
  semi_join(within_3mo_ids, by = "record_id") %>%
  pull(record_id)

n_cohort <- length(analytic_ids)
cat("Analytic cohort n =", n_cohort, "\n\n")

# All ctDNA draws for cohort patients, sorted by patient and time
cohort_draws <- ctdna_raw %>%
  filter(record_id %in% analytic_ids) %>%
  arrange(record_id, time_from_evp)


# =============================================================================
# 4. Panel A — Baseline ctDNA values (log10 scale)
#    One value per patient: the last pre-EVP draw (cohort definition)
# =============================================================================

panel_a_data <- baseline_per_pt %>%
  filter(record_id %in% analytic_ids) %>%
  mutate(log10_val = log10(baseline_ctdna))

med_baseline <- median(panel_a_data$baseline_ctdna, na.rm = TRUE)
q1_baseline  <- quantile(panel_a_data$baseline_ctdna, 0.25, na.rm = TRUE)
q3_baseline  <- quantile(panel_a_data$baseline_ctdna, 0.75, na.rm = TRUE)

panel_a_plot <- ggplot(panel_a_data, aes(x = log10_val)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white",
                 boundary = 0) +
  geom_vline(xintercept = log10(med_baseline), linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text",
           x     = log10(med_baseline), y = Inf,
           label = paste0("Median: ", scales::comma(round(med_baseline, 0)),
                          " [IQR: ", scales::comma(round(q1_baseline, 0)),
                          "\u2013", scales::comma(round(q3_baseline, 0)), "]"),
           hjust = -0.07, vjust = 1.6,
           size  = 2.8, fontface = "italic", color = "firebrick") +
  scale_x_continuous(
    name   = "Baseline ctDNA (copies/mL, log\u2081\u2080 scale)",
    labels = function(x) scales::comma(10^x)
  ) +
  scale_y_continuous(
    name   = "Number of Patients",
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

panel_a <- add_panel_title(
  panel_a_plot,
  paste0("A   Baseline ctDNA Values (n=", n_cohort, ")")
)


# =============================================================================
# 5. Panel B — Number of ctDNA draws per patient (all draws)
# =============================================================================

draws_per_pt <- cohort_draws %>%
  count(record_id, name = "n_draws")

med_draws <- median(draws_per_pt$n_draws)
q1_draws  <- quantile(draws_per_pt$n_draws, 0.25)
q3_draws  <- quantile(draws_per_pt$n_draws, 0.75)
max_draws <- max(draws_per_pt$n_draws)

# Use integer breaks; thin out if max_draws is large
x_breaks_b <- if (max_draws <= 20) {
  seq(1, max_draws, by = 1)
} else {
  seq(0, max_draws, by = 2)
}

panel_b_plot <- ggplot(draws_per_pt, aes(x = n_draws)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white",
                 boundary = 0.5) +
  geom_vline(xintercept = med_draws, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text",
           x     = med_draws, y = Inf,
           label = paste0("Median: ", med_draws,
                          " [IQR: ", q1_draws, "\u2013", q3_draws, "]"),
           hjust = -0.1, vjust = 1.6,
           size  = 2.8, fontface = "italic", color = "firebrick") +
  scale_x_continuous(
    name   = "Number of ctDNA Draws per Patient",
    breaks = x_breaks_b
  ) +
  scale_y_continuous(
    name   = "Number of Patients",
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

panel_b <- add_panel_title(
  panel_b_plot,
  paste0("B   ctDNA Draws per Patient (n=", n_cohort, ")")
)


# =============================================================================
# 6. Panel C — Time intervals between consecutive ctDNA draws
#    Computed across all consecutive draw pairs per patient
# =============================================================================

intervals <- cohort_draws %>%
  group_by(record_id) %>%
  mutate(interval_months = time_from_evp - lag(time_from_evp)) %>%
  ungroup() %>%
  filter(!is.na(interval_months), interval_months > 0)

n_intervals  <- nrow(intervals)
med_interval <- round(median(intervals$interval_months), 2)
q1_interval  <- round(quantile(intervals$interval_months, 0.25), 2)
q3_interval  <- round(quantile(intervals$interval_months, 0.75), 2)
max_interval <- ceiling(max(intervals$interval_months))

x_breaks_c <- seq(0, max_interval, by = if (max_interval <= 12) 2 else 6)

panel_c_plot <- ggplot(intervals, aes(x = interval_months)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white",
                 boundary = 0) +
  geom_vline(xintercept = med_interval, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text",
           x     = med_interval, y = Inf,
           label = paste0("Median: ", med_interval,
                          " months [IQR: ", q1_interval,
                          "\u2013", q3_interval, "]"),
           hjust = -0.07, vjust = 1.6,
           size  = 2.8, fontface = "italic", color = "firebrick") +
  scale_x_continuous(
    name   = "Time Between Consecutive ctDNA Draws (months)",
    breaks = x_breaks_c,
    limits = c(0, max_interval)
  ) +
  scale_y_continuous(
    name   = "Number of Intervals",
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

panel_c <- add_panel_title(
  panel_c_plot,
  paste0("C   Inter-Draw Interval (", n_intervals, " consecutive pairs)")
)


# =============================================================================
# 7. Panel D — Time from EVP initiation to first post-EVP ctDNA draw
# =============================================================================

first_post_evp <- cohort_draws %>%
  filter(time_from_evp > 0) %>%
  group_by(record_id) %>%
  filter(time_from_evp == min(time_from_evp)) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, time_to_first = time_from_evp)

med_first <- round(median(first_post_evp$time_to_first), 2)
q1_first  <- round(quantile(first_post_evp$time_to_first, 0.25), 2)
q3_first  <- round(quantile(first_post_evp$time_to_first, 0.75), 2)

panel_d_plot <- ggplot(first_post_evp, aes(x = time_to_first)) +
  geom_histogram(binwidth = 0.25, fill = "steelblue", color = "white",
                 boundary = 0) +
  geom_vline(xintercept = med_first, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text",
           x     = med_first, y = Inf,
           label = paste0("Median: ", med_first,
                          " months [IQR: ", q1_first,
                          "\u2013", q3_first, "]"),
           hjust = -0.07, vjust = 1.6,
           size  = 2.8, fontface = "italic", color = "firebrick") +
  scale_x_continuous(
    name   = "Time from EVP Initiation to First ctDNA Draw (months)",
    breaks = seq(0, 3, by = 0.5),
    limits = c(0, 3.1)
  ) +
  scale_y_continuous(
    name   = "Number of Patients",
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

panel_d <- add_panel_title(
  panel_d_plot,
  paste0("D   Time to First Post-EVP ctDNA Draw (n=", nrow(first_post_evp), ")")
)


# =============================================================================
# 8. Assemble Supplementary Figure 3
# =============================================================================

supp_fig2 <- plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2
)

print(supp_fig2)


# =============================================================================
# 9. Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Supplementary_Figure3_final.pdf"),
  plot     = supp_fig2,
  width    = 10, height = 8,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Supplementary_Figure3_final.png"),
  plot     = supp_fig2,
  width    = 10, height = 8,
  dpi      = 300
)

cat("Saved to:", output_dir, "\n")
