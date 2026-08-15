# ============================================================================
# Figure 3 (Modified) — ctDNA Dynamics After EVP Initiation
# EVP ctDNA Project
#
# Cohort: First-line EVP patients with non-zero baseline ctDNA AND
#         ≥1 post-EVP draw within 3 months of EVP initiation
#
# Panel A: Overall spaghetti — all patients, 0–12 months
# Panel B: Nadir ctDNA distribution (5-bucket bar chart)
# Panel C: First post-EVP ctDNA (% baseline) vs time to first draw,
#           colored by eventual clearance status
# Panel D: Time-to-ctDNA-clearance curve + risk table (clearers only)
#
# Layout (10 × 6 inches):
#   Row 1 : [A: spaghetti, 1.5×] [B: nadir bar, 1×]
#   Row 2 : [C: scatter, 1×]     [D: TTE curve, 1×]
# ============================================================================

library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(lubridate)
library(survival)
library(broom)
library(cowplot)

# ---- USER: Update these paths -----------------------------------------------
data_dir      <- "data"
outcomes_path <- file.path(data_dir, "ev_ctdna_survival_tab.csv")
ctdna_path    <- file.path(data_dir, "ctdna_data.csv")
output_dir    <- "output"
# -----------------------------------------------------------------------------


# =============================================================================
# Shared theme and color palettes
# =============================================================================

theme_pub <- function() {
  theme_bw(base_size = 8) +
    theme(
      plot.title       = element_text(face = "bold", size = 8,   hjust = 0),
      plot.subtitle    = element_text(size = 5.5, color = "gray40", hjust = 0),
      axis.title       = element_text(face = "bold", size = 7),
      axis.text        = element_text(size = 6,   color = "black"),
      legend.title     = element_text(face = "bold", size = 6.5),
      legend.text      = element_text(size = 6),
      legend.key.size  = unit(0.28, "cm"),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 3, r = 4, b = 3, l = 4, unit = "pt")
    )
}

nadir_colors <- c(
  "Clearance (0%)"           = "#2166AC",
  ">0\u201310% of Baseline"  = "#74C476",
  "11\u201350% of Baseline"  = "#F59E0B",
  "51\u2013100% of Baseline" = "#F97316",
  ">100% of Baseline"        = "#B2182B"
)
nadir_levels <- names(nadir_colors)

# "Ever Cleared" = blue; "Never Cleared" = purple
clear_palette <- c("Ever Cleared" = "#2166AC", "Never Cleared" = "#762A83")


# =============================================================================
# External title helper
# =============================================================================

add_panel_title <- function(panel_grid, title_text, subtitle_text = NULL) {
  t_strip <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 8,
               x = 0.02, hjust = 0, vjust = 0.5)
  if (!is.null(subtitle_text) && nchar(subtitle_text) > 0) {
    s_strip <- ggdraw() +
      draw_label(subtitle_text, color = "gray40", size = 5.5,
                 x = 0.02, hjust = 0, vjust = 0.5)
    plot_grid(t_strip, s_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 0.75, 9))
  } else {
    plot_grid(t_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 10))
  }
}


# =============================================================================
# 1. Load outcomes and define first-line EVP cohort
# =============================================================================

outcomes <- read_csv(outcomes_path) %>%
  mutate(evp_c1d1 = as.Date(evp_c1d1, format = "%m/%d/%y"))

cohort_base <- outcomes %>%
  filter(drug == 1, tx_line == 1)

cat("First-line EVP patients:", nrow(cohort_base), "\n")


# =============================================================================
# 2. Load ctDNA data
# =============================================================================

ctdna_raw <- read_csv(ctdna_path) %>%
  left_join(cohort_base %>% select(record_id, evp_c1d1), by = "record_id") %>%
  filter(record_id %in% cohort_base$record_id) %>%
  mutate(
    ctdna_date    = as.Date(ctdna_date, format = "%m/%d/%y"),
    evp_c1d1      = as.Date(evp_c1d1,  format = "%m/%d/%y"),
    time_from_evp = as.numeric(ctdna_date - evp_c1d1) / 30.44
  ) %>%
  filter(!is.na(redcap_repeat_instrument))


# =============================================================================
# 3. Define analytic cohort
#    - ≥1 post-EVP draw within 3 months
#    - Baseline ctDNA available (last pre-EVP draw) and non-zero
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

baseline_info <- ctdna_raw %>%
  filter(time_from_evp <= 0) %>%
  group_by(record_id) %>%
  filter(time_from_evp == max(time_from_evp)) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, baseline_value = ctdna_value, baseline_time = time_from_evp)

analytic_ids <- baseline_info %>%
  filter(!is.na(baseline_value), baseline_value != 0) %>%
  semi_join(within_3mo_ids, by = "record_id") %>%
  pull(record_id)

n_cohort <- length(analytic_ids)
cat("Analytic cohort n =", n_cohort, "\n\n")


# =============================================================================
# 4. Compute ctDNA % of baseline for all draws (from baseline draw onwards)
# =============================================================================

plot_data_all <- ctdna_raw %>%
  filter(record_id %in% analytic_ids) %>%
  group_by(record_id) %>%
  mutate(
    baseline_ctdna_val = {
      idx <- which(time_from_evp <= 0)
      if (length(idx) == 0) NA_real_
      else ctdna_value[idx][which.max(time_from_evp[idx])]
    },
    ctdna_pct_baseline = ifelse(
      is.na(baseline_ctdna_val) | baseline_ctdna_val == 0,
      NA_real_,
      100 * ctdna_value / baseline_ctdna_val
    ),
    baseline_time_val = {
      pre <- time_from_evp[time_from_evp <= 0]
      if (length(pre) == 0) NA_real_ else max(pre)
    }
  ) %>%
  filter(!is.na(baseline_time_val),
         time_from_evp >= baseline_time_val) %>%
  ungroup()

x_min_spag <- floor(min(plot_data_all$baseline_time_val, na.rm = TRUE) * 2) / 2
x_max_spag <- 12


# =============================================================================
# 5. Nadir bucket per patient
# =============================================================================

nadir_per_pt <- plot_data_all %>%
  filter(time_from_evp > 0) %>%
  group_by(record_id) %>%
  summarise(
    min_pct      = min(ctdna_pct_baseline, na.rm = TRUE),
    reached_zero = any(ctdna_value == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    nadir_bucket = factor(case_when(
      reached_zero   ~ "Clearance (0%)",
      min_pct <= 10  ~ ">0\u201310% of Baseline",
      min_pct <= 50  ~ "11\u201350% of Baseline",
      min_pct <= 100 ~ "51\u2013100% of Baseline",
      TRUE           ~ ">100% of Baseline"
    ), levels = nadir_levels),
    cleared = reached_zero
  )

plot_data_all <- plot_data_all %>%
  left_join(nadir_per_pt %>% select(record_id, nadir_bucket, cleared),
            by = "record_id")

n_spaghetti <- length(unique(plot_data_all$record_id))
cat("Spaghetti cohort n =", n_spaghetti, "\n")
print(nadir_per_pt %>% count(nadir_bucket))
cat("\n")


# =============================================================================
# 6. First post-EVP draw summary (Panel C)
# =============================================================================

first_draw <- ctdna_raw %>%
  filter(record_id %in% analytic_ids, time_from_evp > 0) %>%
  group_by(record_id) %>%
  filter(time_from_evp == min(time_from_evp)) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(baseline_info, by = "record_id") %>%
  mutate(
    first_ctdna_pct_baseline = 100 * ctdna_value / baseline_value,
    time_to_first            = time_from_evp
  ) %>%
  left_join(nadir_per_pt %>% select(record_id, cleared), by = "record_id") %>%
  mutate(
    clearance_group = factor(
      ifelse(cleared, "Ever Cleared", "Never Cleared"),
      levels = c("Ever Cleared", "Never Cleared")
    )
  )

n_first_draw    <- nrow(first_draw)
n_ever_cleared  <- sum(first_draw$clearance_group == "Ever Cleared",  na.rm = TRUE)
n_never_cleared <- sum(first_draw$clearance_group == "Never Cleared", na.rm = TRUE)

cat("First-draw N =", n_first_draw,
    "| Ever Cleared:", n_ever_cleared,
    "| Never Cleared:", n_never_cleared, "\n\n")


# =============================================================================
# 7. Panel A — Overall spaghetti (all patients, 0–12 months)
# =============================================================================

panel_a_plot <- ggplot(
  plot_data_all,
  aes(x = time_from_evp, y = ctdna_pct_baseline, group = record_id)
) +
  geom_line(alpha = 0.28, color = "steelblue", linewidth = 0.38) +
  geom_point(alpha = 0.35, size = 0.5, color = "steelblue") +
  geom_smooth(
    data    = subset(plot_data_all, time_from_evp <= x_max_spag),
    aes(group = 1),
    method  = "loess", formula = y ~ x,
    color   = "firebrick", fill = "mistyrose",
    se      = TRUE, linewidth = 0.75, alpha = 0.25
  ) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black",
             alpha = 0.55, linewidth = 0.35) +
  coord_cartesian(xlim = c(x_min_spag, x_max_spag)) +
  scale_x_continuous(name   = "Months from EVP Initiation",
                     breaks = seq(0, x_max_spag, by = 2)) +
  scale_y_continuous(
    name   = "ctDNA (% of Baseline)",
    trans  = "pseudo_log",
    breaks = c(0, 1, 10, 100, 1000),
    labels = label_percent(scale = 1)
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub()

panel_a <- add_panel_title(
  panel_a_plot,
  paste0("A   ctDNA Trajectories Post-EVP (n=", n_spaghetti, ")"),
  "Steelblue: individual patients | Red line: LOESS (95% CI) | Dashed: baseline"
)


# =============================================================================
# 8. Panel B — Nadir ctDNA 5-bucket bar chart
# =============================================================================

nadir_summary <- nadir_per_pt %>%
  count(nadir_bucket) %>%
  mutate(pct = n / sum(n))

panel_b_plot <- ggplot(nadir_summary,
                       aes(x = nadir_bucket, y = n, fill = nadir_bucket)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.72) +
  geom_text(
    aes(label = paste0("n=", n, "\n(",
                       scales::percent(pct, accuracy = 1), ")")),
    vjust = -0.3, size = 2.2, fontface = "bold"
  ) +
  scale_fill_manual(values = nadir_colors, guide = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(name   = "Number of Patients",
                     expand = expansion(mult = c(0, 0.22))) +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme_pub() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(size = 5.5, lineheight = 1.1)
  )

panel_b <- add_panel_title(
  panel_b_plot,
  paste0("B   Nadir ctDNA Distribution (n=", n_spaghetti, ")"))


# =============================================================================
# 9. Panel C — Scatter: first ctDNA (% baseline) vs time to first draw
#               colored by eventual clearance status
#               "Never Cleared" shown in purple
# =============================================================================

panel_c_plot <- ggplot(first_draw,
                       aes(x = time_to_first, y = first_ctdna_pct_baseline,
                           color = clearance_group)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black",
             alpha = 0.55, linewidth = 0.35) +
  geom_point(size = 1.6, alpha = 0.75) +
  scale_color_manual(
    values = clear_palette,
    labels = c(paste0("Ever Cleared (n=",  n_ever_cleared,  ")"),
               paste0("Never Cleared (n=", n_never_cleared, ")")),
    name   = "ctDNA Clearance"
  ) +
  scale_x_continuous(
    name   = "Time to First Post-EVP ctDNA Draw (months)",
    breaks = seq(0, 3, by = 0.5),
    limits = c(0, NA)
  ) +
  scale_y_continuous(
    name   = "First Post-EVP ctDNA (% of Baseline)",
    trans  = "pseudo_log",
    breaks = c(0, 1, 10, 100, 1000),
    labels = label_percent(scale = 1)
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub() +
  theme(
    legend.position  = "bottom",
    legend.spacing.x = unit(0.08, "cm")
  )

panel_c <- add_panel_title(
  plot_grid(panel_c_plot, ggdraw(),
            ncol = 1, rel_heights = c(4, 1)),
  paste0("C   First ctDNA vs Time to Draw (n=", n_first_draw, ")"))


# =============================================================================
# 10. Panel D — Time-to-ctDNA-clearance curve (pure ggplot2, clearers only)
#               Risk table: patients yet to clear at each time point
# =============================================================================

clearance_times <- ctdna_raw %>%
  filter(record_id %in% analytic_ids,
         time_from_evp > 0, ctdna_value == 0) %>%
  group_by(record_id) %>%
  summarise(time_to_clear = min(time_from_evp), cleared = 1L, .groups = "drop")

n_cleared_tte <- nrow(clearance_times)
med_clear     <- round(median(clearance_times$time_to_clear), 2)
q1_clear      <- round(quantile(clearance_times$time_to_clear, 0.25), 2)
q3_clear      <- round(quantile(clearance_times$time_to_clear, 0.75), 2)

cat("TTE clearance: N =", n_cleared_tte,
    "| Median:", med_clear, "months",
    "| IQR:", q1_clear, "-", q3_clear, "months\n\n")

fit_clear  <- survfit(Surv(time_to_clear, cleared) ~ 1, data = clearance_times)
x_max_tte  <- ceiling(max(clearance_times$time_to_clear) + 0.5)
x_breaks_d <- seq(0, x_max_tte, by = max(1, round(x_max_tte / 5)))

tte_df <- broom::tidy(fit_clear) %>%
  mutate(
    cum_cleared = 1 - estimate,
    cum_lo      = 1 - conf.high,
    cum_hi      = 1 - conf.low
  )

tte_curve <- ggplot(tte_df, aes(x = time, y = cum_cleared)) +
  geom_ribbon(aes(ymin = cum_lo, ymax = cum_hi),
              fill = "#2166AC", alpha = 0.15, color = NA) +
  geom_step(color = "#2166AC", linewidth = 0.65) +
  geom_point(data  = tte_df %>% filter(n.censor > 0),
             aes(x = time, y = cum_cleared),
             shape = 3, size = 1.0, color = "#2166AC") +
  annotate("text", x = Inf, y = -Inf,
           label  = paste0("Median: ", med_clear,
                           " mo [IQR: ", q1_clear, "\u2013", q3_clear, "]"),
           hjust  = 1.05, vjust = -0.4,
           size   = 2.0, fontface = "italic") +
  scale_x_continuous(limits = c(0, x_max_tte), breaks = x_breaks_d,
                     name = "Time from EVP Initiation (months)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
                     name   = "Cumulative Proportion Cleared",
                     labels = label_percent(accuracy = 1)) +
  labs(title = NULL, subtitle = NULL) +
  theme_pub()

tte_risk_sum <- summary(fit_clear, times = x_breaks_d, extend = TRUE)
tte_risk_df  <- data.frame(
  time   = tte_risk_sum$time,
  n.risk = tte_risk_sum$n.risk,
  group  = "Yet to clear"
)

tte_risk_tbl <- ggplot(tte_risk_df, aes(x = time, y = group, label = n.risk)) +
  geom_text(size = 1.8, color = "black") +
  scale_x_continuous(limits = c(0, x_max_tte), breaks = x_breaks_d) +
  scale_y_discrete() +
  labs(x = NULL, y = NULL, title = "No. yet to clear") +
  theme_pub() +
  theme(
    plot.title   = element_text(size = 5, face = "bold", hjust = 0),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y  = element_text(size = 4.5, color = "gray40"),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    plot.margin  = margin(t = 0, r = 4, b = 4, l = 4, unit = "pt")
  )

panel_d <- add_panel_title(
  plot_grid(tte_curve, tte_risk_tbl,
            ncol = 1, rel_heights = c(4, 1)),
  paste0("D   Time to ctDNA Clearance (n=", n_cleared_tte, ")"))


# =============================================================================
# 11. Assemble Figure 3 (Modified)
# =============================================================================

row1 <- plot_grid(
  panel_a, panel_b,
  nrow       = 1,
  rel_widths = c(1.5, 1)
)

row2 <- plot_grid(
  panel_c, panel_d,
  nrow = 1
)

figure_3_mod <- plot_grid(
  row1, row2,
  ncol        = 1,
  rel_heights = c(1, 1)
)

print(figure_3_mod)


# =============================================================================
# 12. Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Figure3_publication_modified.pdf"),
  plot     = figure_3_mod,
  width    = 10, height = 6,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Figure3_publication_modified.png"),
  plot     = figure_3_mod,
  width    = 10, height = 6,
  dpi      = 300
)

cat("Saved to:", output_dir, "\n")


# =============================================================================
# 13. Post-clearance ctDNA status tabulation
#     Among patients who achieved ctDNA clearance (N = n_cleared_tte):
#       - Sustained clearance: all draws after clearance remain undetectable
#       - Rebound:             ≥1 post-clearance draw with detectable ctDNA
#       - No further follow-up: no ctDNA draws recorded after clearance event
# =============================================================================

post_clearance_draws <- ctdna_raw %>%
  filter(record_id %in% clearance_times$record_id) %>%
  left_join(clearance_times %>% select(record_id, time_to_clear),
            by = "record_id") %>%
  filter(time_from_evp > time_to_clear)

post_clearance_summary <- clearance_times %>%
  select(record_id, time_to_clear) %>%
  left_join(
    post_clearance_draws %>%
      group_by(record_id) %>%
      summarise(
        n_post_draws   = n(),
        any_detectable = any(ctdna_value > 0, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "record_id"
  ) %>%
  mutate(
    post_clearance_status = case_when(
      is.na(n_post_draws)  ~ "No further follow-up after clearance",
      any_detectable       ~ "Rebound (ctDNA re-emerged after clearance)",
      TRUE                 ~ "Sustained clearance (remained undetectable)"
    )
  )

cat("\n=== Post-clearance ctDNA status (N =", n_cleared_tte, "who cleared) ===\n")
post_tab <- post_clearance_summary %>%
  count(post_clearance_status) %>%
  mutate(pct = paste0(round(100 * n / n_cleared_tte, 1), "%")) %>%
  rename(Status = post_clearance_status, N = n, `%` = pct)
print(post_tab, row.names = FALSE)
cat("\n")
