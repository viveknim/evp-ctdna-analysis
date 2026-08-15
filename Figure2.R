# ============================================================================
# Figure 2 — Publication Quality Combined Figure
# EVP ctDNA Project
#
# 3×2 grid layout:
# Row 1 (OS):  A: 1-mo Landmark OS  | B: 2-mo Landmark OS  | C: 3-mo Landmark OS
# Row 2 (PFS): D: 1-mo Landmark PFS | E: 2-mo Landmark PFS | F: 3-mo Landmark PFS
#
# Target size: 8.5 × 10 inches
# ============================================================================

library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(lubridate)
library(ggplot2)
library(cowplot)
library(gt)

# ---- USER: Update these paths -----------------------------------------------
data_dir      <- "data"
outcomes_path <- file.path(data_dir, "ev_ctdna_survival_tab.csv")
ctdna_path    <- file.path(data_dir, "ctdna_data.csv")
output_dir    <- "output"
# -----------------------------------------------------------------------------


# =============================================================================
# Shared publication theme
# =============================================================================
#
# NOTE ON TITLE ALIGNMENT
# All panels use pure ggplot2 KM curves (no ggsurvplot), so titles are stripped
# from every ggplot and added as fixed-height strips via add_panel_title().
# All 6 panels share the same absolute title y-position in the assembled figure.
# =============================================================================

theme_pub <- function() {
  theme_classic(base_size = 7) +
    theme(
      plot.title      = element_text(face = "bold", size = 7,   hjust = 0),
      plot.subtitle   = element_text(size = 5.5, color = "gray40", hjust = 0),
      axis.title      = element_text(face = "bold", size = 6.5),
      axis.text       = element_text(size = 6,   color = "black"),
      legend.title    = element_text(face = "bold", size = 6),
      legend.text     = element_text(size = 5.5),
      legend.key.size = unit(0.28, "cm"),
      plot.margin     = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt")
    )
}

add_panel_title <- function(panel_grid, title_text, subtitle_text = NULL) {
  t_strip <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 7,
               x = 0.02, hjust = 0, vjust = 0.5)

  if (!is.null(subtitle_text) && nchar(subtitle_text) > 0) {
    s_strip <- ggdraw() +
      draw_label(subtitle_text, color = "gray40", size = 5,
                 x = 0.02, hjust = 0, vjust = 0.5)
    plot_grid(t_strip, s_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 0.75, 9))
  } else {
    plot_grid(t_strip, panel_grid,
              ncol = 1, rel_heights = c(1, 10))
  }
}


# =============================================================================
# 1. Load outcomes
# =============================================================================

outcomes <- read_csv(outcomes_path) %>%
  mutate(
    evp_c1d1      = as.Date(evp_c1d1,      format = "%m/%d/%y"),
    os_event_date = as.Date(os_event_date,  format = "%m/%d/%y"),
    pfs_date      = as.Date(pfs_date,       format = "%m/%d/%y"),
    os_months     = as.numeric(os_event_date - evp_c1d1) / 30.44,
    os_status     = ifelse(os_event == 2, 1, 0),
    pfs_months    = as.numeric(pfs_date - evp_c1d1) / 30.44,
    pfs_status    = ifelse(pfs_event == 1 | pfs_event == 2, 1, 0)
  )


# =============================================================================
# 2. Load ctDNA data
# =============================================================================

ctdna_raw <- read_csv(ctdna_path) %>%
  left_join(outcomes %>% select(record_id, evp_c1d1), by = "record_id") %>%
  mutate(
    ctdna_date    = as.Date(ctdna_date, format = "%m/%d/%y"),
    evp_c1d1      = as.Date(evp_c1d1,  format = "%m/%d/%y"),
    time_from_evp = as.numeric(ctdna_date - evp_c1d1) / 30.44
  ) %>%
  filter(!is.na(redcap_repeat_instrument))


# =============================================================================
# 3. Define cohort: first-line EVP, ≥1 ctDNA draw within 3 months
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

cohort <- outcomes %>%
  filter(drug == 1, tx_line == 1) %>%
  semi_join(within_3mo_ids, by = "record_id")

cat("Cohort n =", nrow(cohort), "\n")


# =============================================================================
# 4. First clearance per patient
# =============================================================================

first_clearance <- ctdna_raw %>%
  filter(record_id %in% cohort$record_id,
         time_from_evp > 0, ctdna_value == 0) %>%
  group_by(record_id) %>%
  arrange(time_from_evp) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, time_to_clearance = time_from_evp)

cat("Cleared:", nrow(first_clearance),
    "| Never cleared:", nrow(cohort) - nrow(first_clearance), "\n\n")


# =============================================================================
# 5. Landmark KM panel helper (pure ggplot2 — no ggsurvplot theme conflicts)
#
#    Landmark analysis at time t*:
#      - ≥1 ctDNA draw within landmark_t months (clearance assessable)
#      - pfs_months > landmark_t used for BOTH OS and PFS analyses so the
#        same patients (alive AND progression-free at landmark) appear in both
#      - Time axis starts at 0 from the landmark (outcome - t*)
#      - Risk table built from summary(km_fit) — full theme control
#      - Log-rank p-value displayed in upper-right of KM panel
# =============================================================================

make_lm_panel <- function(landmark_t, time_col, status_col,
                           outcome_short, panel_label) {

  # Patients with ≥1 ctDNA draw within landmark_t months of EVP
  within_lm_ids <- ctdna_raw %>%
    filter(record_id %in% cohort$record_id,
           time_from_evp > 0, time_from_evp <= landmark_t) %>%
    distinct(record_id)

  # Filter by pfs_months > landmark_t for BOTH OS and PFS panels so cohort N
  # matches (patients must be alive AND progression-free at the landmark time)
  lm_data <- cohort %>%
    semi_join(within_lm_ids, by = "record_id") %>%
    mutate(
      time_outcome = .data[[time_col]],
      event_status = .data[[status_col]]
    ) %>%
    filter(pfs_months > landmark_t) %>%
    left_join(
      first_clearance %>%
        filter(time_to_clearance <= landmark_t) %>%
        transmute(record_id, cleared_by_lm = TRUE),
      by = "record_id"
    ) %>%
    mutate(
      cleared_by_lm   = ifelse(is.na(cleared_by_lm), FALSE, cleared_by_lm),
      clearance_group = factor(
        ifelse(cleared_by_lm, "Cleared", "Not Cleared"),
        levels = c("Not Cleared", "Cleared")
      ),
      time_from_lm = time_outcome - landmark_t
    )

  n_cleared <- sum(lm_data$cleared_by_lm)
  n_noclear <- sum(!lm_data$cleared_by_lm)
  n_total   <- nrow(lm_data)

  cat(panel_label, "| landmark =", landmark_t, "mo |",
      outcome_short, "| N =", n_total,
      "(Cleared =", n_cleared, ", Not Cleared =", n_noclear, ")\n")

  km_fit <- survfit(Surv(time_from_lm, event_status) ~ clearance_group,
                    data = lm_data)

  lrt  <- survdiff(Surv(time_from_lm, event_status) ~ clearance_group,
                   data = lm_data)
  pval <- 1 - pchisq(lrt$chisq, df = 1)
  plab <- ifelse(pval < 0.001, "p<0.001",
                 paste0("p=", sprintf("%.3f", pval)))

  # ---- KM curve (pure ggplot2 via broom::tidy) ------------------------------
  km_df <- broom::tidy(km_fit) %>%
    mutate(
      strata = gsub("clearance_group=", "", strata),
      strata = factor(strata, levels = c("Not Cleared", "Cleared"))
    )

  x_max    <- ceiling(max(lm_data$time_from_lm, na.rm = TRUE))
  x_breaks <- seq(0, x_max, by = max(1, round(x_max / 4)))

  km_plot <- ggplot(km_df, aes(x = time, y = estimate,
                               color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.15, color = NA) +
    geom_step(linewidth = 0.6) +
    geom_point(
      data  = km_df %>% filter(n.censor > 0),
      shape = 3, size = 1.0
    ) +
    annotate("text", x = Inf, y = Inf, label = plab,
             hjust = 1.05, vjust = 1.8, size = 2.2, fontface = "italic") +
    scale_color_manual(
      values = c("Not Cleared" = "#B2182B", "Cleared" = "#2166AC"),
      labels = c(paste0("Not Cleared (n=", n_noclear, ")"),
                 paste0("Cleared (n=",     n_cleared, ")")),
      name   = "ctDNA at Landmark"
    ) +
    scale_fill_manual(
      values = c("Not Cleared" = "#B2182B", "Cleared" = "#2166AC"),
      guide  = "none"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      name   = paste0(outcome_short, " Probability")
    ) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = x_breaks,
      name   = paste0("Months from ", landmark_t, "-Month Landmark")
    ) +
    labs(title = NULL, subtitle = NULL) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = 6),
      legend.text     = element_text(size = 5),
      legend.key.size = unit(0.28, "cm"),
      plot.margin     = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt")
    )

  # ---- Risk table (pure ggplot2) --------------------------------------------
  risk_sum <- summary(km_fit, times = x_breaks, extend = TRUE)
  risk_df  <- data.frame(
    time   = risk_sum$time,
    n.risk = risk_sum$n.risk,
    strata = factor(
      gsub("clearance_group=", "", as.character(risk_sum$strata)),
      levels = c("Cleared", "Not Cleared")   # reversed so Not Cleared plots on top
    )
  )

  risk_tbl <- ggplot(risk_df, aes(x = time, y = strata, label = n.risk)) +
    geom_text(size = 1.8, color = "black") +
    scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
    scale_y_discrete(labels = c("Not Cleared" = "Not Clr",
                                 "Cleared"     = "Cleared")) +
    labs(x = NULL, y = NULL, title = "No. at risk") +
    theme_pub() +
    theme(
      plot.title         = element_text(size = 5, face = "bold", hjust = 0),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.line.x        = element_blank(),
      axis.text.y        = element_text(size = 5, color = "black"),
      panel.grid.major   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(t = 0, r = 4, b = 4, l = 4, unit = "pt")
    )

  # Inner grid: KM curve (top) + risk table (bottom)
  panel_raw <- plot_grid(
    km_plot, risk_tbl,
    ncol        = 1,
    rel_heights = c(4, 1)
  )

  # Attach external title/subtitle strips at fixed height
  add_panel_title(
    panel_raw,
    paste0(panel_label, "   ", landmark_t, "-Month Landmark: ", outcome_short)
  )
}

# OS landmark panels
panel_a <- make_lm_panel(1, "os_months",  "os_status",  "OS",  "A")
panel_b <- make_lm_panel(2, "os_months",  "os_status",  "OS",  "B")
panel_c <- make_lm_panel(3, "os_months",  "os_status",  "OS",  "C")

# PFS landmark panels
panel_d <- make_lm_panel(1, "pfs_months", "pfs_status", "PFS", "D")
panel_e <- make_lm_panel(2, "pfs_months", "pfs_status", "PFS", "E")
panel_f <- make_lm_panel(3, "pfs_months", "pfs_status", "PFS", "F")


# =============================================================================
# 6. Assemble Figure 2 (3×2 grid)
#
#     Row 1: A (1-mo OS)  | B (2-mo OS)  | C (3-mo OS)
#     Row 2: D (1-mo PFS) | E (2-mo PFS) | F (3-mo PFS)
# =============================================================================

row1 <- plot_grid(panel_a, panel_b, panel_c, nrow = 1)
row2 <- plot_grid(panel_d, panel_e, panel_f, nrow = 1)

figure_2 <- plot_grid(row1, row2, ncol = 1)

print(figure_2)


# =============================================================================
# 7. Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Figure2_publication_final.pdf"),
  plot     = figure_2,
  width    = 8.5, height = 10,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Figure2_publication_final.png"),
  plot     = figure_2,
  width    = 8.5, height = 10,
  dpi      = 300
)

cat("Saved to:", output_dir, "\n")


# =============================================================================
# 8. Median OS and PFS by clearance group for each landmark analysis
# =============================================================================

get_lm_medians <- function(landmark_t, time_col, status_col, outcome_short) {

  within_lm_ids <- ctdna_raw %>%
    filter(record_id %in% cohort$record_id,
           time_from_evp > 0, time_from_evp <= landmark_t) %>%
    distinct(record_id)

  lm_data <- cohort %>%
    semi_join(within_lm_ids, by = "record_id") %>%
    mutate(
      time_outcome = .data[[time_col]],
      event_status = .data[[status_col]]
    ) %>%
    filter(pfs_months > landmark_t) %>%
    left_join(
      first_clearance %>%
        filter(time_to_clearance <= landmark_t) %>%
        transmute(record_id, cleared_by_lm = TRUE),
      by = "record_id"
    ) %>%
    mutate(
      cleared_by_lm   = ifelse(is.na(cleared_by_lm), FALSE, cleared_by_lm),
      clearance_group = factor(
        ifelse(cleared_by_lm, "Cleared", "Not Cleared"),
        levels = c("Not Cleared", "Cleared")
      ),
      time_from_lm = time_outcome - landmark_t
    )

  km_fit <- survfit(Surv(time_from_lm, event_status) ~ clearance_group,
                    data = lm_data)

  tbl <- summary(km_fit)$table

  fmt_mo <- function(x) ifelse(is.na(x), "NR", sprintf("%.1f", x))

  data.frame(
    Landmark  = paste0(landmark_t, "-Month"),
    Outcome   = outcome_short,
    Group     = gsub("clearance_group=", "", rownames(tbl)),
    N         = as.integer(tbl[, "records"]),
    Events    = as.integer(tbl[, "events"]),
    Median_mo = fmt_mo(tbl[, "median"]),
    CI_95     = paste0(fmt_mo(tbl[, "0.95LCL"]), "\u2013", fmt_mo(tbl[, "0.95UCL"])),
    stringsAsFactors = FALSE
  )
}

median_rows <- bind_rows(
  get_lm_medians(1, "os_months",  "os_status",  "OS"),
  get_lm_medians(2, "os_months",  "os_status",  "OS"),
  get_lm_medians(3, "os_months",  "os_status",  "OS"),
  get_lm_medians(1, "pfs_months", "pfs_status", "PFS"),
  get_lm_medians(2, "pfs_months", "pfs_status", "PFS"),
  get_lm_medians(3, "pfs_months", "pfs_status", "PFS")
)

cat("\nMedian survival by landmark and clearance group:\n")
print(median_rows)

median_tbl <- median_rows %>%
  gt() %>%
  tab_header(
    title    = md("**Median Survival by Landmark Time Point and ctDNA Clearance Status**"),
    subtitle = md("Months from landmark | NR = not reached | 95% CI shown")
  ) %>%
  cols_label(
    Landmark  = "Landmark",
    Outcome   = "Outcome",
    Group     = "ctDNA Status",
    N         = "N",
    Events    = "Events",
    Median_mo = "Median (months)",
    CI_95     = "95% CI"
  ) %>%
  cols_align(align = "center", columns = c(N, Events, Median_mo, CI_95)) %>%
  cols_align(align = "left",   columns = c(Landmark, Outcome, Group)) %>%
  tab_row_group(label = "3-Month Landmark", rows = Landmark == "3-Month") %>%
  tab_row_group(label = "2-Month Landmark", rows = Landmark == "2-Month") %>%
  tab_row_group(label = "1-Month Landmark", rows = Landmark == "1-Month") %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  opt_stylize(style = 1, color = "gray") %>%
  tab_options(table.font.size = px(11))

median_tbl

gtsave(median_tbl,
       file.path(output_dir, "Figure2_median_survival_table.docx"))
gtsave(median_tbl,
       file.path(output_dir, "Figure2_median_survival_table.png"),
       vwidth = 900, zoom = 2)

cat("Median survival table saved to:", output_dir, "\n")
