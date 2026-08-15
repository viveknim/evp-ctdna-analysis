# ============================================================================
# Supplementary Figure — 6-Month Landmark Analysis
# EVP ctDNA Project
#
# 1×2 layout:
# Panel A: 6-mo Landmark OS  | Panel B: 6-mo Landmark PFS
#
# Target size: 8.5 × 5 inches
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

theme_pub <- function() {
  theme_classic(base_size = 8) +
    theme(
      plot.title      = element_text(face = "bold", size = 9,   hjust = 0),
      plot.subtitle   = element_text(size = 7, color = "gray40", hjust = 0),
      axis.title      = element_text(face = "bold", size = 8),
      axis.text       = element_text(size = 7,   color = "black"),
      legend.title    = element_text(face = "bold", size = 7),
      legend.text     = element_text(size = 6.5),
      legend.key.size = unit(0.35, "cm"),
      plot.margin     = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt")
    )
}

add_panel_title <- function(panel_grid, title_text, subtitle_text = NULL) {
  t_strip <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 9,
               x = 0.02, hjust = 0, vjust = 0.5)
  
  if (!is.null(subtitle_text) && nchar(subtitle_text) > 0) {
    s_strip <- ggdraw() +
      draw_label(subtitle_text, color = "gray40", size = 6.5,
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
#    (Study inclusion criterion — unchanged)
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

cohort <- outcomes %>%
  filter(drug == 1, tx_line == 1) %>%
  semi_join(within_3mo_ids, by = "record_id")

cat("Cohort n =", nrow(cohort), "\n")


# =============================================================================
# 4. First clearance per patient (no time restriction)
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
# 5. Landmark KM panel helper (pure ggplot2)
# =============================================================================

make_lm_panel <- function(landmark_t, time_col, status_col,
                          outcome_short, panel_label) {
  
  
  # Patients with ≥1 ctDNA draw within landmark_t months of EVP
  # (capped at 3 months per study inclusion)
  within_lm_ids <- ctdna_raw %>%
    filter(record_id %in% cohort$record_id,
           time_from_evp > 0, time_from_evp <= min(landmark_t, 3)) %>%
    distinct(record_id)
  
  # Filter by pfs_months > landmark_t for BOTH OS and PFS panels
  # (patients must be alive AND progression-free at the landmark time)
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
  x_breaks <- seq(0, x_max, by = max(1, round(x_max / 5)))
  
  km_plot <- ggplot(km_df, aes(x = time, y = estimate,
                               color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                alpha = 0.15, color = NA) +
    geom_step(linewidth = 0.7) +
    geom_point(
      data  = km_df %>% filter(n.censor > 0),
      shape = 3, size = 1.2
    ) +
    annotate("text", x = Inf, y = Inf, label = plab,
             hjust = 1.05, vjust = 1.8, size = 2.8, fontface = "italic") +
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
      legend.title    = element_text(face = "bold", size = 7),
      legend.text     = element_text(size = 6),
      legend.key.size = unit(0.35, "cm"),
      plot.margin     = margin(t = 4, r = 6, b = 4, l = 6, unit = "pt")
    )
  
  # ---- Risk table (pure ggplot2) --------------------------------------------
  risk_sum <- summary(km_fit, times = x_breaks, extend = TRUE)
  risk_df  <- data.frame(
    time   = risk_sum$time,
    n.risk = risk_sum$n.risk,
    strata = factor(
      gsub("clearance_group=", "", as.character(risk_sum$strata)),
      levels = c("Cleared", "Not Cleared")
    )
  )
  
  risk_tbl <- ggplot(risk_df, aes(x = time, y = strata, label = n.risk)) +
    geom_text(size = 2.2, color = "black") +
    scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
    scale_y_discrete(labels = c("Not Cleared" = "Not Clr",
                                "Cleared"     = "Cleared")) +
    labs(x = NULL, y = NULL, title = "No. at risk") +
    theme_pub() +
    theme(
      plot.title         = element_text(size = 6, face = "bold", hjust = 0),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.line.x        = element_blank(),
      axis.text.y        = element_text(size = 6, color = "black"),
      panel.grid.major   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(t = 0, r = 6, b = 6, l = 6, unit = "pt")
    )
  
  # Inner grid: KM curve (top) + risk table (bottom)
  panel_raw <- plot_grid(
    km_plot, risk_tbl,
    ncol        = 1,
    rel_heights = c(4, 1)
  )
  
  # Attach external title/subtitle strips
  add_panel_title(
    panel_raw,
    paste0(panel_label, "   ", landmark_t, "-Month Landmark: ", outcome_short)
  )
}


# =============================================================================
# 6. Create 6-month landmark panels
# =============================================================================

panel_a <- make_lm_panel(6, "os_months",  "os_status",  "OS",  "A")
panel_b <- make_lm_panel(6, "pfs_months", "pfs_status", "PFS", "B")


# =============================================================================
# 7. Assemble Supplementary Figure (1×2 layout)
# =============================================================================

supp_figure <- plot_grid(
  panel_a, panel_b,
  nrow = 1,
  align = "h",
  axis = "tb"
)

print(supp_figure)


# =============================================================================
# 8. Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Supplementary_Figure2_final.pdf"),
  plot     = supp_figure,
  width    = 8.5, height = 5,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Supplementary_Figure2_final.png"),
  plot     = supp_figure,
  width    = 8.5, height = 5,
  dpi      = 300
)

cat("Supplementary figure saved to:", output_dir, "\n")


# =============================================================================
# 9. Median OS and PFS at 6-month landmark
# =============================================================================

get_lm_medians <- function(landmark_t, time_col, status_col, outcome_short) {
  
  within_lm_ids <- ctdna_raw %>%
    filter(record_id %in% cohort$record_id,
           time_from_evp > 0, time_from_evp <= min(landmark_t, 3)) %>%
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
  get_lm_medians(6, "os_months",  "os_status",  "OS"),
  get_lm_medians(6, "pfs_months", "pfs_status", "PFS")
)

cat("\nMedian survival at 6-month landmark by clearance group:\n")
print(median_rows)

# Create gt table
median_tbl <- median_rows %>%
  gt() %>%
  tab_header(
    title    = md("**Median Survival at 6-Month Landmark by ctDNA Clearance Status**"),
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
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  opt_stylize(style = 1, color = "gray") %>%
  tab_options(table.font.size = px(11))

median_tbl

gtsave(median_tbl,
       file.path(output_dir, "Supplementary_6mo_Landmark_median_table.docx"))
gtsave(median_tbl,
       file.path(output_dir, "Supplementary_6mo_Landmark_median_table.png"),
       vwidth = 700, zoom = 2)

cat("Median survival table saved to:", output_dir, "\n")