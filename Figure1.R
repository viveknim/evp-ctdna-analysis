# ============================================================================
# Figure 1 (Modified) — Publication Quality Combined Figure
# EVP ctDNA Project
#
# Panel A: Swimmer's plot — first-line EVP, ctDNA within 3 months
# Panel B: Multivariable TD Cox forest plot — Overall Survival
# Panel C: Multivariable TD Cox forest plot — Progression-Free Survival
#
# Target size: 8.5 × 7 inches (~2/3 of a letter page)
# Layout: Panel A (left, full height) | Panels B + C (right, stacked)
# ============================================================================

library(dplyr)
library(tidyverse)
library(survival)
library(broom)
library(lubridate)
library(ggplot2)
library(cowplot)
library(scales)
library(gt)

# ---- USER: Update these paths -----------------------------------------------
data_dir      <- "data"
outcomes_path <- file.path(data_dir, "ev_ctdna_survival_tab.csv")
ctdna_path    <- file.path(data_dir, "ctdna_data.csv")
output_dir    <- "output"
# -----------------------------------------------------------------------------

fmt_p <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))


# =============================================================================
# Shared publication theme (base_size = 7 for compact multi-panel figure)
# =============================================================================

theme_pub <- function() {
  theme_classic(base_size = 7) +
    theme(
      plot.title         = element_text(face = "bold", size = 7, hjust = 0),
      plot.subtitle      = element_text(size = 6, color = "gray40", hjust = 0),
      axis.title         = element_text(face = "bold", size = 6.5),
      axis.text          = element_text(size = 6, color = "black"),
      legend.title       = element_text(face = "bold", size = 6),
      legend.text        = element_text(size = 5.5),
      legend.key.size    = unit(0.35, "cm"),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.25),
      plot.margin        = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt")
    )
}


# =============================================================================
# 1. Load outcomes
# =============================================================================

outcomes <- read_csv(outcomes_path) %>%
  mutate(
    evp_c1d1      = as.Date(evp_c1d1,     format = "%m/%d/%y"),
    os_event_date = as.Date(os_event_date, format = "%m/%d/%y"),
    pfs_date      = as.Date(pfs_date,      format = "%m/%d/%y"),
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

cohort_base <- outcomes %>%
  filter(drug == 1, tx_line == 1) %>%
  semi_join(within_3mo_ids, by = "record_id")

cat("Cohort n =", nrow(cohort_base), "\n")


# =============================================================================
# 4. Derive clinical covariates (for forest plots)
# =============================================================================

cohort_cov <- cohort_base %>%
  mutate(
    dob        = as.Date(dob, format = "%m/%d/%y"),
    dob        = if_else(year(dob) > year(Sys.Date()), dob - years(100), dob),
    age_at_evp = floor(as.numeric(difftime(evp_c1d1, dob, units = "days")) / 365.25),
    age_cat    = if_else(age_at_evp < 75, "< 75 years", "\u226575 years"),
    sex = case_when(gender == 1 ~ "Male", gender == 2 ~ "Female"),
    race_var = case_when(race == 1 ~ "White", race == 2 ~ "Black", race == 3 ~ "Other"),
    var_histology = case_when(
      variant_histology == 0 ~ "No Variant Histology",
      variant_histology == 1 ~ "Variant Histology"),
    primary_site = case_when(
      primary_location___1 == 1 ~ "Lower Tract",
      primary_location___1 == 0 ~ "Upper Tract"),
    metastatic_burden = case_when(
      metastasis_present == 2 ~ "Visceral",
      TRUE                    ~ "No Visceral Disease Burden"),
    ecog = case_when(ps == 0 ~ "0-1", ps == 1 ~ "0-1", TRUE ~ "2+"),
    bmi_cat = case_when(
      bmi < 20              ~ "BMI < 20",
      bmi >= 20 & bmi <= 30 ~ "BMI 20-30",
      bmi > 30              ~ "BMI > 30"),
    renal_function = case_when(
      egfr > 60               ~ "eGFR > 60",
      egfr >= 30 & egfr <= 60 ~ "eGFR 30-60",
      egfr < 30               ~ "eGFR < 30"),
  ) %>%
  mutate(
    age_cat           = relevel(as.factor(age_cat),           ref = "< 75 years"),
    sex               = relevel(as.factor(sex),               ref = "Male"),
    race_var          = relevel(as.factor(race_var),          ref = "White"),
    var_histology     = relevel(as.factor(var_histology),     ref = "No Variant Histology"),
    primary_site      = relevel(as.factor(primary_site),      ref = "Lower Tract"),
    metastatic_burden = relevel(as.factor(metastatic_burden), ref = "No Visceral Disease Burden"),
    ecog              = relevel(as.factor(ecog),              ref = "0-1"),
    bmi_cat           = relevel(as.factor(bmi_cat),           ref = "BMI 20-30"),
    renal_function    = relevel(as.factor(renal_function),    ref = "eGFR 30-60")
  )


# =============================================================================
# 5. Identify first ctDNA clearance
# =============================================================================

first_clearance <- ctdna_raw %>%
  filter(record_id %in% cohort_base$record_id,
         time_from_evp > 0, ctdna_value == 0) %>%
  group_by(record_id) %>%
  arrange(time_from_evp) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, time_to_clearance = time_from_evp)


# =============================================================================
# 6. Build time-dependent datasets and fit Cox models
# =============================================================================

covariates_all <- c("age_cat", "sex", "race_var", "var_histology", "primary_site",
                    "metastatic_burden", "ecog", "bmi_cat", "renal_function")

make_td_dataset <- function(outcome_col, status_col) {
  df_base <- cohort_cov %>%
    select(record_id, time = all_of(outcome_col), status = all_of(status_col),
           all_of(covariates_all))

  df_td <- tmerge(df_base, df_base,
                  id = record_id,
                  ev = event(time, status))
  df_td <- tmerge(df_td, first_clearance,
                  id      = record_id,
                  cleared = tdc(time_to_clearance))
  df_td %>%
    mutate(
      clearance_status = factor(
        ifelse(cleared == 1, "Cleared", "Never Cleared"),
        levels = c("Never Cleared", "Cleared")
      )
    ) %>%
    filter(if_all(all_of(covariates_all), ~ !is.na(.)))
}

df_os  <- make_td_dataset("os_months",  "os_status")
df_pfs <- make_td_dataset("pfs_months", "pfs_status")

fit_mv <- function(df) {
  coxph(
    Surv(tstart, tstop, ev) ~ clearance_status +
      age_cat + sex + race_var + var_histology + primary_site +
      metastatic_burden + ecog + bmi_cat + renal_function,
    data = df
  )
}

model_os  <- fit_mv(df_os)
model_pfs <- fit_mv(df_pfs)

cat("OS model  — N:", length(unique(df_os$record_id)),
    "| Events:", model_os$nevent, "\n")
cat("PFS model — N:", length(unique(df_pfs$record_id)),
    "| Events:", model_pfs$nevent, "\n\n")


# =============================================================================
# 7. Helper: build compact forest plot panel
# =============================================================================

term_labels <- c(
  "clearance_statusCleared"               = "ctDNA Cleared vs Never Cleared",
  "age_cat\u226575 years"                 = "Age \u226575 vs <75 years",
  "sexFemale"                             = "Sex: Female vs Male",
  "race_varBlack"                         = "Race: Black vs White",
  "race_varOther"                         = "Race: Other vs White",
  "var_histologyVariant Histology"        = "Variant Histology vs None",
  "primary_siteUpper Tract"               = "Primary Site: Upper vs Lower",
  "metastatic_burdenVisceral"             = "Visceral vs No Visceral Mets",
  "ecog2+"                                = "ECOG PS: 2+ vs 0\u20131",
  "bmi_catBMI < 20"                       = "BMI <20 vs 20\u201330",
  "bmi_catBMI > 30"                       = "BMI >30 vs 20\u201330",
  "renal_functioneGFR > 60"               = "eGFR >60 vs 30\u201360",
  "renal_functioneGFR < 30"               = "eGFR <30 vs 30\u201360"
)

make_forest_panel <- function(model, df, panel_title) {

  n_pts   <- length(unique(df$record_id))
  n_ev    <- model$nevent
  cindex  <- round(summary(model)$concordance[1], 3)

  fp <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      label = recode(term, !!!term_labels),
      label = ifelse(label == term, term, label),
      sig_group = ifelse(p.value < 0.05, "p < 0.05", "p \u2265 0.05"),
      # Clamp CI bounds for log scale display only
      conf.low  = pmax(conf.low,  0.01),
      conf.high = pmin(conf.high, 200)
    )

  fp$label <- factor(fp$label, levels = rev(fp$label))

  x_lo <- min(0.05, min(fp$conf.low,  na.rm = TRUE) * 0.7)
  x_hi <- max(10,   max(fp$conf.high, na.rm = TRUE) * 1.3)

  ggplot(fp, aes(y = label, x = estimate, xmin = conf.low, xmax = conf.high)) +

    geom_vline(xintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.35) +

    geom_errorbarh(aes(color = sig_group), height = 0.22, linewidth = 0.5) +

    geom_point(aes(color = sig_group, shape = sig_group), size = 1.8) +

    scale_color_manual(
      values = c("p < 0.05" = "#C0392B", "p \u2265 0.05" = "#2C3E50"),
      name   = "Significance"
    ) +
    scale_shape_manual(
      values = c("p < 0.05" = 18, "p \u2265 0.05" = 15),
      name   = "Significance"
    ) +
    scale_x_log10(
      name   = "Hazard Ratio (log scale)",
      breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
      labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8"),
      limits = c(x_lo, x_hi * 1.1)
    ) +
    labs(
      y        = NULL,
      title    = panel_title
    ) +
    theme_pub() +
    theme(
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      plot.margin        = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"),
      panel.grid.major.y = element_line(color = "gray94", linewidth = 0.2),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 5.5)
    )
}


# =============================================================================
# PANEL B: OS forest plot
# =============================================================================

panel_b <- make_forest_panel(
  model       = model_os,
  df          = df_os,
  panel_title = "B   Time dependent CPH model of OS"
)


# =============================================================================
# PANEL C: PFS forest plot
# =============================================================================

panel_c <- make_forest_panel(
  model       = model_pfs,
  df          = df_pfs,
  panel_title = "C   Time dependent CPH model of PFS"
)


# =============================================================================
# PANEL A: Swimmer's plot (compact publication version)
# =============================================================================

cohort_sw <- cohort_base %>%
  arrange(desc(os_months)) %>%
  mutate(patient_rank = row_number())

n_pts_sw <- nrow(cohort_sw)

ctdna_pts <- ctdna_raw %>%
  inner_join(cohort_sw %>% select(record_id, patient_rank, os_months),
             by = "record_id") %>%
  filter(time_from_evp <= os_months) %>%
  group_by(record_id) %>%
  filter(
    time_from_evp > 0 |
    (time_from_evp <= 0 &
       time_from_evp == max(time_from_evp[time_from_evp <= 0], na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(ctdna_status = ifelse(ctdna_value == 0,
                               "Negative (undetectable)",
                               "Positive (detectable)"))

bar_starts <- ctdna_pts %>%
  group_by(record_id) %>%
  summarise(min_time = min(time_from_evp, na.rm = TRUE), .groups = "drop") %>%
  mutate(bar_start = ifelse(min_time < 0, min_time, 0))

cohort_sw <- cohort_sw %>%
  left_join(bar_starts %>% select(record_id, bar_start), by = "record_id") %>%
  mutate(bar_start = replace_na(bar_start, 0))

draw_segments <- ctdna_pts %>%
  arrange(record_id, time_from_evp) %>%
  group_by(record_id) %>%
  mutate(x_end = lead(time_from_evp)) %>%
  ungroup() %>%
  mutate(
    x_start   = time_from_evp,
    x_end     = ifelse(is.na(x_end), os_months, x_end),
    seg_color = ifelse(ctdna_value == 0,
                       "Negative (undetectable)",
                       "Positive (detectable)")
  ) %>%
  filter(x_start < x_end) %>%
  select(record_id, patient_rank, x_start, x_end, seg_color)

pre_draw_segs <- cohort_sw %>%
  select(record_id, patient_rank, bar_start) %>%
  left_join(
    ctdna_pts %>%
      group_by(record_id) %>%
      summarise(first_draw = min(time_from_evp), .groups = "drop"),
    by = "record_id"
  ) %>%
  filter(bar_start == 0, first_draw > 0) %>%
  transmute(record_id, patient_rank,
            x_start = 0, x_end = first_draw,
            seg_color = "No Preceding ctDNA")

all_segments <- bind_rows(draw_segments, pre_draw_segs)

pfs_events <- cohort_sw %>%
  filter(pfs_status == 1, os_status == 0) %>%
  select(patient_rank, pfs_months)

os_events <- cohort_sw %>%
  filter(os_status == 1) %>%
  select(patient_rank, os_months)

x_max_sw <- ceiling(max(cohort_sw$os_months, na.rm = TRUE) / 6) * 6
x_min_sw <- floor(min(cohort_sw$bar_start, na.rm = TRUE))

# Panel A: embed "A" into title to avoid cowplot label overlap
panel_a <- ggplot() +

  geom_segment(
    data = all_segments,
    aes(x = x_start, xend = x_end,
        y = patient_rank, yend = patient_rank,
        color = seg_color),
    linewidth = 1.2, lineend = "butt"
  ) +

  geom_point(
    data = ctdna_pts,
    aes(x = time_from_evp, y = patient_rank, fill = ctdna_status),
    shape = 21, size = 1.3, stroke = 0.5, color = "black"
  ) +

  geom_point(
    data = cohort_sw,
    aes(x = 0, y = patient_rank),
    shape = 17, size = 1.6, color = "#2166AC"
  ) +

  geom_point(
    data = pfs_events,
    aes(x = pfs_months, y = patient_rank),
    shape = 25, size = 1.8, fill = "#D32F2F", color = "black"
  ) +

  geom_point(
    data = os_events,
    aes(x = os_months, y = patient_rank),
    shape = 4, size = 1.8, color = "black", stroke = 0.9
  ) +

  scale_fill_manual(
    values = c("Positive (detectable)"   = "black",
               "Negative (undetectable)" = "white"),
    labels = c("Positive (detectable)"   = "Detectable",
               "Negative (undetectable)" = "Undetectable"),
    name  = "ctDNA",
    guide = guide_legend(
      override.aes = list(shape = 21, size = 3, stroke = 0.7, color = "black")
    )
  ) +

  scale_color_manual(
    values = c("Positive (detectable)"   = "#FF9800",
               "Negative (undetectable)" = "#4CAF50",
               "No Preceding ctDNA"      = "gray80"),
    labels = c("Positive (detectable)"   = "Detectable",
               "Negative (undetectable)" = "Undetectable",
               "No Preceding ctDNA"      = "No Preceding ctDNA"),
    name  = "Preceding ctDNA",
    guide = guide_legend(override.aes = list(linewidth = 2))
  ) +

  scale_x_continuous(
    name   = "Time from EVP Initiation (months)",
    breaks = seq(0, x_max_sw, by = 6),
    limits = c(x_min_sw - 0.3, x_max_sw + 0.5),
    expand = expansion(add = c(0, 0))
  ) +

  scale_y_reverse(
    breaks = NULL,
    expand = expansion(add = c(0.8, 0.8))
  ) +

  labs(
    title = paste0("A   Cohort of First Line EVP patients with ctDNA testing (n = ",
                   n_pts_sw, ")"),
    y     = NULL
  ) +

  theme_pub() +
  theme(
    axis.title.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line.y     = element_blank(),
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.justification = "left",
    legend.box.just = "left",
    legend.spacing  = unit(0.1, "cm"),
    plot.margin     = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt")
  )

# Marker legend (EVP start, PFS event, OS death)
marker_df <- data.frame(
  x     = c(1, 1, 1),
  y     = c(3, 2, 1),
  label = c("EVP Start (C1D1)", "PFS Event", "OS Death")
)

marker_legend <- ggplot(marker_df, aes(x = x, y = y)) +
  geom_point(aes(shape = label, fill = label, color = label),
             size = 2.2, stroke = 0.8) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.06, size = 1.9) +
  scale_shape_manual(values = c(
    "EVP Start (C1D1)" = 17, "PFS Event" = 25, "OS Death" = 4)) +
  scale_fill_manual(values = c(
    "EVP Start (C1D1)" = "#2166AC", "PFS Event" = "#D32F2F", "OS Death" = "black")) +
  scale_color_manual(values = c(
    "EVP Start (C1D1)" = "#2166AC", "PFS Event" = "black", "OS Death" = "black")) +
  scale_x_continuous(limits = c(0.9, 2)) +
  scale_y_continuous(limits = c(0.5, 3.5)) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0))

panel_a_combined <- plot_grid(
  panel_a,
  marker_legend,
  ncol        = 1,
  rel_heights = c(1, 0.07)
)


# =============================================================================
# ASSEMBLE FIGURE 1
# Panel labels (A/B/C) embedded in plot titles to prevent overlap
# =============================================================================

right_col <- plot_grid(
  panel_b,
  panel_c,
  ncol        = 1,
  rel_heights = c(1, 1),
  align       = "v",
  axis        = "lr"
)

figure_1 <- plot_grid(
  panel_a_combined,
  right_col,
  ncol       = 2,
  rel_widths = c(0.38, 0.62),
  align      = "h",
  axis       = "tb"
)

print(figure_1)


# =============================================================================
# Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Figure1_publication_modified_final.pdf"),
  plot     = figure_1,
  width    = 8.5, height = 7,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Figure1_publication_modified_final.png"),
  plot     = figure_1,
  width    = 8.5, height = 7,
  dpi      = 300
)

cat("Saved to:", output_dir, "\n")


# =============================================================================
# Supplementary: HR tables corresponding to forest plots (Panels B & C)
# =============================================================================

build_hr_table <- function(model, df, outcome_label) {
  n_pts    <- length(unique(df$record_id))
  n_events <- model$nevent
  cindex   <- round(summary(model)$concordance[1], 3)
  cindex_se <- round(summary(model)$concordance[2], 3)

  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value) %>%
    mutate(
      label   = recode(term, !!!term_labels, .default = term),
      hr_ci   = paste0(sprintf("%.2f", estimate),
                       " (", sprintf("%.2f", conf.low),
                       "\u2013", sprintf("%.2f", conf.high), ")"),
      p_fmt   = fmt_p(p.value)
    ) %>%
    select(label, hr_ci, p_fmt) %>%
    gt() %>%
    tab_header(
      title    = md(paste0("**Multivariable Time-Dependent Cox Model — ", outcome_label, "**"))
    ) %>%
    cols_label(
      label  = "Variable",
      hr_ci  = "HR (95% CI)",
      p_fmt  = "p-value"
    ) %>%
    cols_align(align = "left",  columns = label) %>%
    cols_align(align = "center", columns = c(hr_ci, p_fmt)) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    opt_stylize(style = 1, color = "gray") %>%
    tab_options(table.font.size = px(11))
}

tbl_os  <- build_hr_table(model_os,  df_os,  "Overall Survival")
tbl_pfs <- build_hr_table(model_pfs, df_pfs, "Progression-Free Survival")

tbl_os
tbl_pfs

gtsave(tbl_os,
       file.path(output_dir, "Figure1_modified_HR_table_OS.docx"))
gtsave(tbl_pfs,
       file.path(output_dir, "Figure1_modified_HR_table_PFS.docx"))

gtsave(tbl_os,
       file.path(output_dir, "Figure1_modified_HR_table_OS.png"),
       vwidth = 900, zoom = 2)
gtsave(tbl_pfs,
       file.path(output_dir, "Figure1_modified_HR_table_PFS.png"),
       vwidth = 900, zoom = 2)

cat("HR tables saved to:", output_dir, "\n")
