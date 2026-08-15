# ============================================================================
# Figure 4 — Clinical Covariates Associated with ctDNA Clearance
# EVP ctDNA Project
#
# Fixes applied vs. earlier draft:
#   - Panel A title: wrapped to 2 lines (was cut off)
#   - Panel B subtitle: shortened (was cut off)
#   - Panel C title: shortened (was cut off)
#   - Panel D title: shortened (was cut off)
#   - Panel D legend: 2-row layout (was cut off on 1 row)
#   - Panel labels embedded in titles (avoids cowplot label overlap)
#   - plot.title size reduced 13 → 11 for comfortable fit in each half-width panel
#
# Panels:
#   A) Forest plot — univariable logistic ORs for ctDNA clearance
#   B) Boxplot of baseline ctDNA by clearance status
#   C) Skin toxicity by clearance status (stacked bar)
#   D) Best radiologic response by clearance status (stacked bar)
# ============================================================================

library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
library(broom)
library(purrr)
library(lubridate)
library(gt)

# ---- USER: Update these paths -----------------------------------------------
data_dir      <- "data"
outcomes_path <- file.path(data_dir, "ev_ctdna_survival_tab.csv")
ctdna_path    <- file.path(data_dir, "ctdna_data.csv")
output_dir    <- "output"
# -----------------------------------------------------------------------------

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))


# =============================================================================
# Shared theme
# plot.title reduced to size 11 (from 13) so titles fit in 6"-wide panels
# =============================================================================

theme_manuscript <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, color = "gray40"),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(size = 10, color = "black"),
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

COL_CLEARANCE    <- "#2166AC"
COL_NO_CLEARANCE <- "#B2182B"


# =============================================================================
# 1. Load outcomes and define first-line EVP cohort
# =============================================================================

outcomes <- read_csv(outcomes_path) %>%
  mutate(evp_c1d1 = as.Date(evp_c1d1, format = "%m/%d/%y"))

fl_evp <- outcomes %>%
  filter(drug == 1, tx_line == 1)

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
# 3. Define analytic cohort
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

baseline_ctdna_per_pt <- ctdna_raw %>%
  group_by(record_id) %>%
  summarise(
    baseline_ctdna = {
      idx <- which(time_from_evp <= 0)
      if (length(idx) == 0) NA_real_
      else ctdna_value[idx][which.max(time_from_evp[idx])]
    },
    .groups = "drop"
  )

analytic_ids <- baseline_ctdna_per_pt %>%
  filter(!is.na(baseline_ctdna), baseline_ctdna != 0) %>%
  semi_join(within_3mo_ids, by = "record_id") %>%
  pull(record_id)

cat("Analytic cohort n =", length(analytic_ids), "\n\n")


# =============================================================================
# 4. ctDNA clearance status
# =============================================================================

clearance_status <- ctdna_raw %>%
  filter(record_id %in% analytic_ids, time_from_evp > 0) %>%
  group_by(record_id) %>%
  summarise(reached_zero = any(ctdna_value == 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(clearance = ifelse(reached_zero, "ctDNA Clearance", "No ctDNA Clearance"))

cat("Clearance breakdown:\n")
clearance_status %>% count(clearance) %>% print()
cat("\n")


# =============================================================================
# 5. Build patient-level analytic dataset
# =============================================================================

clin_characteristics <- fl_evp %>%
  filter(record_id %in% analytic_ids) %>%
  left_join(baseline_ctdna_per_pt, by = "record_id") %>%
  left_join(clearance_status %>% select(record_id, clearance), by = "record_id") %>%
  mutate(
    dob        = as.Date(dob, format = "%m/%d/%y"),
    dob        = if_else(year(dob) > year(Sys.Date()), dob - years(100), dob),
    age_at_evp = floor(as.numeric(difftime(evp_c1d1, dob, units = "days")) / 365.25),
    age_cat    = if_else(age_at_evp < 75, "< 75 years", "\u226575 years"),
    sex = case_when(gender == 1 ~ "Male", gender == 2 ~ "Female"),
    race_var = case_when(
      race == 1 ~ "White", race == 2 ~ "Black", race == 3 ~ "Other"),
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
    skin_toxicity_identified = case_when(
      skin_toxicity == 1 ~ "Skin Toxicity",
      TRUE               ~ "No Skin Toxicity"),
    best_radiologic_response = case_when(
      best_recist_respo == 1 ~ "CR",
      best_recist_respo == 2 ~ "PR",
      best_recist_respo == 3 ~ "SD",
      best_recist_respo == 4 ~ "PD",
      best_recist_respo == 5 ~ "No Scans"),
    clearance_numeric = ifelse(clearance == "ctDNA Clearance", 1, 0)
  ) %>%
  mutate(
    age_cat        = relevel(as.factor(age_cat),        ref = "< 75 years"),
    sex            = relevel(as.factor(sex),            ref = "Male"),
    race_var       = relevel(as.factor(race_var),       ref = "White"),
    ecog           = relevel(as.factor(ecog),           ref = "0-1"),
    bmi_cat        = relevel(as.factor(bmi_cat),        ref = "BMI 20-30"),
    renal_function = relevel(as.factor(renal_function), ref = "eGFR 30-60")
  )

n_pts     <- nrow(clin_characteristics)
n_cleared <- sum(clin_characteristics$clearance == "ctDNA Clearance",  na.rm = TRUE)
n_noclr   <- sum(clin_characteristics$clearance == "No ctDNA Clearance", na.rm = TRUE)

cat("N =", n_pts, "| Cleared =", n_cleared, "| Not Cleared =", n_noclr, "\n\n")


# =============================================================================
# PANEL A — Forest plot
# FIX: title wrapped to 2 lines; label embedded in title; font size reduced
# =============================================================================

models_list <- list(
  age_cat           = glm(clearance_numeric ~ age_cat,           data = clin_characteristics, family = binomial),
  sex               = glm(clearance_numeric ~ sex,               data = clin_characteristics, family = binomial),
  race_var          = glm(clearance_numeric ~ race_var,          data = clin_characteristics, family = binomial),
  var_histology     = glm(clearance_numeric ~ var_histology,     data = clin_characteristics, family = binomial),
  primary_site      = glm(clearance_numeric ~ primary_site,      data = clin_characteristics, family = binomial),
  metastatic_burden = glm(clearance_numeric ~ metastatic_burden, data = clin_characteristics, family = binomial),
  ecog              = glm(clearance_numeric ~ ecog,              data = clin_characteristics, family = binomial),
  bmi_cat           = glm(clearance_numeric ~ bmi_cat,           data = clin_characteristics, family = binomial),
  renal_function    = glm(clearance_numeric ~ renal_function,    data = clin_characteristics, family = binomial)
)

forest_data <- map_dfr(names(models_list), function(var_name) {
  mod <- models_list[[var_name]]

  # Wald estimates and p-values from tidy()
  est <- broom::tidy(mod, exponentiate = TRUE) %>%
    filter(term != "(Intercept)")

  # Wald CIs via confint.default() — guaranteed Wald (estimate ± 1.96*SE),
  # not profile likelihood; exponentiate to OR scale
  ci <- confint.default(mod) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term != "(Intercept)") %>%
    mutate(conf.low  = exp(`2.5 %`),
           conf.high = exp(`97.5 %`)) %>%
    select(term, conf.low, conf.high)

  est %>%
    left_join(ci, by = "term") %>%
    mutate(variable = var_name)
}) %>%
  mutate(
    clean_label = case_when(
      term == "age_cat\u226575 years"                   ~ "\u226575 vs <75 years",
      term == "sexFemale"                               ~ "Female vs Male",
      term == "race_varBlack"                           ~ "Black vs White",
      term == "race_varOther"                           ~ "Other Race vs White",
      term == "var_histologyVariant Histology"          ~ "Variant Histology",
      term == "primary_siteUpper Tract"                 ~ "Upper vs Lower Tract",
      term == "metastatic_burdenVisceral"               ~ "Visceral vs No Visceral",
      term == "ecog2+"                                  ~ "ECOG 2+ vs 0-1",
      term == "bmi_catBMI < 20"                         ~ "BMI <20 vs 20-30",
      term == "bmi_catBMI > 30"                         ~ "BMI >30 vs 20-30",
      term == "renal_functioneGFR < 30"                 ~ "eGFR <30 vs 30-60",
      term == "renal_functioneGFR > 60"                 ~ "eGFR >60 vs 30-60",
      TRUE ~ term
    ),
    category = case_when(
      variable == "age_cat"           ~ "Age",
      variable == "sex"               ~ "Sex",
      variable == "race_var"          ~ "Race",
      variable == "var_histology"     ~ "Histology",
      variable == "primary_site"      ~ "Primary Site",
      variable == "metastatic_burden" ~ "Metastatic Burden",
      variable == "ecog"              ~ "ECOG PS",
      variable == "bmi_cat"           ~ "BMI",
      variable == "renal_function"    ~ "Renal Function",
      TRUE ~ variable
    )
  ) %>%
  filter(clean_label != "Other Race vs White", !is.na(estimate)) %>%
  arrange(category, clean_label) %>%
  mutate(order = row_number())

panel_a <- ggplot(forest_data, aes(y = reorder(clean_label, -order), x = estimate)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.3, linewidth = 0.7, color = "gray30") +
  geom_point(size = 3.5, shape = 21, fill = COL_CLEARANCE, color = "black") +
  scale_x_log10(
    breaks = c(0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 100),
    labels = c("0.01", "0.1", "0.5", "1", "2", "5", "10", "50", "100")
  ) +
  coord_cartesian(xlim = c(0.01, 100)) +
  labs(
    # FIX: \n splits the title across 2 lines so it fits in a 6"-wide panel
    title    = "A   Clinical Predictors of ctDNA Clearance\n    with EVP",
    # FIX: directional text folded into the x-axis label (2nd line) so it
    # sits immediately below the axis instead of in the distant caption area
    x        = "Odds Ratio (95% CI)\n\u2190 Less likely clearance          More likely clearance \u2192",
    y        = NULL
  ) +
  theme_manuscript() +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    axis.title.x       = element_text(face = "bold", size = 9,
                                      color = "gray30", hjust = 0.5)
  )


# =============================================================================
# PANEL B — Baseline ctDNA boxplot
# FIX: subtitle shortened so it fits without overflow
# =============================================================================

box_data  <- clin_characteristics %>%
  filter(!is.na(baseline_ctdna), !is.na(clearance))
wilcox_p  <- wilcox.test(baseline_ctdna ~ clearance, data = box_data)$p.value

box_stats <- box_data %>%
  group_by(clearance) %>%
  summarise(
    med = median(baseline_ctdna, na.rm = TRUE),
    q1  = quantile(baseline_ctdna, 0.25, na.rm = TRUE),
    q3  = quantile(baseline_ctdna, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fmt_ctdna <- function(x) scales::comma(round(x, 0))

box_subtitle <- paste0(
  "Wilcoxon ", fmt_p(wilcox_p), "\n",
  "Cleared: ", fmt_ctdna(box_stats$med[box_stats$clearance == "ctDNA Clearance"]),
  " [IQR: ",  fmt_ctdna(box_stats$q1[box_stats$clearance == "ctDNA Clearance"]),
  "\u2013",   fmt_ctdna(box_stats$q3[box_stats$clearance == "ctDNA Clearance"]),
  "] MTM/mL | ",
  "Not Cleared: ", fmt_ctdna(box_stats$med[box_stats$clearance == "No ctDNA Clearance"]),
  " [IQR: ",      fmt_ctdna(box_stats$q1[box_stats$clearance == "No ctDNA Clearance"]),
  "\u2013",        fmt_ctdna(box_stats$q3[box_stats$clearance == "No ctDNA Clearance"]),
  "] MTM/mL"
)

panel_b <- ggplot(box_data, aes(x = clearance, y = baseline_ctdna, color = clearance)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, width = 0.5) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.6) +
  annotate("text",
           x     = 1.5,
           y     = max(box_data$baseline_ctdna, na.rm = TRUE) * 1.5,
           label = paste0("Wilcoxon ", fmt_p(wilcox_p)),
           size  = 4, fontface = "italic") +
  scale_color_manual(values = c("ctDNA Clearance"    = COL_CLEARANCE,
                                "No ctDNA Clearance" = COL_NO_CLEARANCE)) +
  # FIX: compact labels ("1K", "1M") reduce tick-label width, which is the
  # dominant source of visual gap between the y-axis title and the plot panel
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(
    title    = "B   Baseline ctDNA by Clearance Status",
    x        = NULL,
    y        = "Baseline ctDNA (MTM/mL, log scale)"
  ) +
  theme_manuscript() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    # FIX: zero margin removes the remaining gap between title and tick labels
    axis.title.y       = element_text(face = "bold", size = 11,
                                      margin = margin(r = 0, unit = "pt"))
  )


# =============================================================================
# PANEL C — Skin toxicity stacked bar
# FIX: title shortened to fit
# =============================================================================

skin_data <- clin_characteristics %>%
  filter(!is.na(skin_toxicity_identified), !is.na(clearance)) %>%
  group_by(clearance, skin_toxicity_identified) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clearance) %>%
  mutate(
    total   = sum(n),
    prop    = n / total,
    percent = prop * 100,
    label   = paste0(n, " (", round(percent, 1), "%)")
  )

skin_p <- fisher.test(
  table(clin_characteristics$clearance,
        clin_characteristics$skin_toxicity_identified)
)$p.value

panel_c <- ggplot(skin_data,
                  aes(x = clearance, y = prop, fill = skin_toxicity_identified)) +
  geom_bar(stat = "identity", position = "fill",
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_fill(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("No Skin Toxicity" = "#E0E0E0", "Skin Toxicity" = "#D32F2F"),
    name   = NULL
  ) +
  labs(
    # FIX: shortened title (original was cut off in 6"-wide panel)
    title    = "C   ctDNA Clearance and Skin Toxicity",
    x        = NULL,
    y        = "Proportion"
  ) +
  theme_manuscript() +
  theme(
    axis.text.x     = element_text(face = "bold"),
    legend.position = "bottom",
    plot.subtitle   = element_text(size = 10, color = "gray30"),
    # FIX: r = 0 eliminates the gap between title and tick labels entirely
    axis.title.y    = element_text(face = "bold", size = 11,
                                   margin = margin(r = 0, unit = "pt"))
  )


# =============================================================================
# PANEL D — Radiologic response stacked bar
# FIX: title shortened; legend placed on 2 rows so it doesn't overflow
# =============================================================================

response_data <- clin_characteristics %>%
  filter(!is.na(best_radiologic_response), !is.na(clearance)) %>%
  group_by(clearance, best_radiologic_response) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(clearance) %>%
  mutate(
    total   = sum(n),
    prop    = n / total,
    percent = prop * 100,
    label   = paste0(n, " (", round(percent, 1), "%)")
  ) %>%
  mutate(best_radiologic_response = factor(
    best_radiologic_response,
    levels = c("CR", "PR", "SD", "PD", "No Scans")
  ))

response_p <- fisher.test(
  table(clin_characteristics$clearance,
        clin_characteristics$best_radiologic_response),
  simulate.p.value = TRUE
)$p.value

panel_d <- ggplot(response_data,
                  aes(x = clearance, y = prop, fill = best_radiologic_response)) +
  geom_bar(stat = "identity", position = "fill",
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_fill(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "CR"       = "#1B5E20",
      "PR"       = "#66BB6A",
      "SD"       = "#FDD835",
      "PD"       = "#FF8A65",
      "No Scans" = "#BDBDBD"
    ),
    name = "Best RECIST Response"
  ) +
  labs(
    # FIX: shortened title (original was cut off in 6"-wide panel)
    title    = "D   ctDNA Clearance and\n    Radiologic Response",
    x        = NULL,
    y        = "Proportion"
  ) +
  theme_manuscript() +
  theme(
    axis.text.x     = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    plot.subtitle   = element_text(size = 10, color = "gray30"),
    # FIX: r = 0 eliminates the gap between title and tick labels entirely
    axis.title.y    = element_text(face = "bold", size = 11,
                                   margin = margin(r = 0, unit = "pt"))
  ) +
  # FIX: 2-row legend so all 5 items fit without horizontal overflow
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


# =============================================================================
# Assemble Figure 4
# Labels embedded in titles — labels argument removed from plot_grid to
# prevent the cowplot label from overlapping with panel title text
# =============================================================================

figure_4 <- plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol  = 2,
  align = "hv"
)

print(figure_4)


# =============================================================================
# Save
# =============================================================================

ggsave(
  filename = file.path(output_dir, "Figure4_publication_final.pdf"),
  plot     = figure_4,
  width    = 12, height = 10,
  device   = "pdf"
)

ggsave(
  filename = file.path(output_dir, "Figure4_publication_final.png"),
  plot     = figure_4,
  width    = 12, height = 10,
  dpi      = 300
)

cat("Saved to:", output_dir, "\n")


# =============================================================================
# OR table corresponding to Panel A forest plot
# =============================================================================

or_table <- forest_data %>%
  select(clean_label, category, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    or_ci = paste0(sprintf("%.2f", estimate),
                   " (", sprintf("%.2f", conf.low),
                   "\u2013", sprintf("%.2f", conf.high), ")"),
    p_fmt = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  select(category, clean_label, or_ci, p_fmt) %>%
  gt(groupname_col = "category") %>%
  tab_header(
    title    = md("**Figure 4A — Clinical Predictors of ctDNA Clearance**")
  ) %>%
  cols_label(
    clean_label = "Comparison",
    or_ci       = "OR (95% CI)",
    p_fmt       = "p-value"
  ) %>%
  cols_align(align = "left",   columns = clean_label) %>%
  cols_align(align = "center", columns = c(or_ci, p_fmt)) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold", color = "gray30"),
    locations = cells_row_groups()
  ) %>%
  tab_footnote(
    footnote = "OR > 1 indicates greater odds of ctDNA clearance. CIs and p-values derived from Wald method."
  ) %>%
  opt_stylize(style = 1, color = "gray") %>%
  tab_options(table.font.size = px(11))

or_table

gtsave(or_table,
       file.path(output_dir, "Figure4A_OR_table.docx"))
gtsave(or_table,
       file.path(output_dir, "Figure4A_OR_table.png"),
       vwidth = 900, zoom = 2)

cat("OR table saved to:", output_dir, "\n")
