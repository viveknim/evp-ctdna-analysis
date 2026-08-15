# =============================================================================
# Supplementary Table 1: Clinical Characteristics by Study Cohort
# EVP ctDNA Project
#
# Cohorts compared:
#   Figure 1 Cohort: First-line EVP patients with ≥1 post-EVP ctDNA draw
#                    (the full included cohort for the swimmer's plot and
#                     time-dependent CPH models in Figure 1)
#
#   Figure 3 Cohort: First-line EVP patients with non-zero baseline ctDNA
#                    AND ≥1 post-EVP draw within 3 months of EVP initiation
#                    (the analytic cohort for ctDNA trajectory analyses in
#                     Figure 3)
#
# Output: Supplementary_Table1.docx  +  Supplementary_Table1.png
# =============================================================================

library(dplyr)
library(tidyverse)
library(lubridate)
library(gtsummary)
library(gt)

# ---- USER: Update these paths -----------------------------------------------
data_dir   <- "data"
output_dir <- "output"
# -----------------------------------------------------------------------------


# =============================================================================
# 1. Load and prepare outcomes data
# =============================================================================

outcomes <- read_csv(file.path(data_dir, "ev_ctdna_survival_tab.csv")) %>%
  mutate(
    evp_c1d1      = as.Date(evp_c1d1,      format = "%m/%d/%y"),
    os_event_date = as.Date(os_event_date,  format = "%m/%d/%y"),
    pfs_date      = as.Date(pfs_date,       format = "%m/%d/%y"),
    dob           = as.Date(dob,            format = "%m/%d/%y"),
    dob           = if_else(year(dob) > year(Sys.Date()), dob - years(100), dob),
    age_at_evp    = floor(as.numeric(difftime(evp_c1d1, dob, units = "days")) / 365.25),
    os_months     = as.numeric(os_event_date - evp_c1d1) / 30.44,
    os_status     = ifelse(os_event == 2, 1, 0),
    pfs_months    = as.numeric(pfs_date - evp_c1d1) / 30.44,
    pfs_status    = ifelse(pfs_event == 1 | pfs_event == 2, 1, 0)
  ) %>%
  mutate(
    age_cat = case_when(
      age_at_evp < 75 ~ "< 75 years",
      TRUE            ~ "\u226575 years"),
    sex = case_when(
      gender == 1 ~ "Male",
      gender == 2 ~ "Female"),
    race_var = case_when(
      race == 1 ~ "White",
      race == 2 ~ "Black",
      race == 3 ~ "Other"),
    var_histology = case_when(
      variant_histology == 0 ~ "No Variant Histology",
      variant_histology == 1 ~ "Variant Histology"),
    primary_site = case_when(
      primary_location___1 == 1 ~ "Lower Tract",
      primary_location___1 == 0 ~ "Upper Tract"),
    metastatic_burden = case_when(
      metastasis_present == 2 ~ "Visceral",
      TRUE                    ~ "Lymph Nodes / No Visceral"),
    ecog = case_when(
      ps == 0 ~ "0-1",
      ps == 1 ~ "0-1",
      TRUE    ~ "2+"),
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
      TRUE                   ~ "No Scans")
  )

# First-line EVP patients
fl_evp <- outcomes %>% filter(drug == 1, tx_line == 1)

cat("First-line EVP patients:", nrow(fl_evp), "\n")


# =============================================================================
# 2. Load ctDNA data and compute timing flags
# =============================================================================

ctdna_raw <- read_csv(file.path(data_dir, "ctdna_data.csv")) %>%
  left_join(fl_evp %>% select(record_id, evp_c1d1), by = "record_id") %>%
  filter(record_id %in% fl_evp$record_id) %>%
  mutate(
    ctdna_date    = as.Date(ctdna_date, format = "%m/%d/%y"),
    evp_c1d1      = as.Date(evp_c1d1,  format = "%m/%d/%y"),
    time_from_evp = as.numeric(ctdna_date - evp_c1d1) / 30.44
  ) %>%
  filter(!is.na(redcap_repeat_instrument))


# =============================================================================
# 3. Define Figure 1 Cohort
#    First-line EVP + ≥1 post-EVP ctDNA draw (any time)
# =============================================================================

within_3mo_ids_fig1 <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)

fig1_cohort <- fl_evp %>%
  semi_join(within_3mo_ids_fig1, by = "record_id")

cat("Figure 1 cohort n =", nrow(fig1_cohort), "\n")


# =============================================================================
# 4. Define Figure 3 Cohort
#    First-line EVP + non-zero baseline ctDNA (last pre-EVP draw) +
#    ≥1 post-EVP draw within 3 months
# =============================================================================

within_3mo_ids <- ctdna_raw %>%
  filter(time_from_evp > 0, time_from_evp <= 3) %>%
  distinct(record_id)  # note: same window as fig1; fig3 additionally requires non-zero baseline

# Baseline ctDNA per patient (last pre-EVP draw)
baseline_per_pt <- ctdna_raw %>%
  filter(time_from_evp <= 0) %>%
  group_by(record_id) %>%
  filter(time_from_evp == max(time_from_evp)) %>%
  slice(1) %>%
  ungroup() %>%
  select(record_id, baseline_ctdna = ctdna_value)

# Figure 3 analytic IDs: non-zero baseline + draw within 3 months
fig3_ids <- baseline_per_pt %>%
  filter(!is.na(baseline_ctdna), baseline_ctdna != 0) %>%
  semi_join(within_3mo_ids, by = "record_id") %>%
  pull(record_id)

fig3_cohort <- fl_evp %>%
  filter(record_id %in% fig3_ids) %>%
  left_join(baseline_per_pt, by = "record_id")

cat("Figure 3 cohort n =", nrow(fig3_cohort), "\n\n")


# =============================================================================
# 5. Apply consistent factor ordering to both cohorts
# =============================================================================

set_factors <- function(df) {
  df %>% mutate(
    age_cat = factor(age_cat,
                     levels = c("< 75 years", "\u226575 years")),
    sex = factor(sex,
                 levels = c("Male", "Female")),
    race_var = factor(race_var,
                      levels = c("White", "Black", "Other")),
    var_histology = factor(var_histology,
                           levels = c("No Variant Histology", "Variant Histology")),
    primary_site = factor(primary_site,
                          levels = c("Lower Tract", "Upper Tract")),
    metastatic_burden = factor(metastatic_burden,
                               levels = c("Lymph Nodes / No Visceral", "Visceral")),
    ecog = factor(ecog,
                  levels = c("0-1", "2+")),
    bmi_cat = factor(bmi_cat,
                     levels = c("BMI < 20", "BMI 20-30", "BMI > 30")),
    renal_function = factor(renal_function,
                            levels = c("eGFR > 60", "eGFR 30-60", "eGFR < 30")),
    skin_toxicity_identified = factor(skin_toxicity_identified,
                                      levels = c("No Skin Toxicity", "Skin Toxicity")),
    best_radiologic_response = factor(best_radiologic_response,
                                      levels = c("CR", "PR", "SD", "PD", "No Scans"))
  )
}

fig1_cohort <- set_factors(fig1_cohort)
fig3_cohort <- set_factors(fig3_cohort)


# =============================================================================
# 6. Build per-cohort gtsummary tables
# =============================================================================

table_vars <- c(
  "age_at_evp", "age_cat", "sex", "race_var",
  "var_histology", "primary_site", "metastatic_burden",
  "ecog", "bmi_cat", "renal_function",
  "skin_toxicity_identified", "best_radiologic_response",
  "os_months", "pfs_months"
)

var_labels <- list(
  age_at_evp               ~ "Age at EVP Initiation (years)",
  age_cat                  ~ "Age Category",
  sex                      ~ "Sex",
  race_var                 ~ "Race",
  var_histology            ~ "Variant Histology",
  primary_site             ~ "Primary Site",
  metastatic_burden        ~ "Metastatic Burden",
  ecog                     ~ "ECOG Performance Status",
  bmi_cat                  ~ "Body Mass Index",
  renal_function           ~ "Renal Function",
  skin_toxicity_identified ~ "Skin Toxicity",
  best_radiologic_response ~ "Best Radiologic Response",
  os_months                ~ "Overall Survival (months)",
  pfs_months               ~ "Progression-Free Survival (months)"
)

make_tbl <- function(df) {
  df %>%
    select(any_of(table_vars)) %>%
    tbl_summary(
      label     = var_labels,
      statistic = list(
        all_continuous()  ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits   = list(
        all_continuous()  ~ 1,
        all_categorical() ~ c(0, 1)
      ),
      missing      = "ifany",
      missing_text = "Missing"
    )
}

tbl_fig1 <- make_tbl(fig1_cohort)
tbl_fig3 <- make_tbl(fig3_cohort)


# =============================================================================
# 7. Merge and format Supplementary Table 1
# =============================================================================

supp_table_1 <- tbl_merge(
  tbls = list(tbl_fig1, tbl_fig3),
  tab_spanner = c(
    paste0("**Time Dependent Survival Analysis**<br>",
           "\u22651 ctDNA draw within 3 months<br>(N=",
           nrow(fig1_cohort), ")"),
    paste0("**Exploring ctDNA Trends**<br>",
           "Non-zero baseline ctDNA +<br>",
           "\u22651 draw within 3 months<br>(N=",
           nrow(fig3_cohort), ")")
  )
) %>%
  modify_caption(
    "**Supplementary Table 1. Baseline Clinical Characteristics by Study Cohort**"
  ) %>%
  bold_labels()

supp_table_1


# =============================================================================
# 8. Save
# =============================================================================

supp_table_1 %>%
  as_gt() %>%
  tab_options(
    table.font.size         = px(12),
    column_labels.font.size = px(11)
  ) %>%
  opt_stylize(style = 1, color = "gray") %>%
  gtsave(file.path(output_dir, "Supplementary_Table1.docx"))

supp_table_1 %>%
  as_gt() %>%
  tab_options(
    table.font.size         = px(12),
    column_labels.font.size = px(11)
  ) %>%
  opt_stylize(style = 1, color = "gray") %>%
  gtsave(file.path(output_dir, "Supplementary_Table1.png"),
         vwidth = 1600, zoom = 2)

cat("Saved to:", output_dir, "\n")
