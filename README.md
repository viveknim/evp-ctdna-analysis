# Real-World ctDNA Dynamics and Survival in Advanced Urothelial Cancer Treated with EV+P

Code repository for:

> **Real-world tumor-informed ctDNA dynamics and survival associations in advanced urothelial cancer treated with first-line enfortumab vedotin and pembrolizumab**
>
> *Published in JCO Precision Oncology*

---

## Overview

This repository contains all R scripts used to generate the figures and supplementary materials in the manuscript. Patient-level data are not included in this repository due to patient privacy requirements (HIPAA). Instructions for requesting access to the underlying data can be provided upon reasonable request to the corresponding author.

---

## Repository Structure

```
.
├── Figure1.R                    # Swimmer plot of ctDNA trajectories; time-dependent Cox
│                                #   forest plots for OS and PFS
├── Figure2.R                    # 1-, 2-, and 3-month landmark Kaplan-Meier analyses
│                                #   of ctDNA clearance vs. OS and PFS
├── Figure3.R                    # Post-EVP ctDNA dynamics: individual trajectories,
│                                #   nadir distribution, time to clearance, post-clearance status
├── Figure4.R                    # Logistic regression forest plot of clinical predictors
│                                #   of ctDNA clearance; baseline ctDNA by clearance group;
│                                #   skin toxicity and radiologic response by clearance group
├── Supplementary_Figure2.R      # 6-month landmark Kaplan-Meier analyses
│                                #   of ctDNA clearance vs. OS and PFS
├── Supplementary_Figure3.R      # ctDNA collection characteristics: baseline ctDNA
│                                #   distribution, draws per patient, inter-draw intervals,
│                                #   time from EVP initiation to first post-EVP draw
└── Supplementary_Table1.R       # Patient characteristics across analytical cohorts
```

---

## Data Requirements

Each script requires the following two input files, placed in a `data/` directory (or update `data_dir` at the top of each script):

| File | Description | Used by |
|------|-------------|---------|
| `ev_ctdna_survival_tab.csv` | Patient-level outcomes and clinical covariates | All scripts |
| `ctdna_data.csv` | Longitudinal ctDNA measurements per patient | All scripts |

Set the paths at the top of each script:

```r
data_dir   <- "data"    # folder containing input CSV files
output_dir <- "output"  # folder where figures and tables will be saved
```

---

## Usage

1. Open any script in R or RStudio.
2. Update `data_dir` and `output_dir` at the top of the script.
3. Run the script. Figures are saved as both `.pdf` and `.png`; tables are saved as both `.docx` and `.png`.

---

## Required R Packages

Install all dependencies with:

```r
install.packages(c(
  "dplyr", "tidyverse", "ggplot2", "cowplot", "scales",
  "survival", "survminer", "gtsummary", "gt",
  "broom", "purrr", "lubridate", "tibble"
))

# Required for exporting gt tables to .png
remotes::install_github("rstudio/webshot2")
```

---

## Session Info

Scripts were developed and tested in R version 4.4. A full session info log can be reproduced by running:

```r
sessionInfo()
```

---

## Contact

For questions regarding the data or analysis, please contact the corresponding author.
