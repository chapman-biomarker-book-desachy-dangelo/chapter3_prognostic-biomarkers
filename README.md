# Chapter 3: Prognostic Biomarkers

## Overview

This repository accompanies Chapter 3: Prognostic Biomarkers, which focuses on statistical methods for identifying and validating prognostic biomarkers using binary and time-to-event endpoints.

The chapter emphasizes:

- Use of proper performance measures (PMs) and performance improvement measures (PIMs)
- A two-step modeling framework using likelihood-based testing and nested models to assess the incremental prognostic value of a biomarker in the presence of known prognostic factors.
- The importance of study design and a pre-specified analysis plan in the identification and validation of prognostic biomarkers.
- Strategies for identifying an “optimal” cut-off for continuous biomarkers, while emphasizing the importance of clinical relevance and data-driven considerations.

Illustrations of the chapter's concepts are based on an oncology clinical trial dataset. Since the original data cannot be released publicly due to stakeholder policies, a synthetic dataset (n = 300 subjects) is provided to recreate the analyses and run the accompanying R scripts. Although the dataset was simulated using a Cox model with independent covariates and parameters chosen to closely reflect the target population, it is important to note that the estimates reported in this chapter and its supplement may not exactly match those obtained when running the R code.

---

## Repository Contents

### 1. Synthetic Dataset

- **File:** `simulated_hnca_300.csv`
- **Description:** Simulated dataset with 300 subjects, including all variables as demonstrated in the chapter. The outcome is overall survival (time-to-event).
- **Purpose:** Enables users to recreate the analyses and outputs shown in the chapter and its supplement.

**Notes:**

- The dataset is synthetic and does not represent real patients.
- All variables are complete (no missing data).
- Results are for demonstration purposes only. Estimates obtained using this data set may not match those presented in the chapter and its supplement.

---

### 2. Data Dictionary

- **File:** `Data Dictionary.xlsx`
- **Description:** Defines all variables included in the dataset, including coding and descriptions.

---

### 3. R Scripts

Three R scripts are provided to produce analyses, tables, and figures similar to those used in the examples. Each script is self-contained and uses the synthetic dataset as input.

- Example of performance measures for a biomarker with a binary outcome: `Table S3.3 Binary outcome.R`
- Example of performance measures for a biomarker with a time-to-event outcome (OS): `Table 3.4 Survival outcome.R`
- Cut-off identification for baseline hemoglobin (the outcome is OS): `Figure S3.6 & Table S3.7 hemoglobin cutoff.R`

All analyses were conducted using R version 4.5.1. The packages and functions described in this chapter and its supplement were available and verified at the time of writing. Future updates to these packages may cause code examples to behave differently.

---

### 4. Supplementary Material

- **File:** `Supplement.docx`
- **Description:** Provides further insights into the statistical methods described in the main chapter. Presents supplementary figures and tables to support the analyses referenced in the main text.


