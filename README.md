# Lifetime Risk of Dementia and Mild Cognitive Impairment

This repository contains analysis code supporting the manuscript:

> **Lifetime risk of incident dementia and incident mild cognitive impairment in older adults**

The code reproduces all primary analyses, sensitivity analyses, tables, and figures
reported in the manuscript using harmonized longitudinal data from the
Rush Alzheimer’s Disease Center (RADC) cohorts.

---

## Overview

- Outcomes: Incident dementia and incident mild cognitive impairment (MCI)
- Time scale: Age (years)
- Method: Nonparametric cumulative incidence functions (Aalen–Johansen)
- Features:
  - Left truncation (delayed entry)
  - Competing risk of death
  - Multiple index ages (55, 65, 75, 85)
  - Sensitivity analyses for alternative MCI definitions
  - Stratification by sex, race, Latino ethnicity, and stroke history

---

## Data Availability

Data used for this study can be requested for research purposes through the Rush Alzheimer’s Disease Center Research Resource Sharing Hub (https://www.radc.rush.edu/), subject to institutional approval and data-use agreements.

---

## Repository Structure

- `R/`  
  Core functions for cumulative incidence estimation, MCI definitions,
  sensitivity analyses, and plotting.

- `analysis/`  
  Scripts that generate the main and supplementary analyses.

- `tables/` and `figures/`  
  Code to generate manuscript tables and figures.

- `output/`  
  Destination folder for rendered tables and figures.

---

## Software Requirements

- R (≥ 4.4.0)
- Packages:
  - `prodlim`
  - `survival`
  - `tidyverse`
  - `gt`
  - `splines`
  - `mgcv`

Exact package versions used in the manuscript are recorded in
`session_info.txt`.

---

## Reproducibility Notes

All cumulative incidence estimates account for:
- Competing risk of death
- Left truncation using age at cohort entry
- Administrative censoring at December 2025

Sensitivity analyses reproduce results under alternative MCI definitions
and censoring assumptions, as described in the manuscript.

---

## Citation

If you use this code, please cite the associated manuscript.
