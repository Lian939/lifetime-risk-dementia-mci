# analysis/mci_sensitivity.R
# Sensitivity suite for MCI definition:
#  1) First MCI (no confirmation)
#  2) Confirmed within 1 year
#  3) Confirmed within 2 years
#  Primary (your confirmed definition)
#
# EXPECTS: mci_set (baseline rows, already filtered mci_bl_new==0 in your pipeline)
#          r3_long_nomiss2 (long visits with dcfdx_3gp, age_at_visit, age_death, age_last, age_bl)
#
# OUTPUT:
#   - output/tables/mci_sensitivity_table.rds
#   - output/figures/mci_sensitivity_overlay.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  source("R/cif_helpers.R")
  source("R/mci_definitions.R")
})

stopifnot(exists("mci_set"), exists("r3_long_nomiss2"))

# Build onset ages per person from long data
long_keep <- r3_long_nomiss2 |>
  filter(projid %in% mci_set$projid) |>
  arrange(projid, age_at_visit)

onsets <- long_keep |>
  group_by(projid) |>
  group_modify(~{
    tibble::tibble(
      age_mci_primary = primary_mci_onset_age(.x),
      age_mci_first   = first_mci_onset_age(.x),
      age_mci_win1    = confirmed_mci_onset_age(.x, window_years = 1),
      age_mci_win2    = confirmed_mci_onset_age(.x, window_years = 2)
    )
  }) |>
  ungroup()

# Join to mci_set to construct alternative status/exit_age
base0 <- mci_set |>
  select(projid, age_bl, age_last, age_death, msex, sex, status, new_age_bl, persontime) |>
  left_join(onsets, by = "projid")

make_alt_surv <- function(base, onset_col) {
  eps <- 1e-6
  base |>
    mutate(
      age_mci = .data[[onset_col]],
      age_cens = ifelse(!is.na(age_death), age_death, age_last),
      exit_age = pmin(age_mci, age_cens, na.rm = TRUE),
      status_num = case_when(
        !is.na(age_mci)   & abs(exit_age - age_mci)   < eps ~ 1L,
        !is.na(age_death) & abs(exit_age - age_death) < eps ~ 2L,
        TRUE ~ 0L
      ),
      new_age_bl = age_bl - 55,
      persontime = exit_age - 55,
      persontime = ifelse(persontime <= new_age_bl, new_age_bl + eps, persontime)
    ) |>
    filter(!is.na(new_age_bl), !is.na(persontime), new_age_bl >= 0)
}

dat_primary <- make_alt_surv(base0, "age_mci_primary")
dat_first   <- make_alt_surv(base0, "age_mci_first")
dat_win1    <- make_alt_surv(base0, "age_mci_win1")
dat_win2    <- make_alt_surv(base0, "age_mci_win2")

# Fit CIFs at age 105
fit_primary <- fit_prodlim_lefttrunc(make_index_subset(dat_primary, 55))
fit_first   <- fit_prodlim_lefttrunc(make_index_subset(dat_first,   55))
fit_win1    <- fit_prodlim_lefttrunc(make_index_subset(dat_win1,    55))
fit_win2    <- fit_prodlim_lefttrunc(make_index_subset(dat_win2,    55))

summ <- dplyr::bind_rows(
  tibble::tibble(def_label = "Primary (confirmed at next evaluation)",
                 get_cif_at_age(fit_primary, 105), n_events = sum(dat_primary$status_num==1, na.rm=TRUE), n_total=nrow(dat_primary)),
  tibble::tibble(def_label = "First MCI (no confirmation)",
                 get_cif_at_age(fit_first, 105), n_events = sum(dat_first$status_num==1, na.rm=TRUE), n_total=nrow(dat_first)),
  tibble::tibble(def_label = "Confirmed within 1 year",
                 get_cif_at_age(fit_win1, 105), n_events = sum(dat_win1$status_num==1, na.rm=TRUE), n_total=nrow(dat_win1)),
  tibble::tibble(def_label = "Confirmed within 2 years",
                 get_cif_at_age(fit_win2, 105), n_events = sum(dat_win2$status_num==1, na.rm=TRUE), n_total=nrow(dat_win2))
) |>
  mutate(
    lifetime_risk_105 = est * 100,
    lcl_105 = lcl * 100,
    ucl_105 = ucl * 100
  )

# Curves for overlay
ages_grid <- 55:105
curves <- list(
  "Primary (confirmed at next evaluation)" = get_curve(fit_primary, ages_grid),
  "First MCI (no confirmation)"            = get_curve(fit_first,   ages_grid),
  "Confirmed within 1 year"                = get_curve(fit_win1,    ages_grid),
  "Confirmed within 2 years"               = get_curve(fit_win2,    ages_grid)
)

saveRDS(list(summaries = summ, curves = curves),
        "output/tables/mci_sensitivity_suite.rds")
