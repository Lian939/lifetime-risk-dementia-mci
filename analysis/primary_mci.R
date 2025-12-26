# analysis/primary_mci.R
# Runs primary MCI CIFs and saves key outputs

suppressPackageStartupMessages({
  library(dplyr)
  source("R/cif_helpers.R")
})

# EXPECTS: mci_set (overall/sex), mci_set1 (race White/Black subset for MCI)
stopifnot(exists("mci_set"))

index_ages <- c(55, 65, 75, 85)

mci_overall <- lapply(index_ages, function(a) {
  estimate_lifetime_risk(mci_set, index_age = a, max_age = 105)
}) |> bind_rows() |>
  left_join(lapply(index_ages, function(a) count_events_after_index(mci_set, a)) |> bind_rows(),
            by = c("index_age","n_at_risk"))

saveRDS(mci_overall, "output/tables/mci_overall_indexage.rds")

# Curves for Figure 1B
fit_mci <- fit_prodlim_lefttrunc(make_index_subset(mci_set, 55))
ages_grid <- 55:105
curve_mci <- get_curve(fit_mci, ages = ages_grid)

saveRDS(curve_mci, "output/figures/curve_mci_overall.rds")
