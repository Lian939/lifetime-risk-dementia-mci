# analysis/primary_dementia.R
# Runs primary dementia CIFs and saves key outputs

suppressPackageStartupMessages({
  library(dplyr)
  source("R/cif_helpers.R")
})

# EXPECTS: dementia_set (overall/sex), dementia_set1 (race White/Black subset)
stopifnot(exists("dementia_set"))

index_ages <- c(55, 65, 75, 85)

# Overall
dem_overall <- lapply(index_ages, function(a) {
  estimate_lifetime_risk(dementia_set, index_age = a, max_age = 105)
}) |> bind_rows() |>
  left_join(lapply(index_ages, function(a) count_events_after_index(dementia_set, a)) |> bind_rows(),
            by = c("index_age","n_at_risk"))

saveRDS(dem_overall, "output/tables/dementia_overall_indexage.rds")

# Curves for Figure 1A
fit_dem <- fit_prodlim_lefttrunc(make_index_subset(dementia_set, 55))
ages_grid <- 55:105
curve_dem <- get_curve(fit_dem, ages = ages_grid)

saveRDS(curve_dem, "output/figures/curve_dementia_overall.rds")
