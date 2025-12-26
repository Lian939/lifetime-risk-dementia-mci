# tables/table2_index_age.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gt)
  source("R/cif_helpers.R")
})

stopifnot(exists("dementia_set"), exists("mci_set"), exists("dementia_set1"), exists("mci_set1"))

index_ages <- c(55,65,75,85)

# Max ages per your note
max_age_dem_default <- 105
max_age_dem_black   <- 102

max_age_mci_default <- 105
max_age_mci_black   <- 100
max_age_mci_male    <- 104

# helper to build block
block_for <- function(dat_overall_sex, dat_race, outcome_name, max_default, max_black, max_male = max_default) {
  dat_overall_sex <- ensure_sex(dat_overall_sex)
  dat_overall_sex <- ensure_status_num(dat_overall_sex)
  dat_race        <- ensure_status_num(dat_race)

  # overall
  overall <- lapply(index_ages, function(a){
    fit <- fit_prodlim_lefttrunc(make_index_subset(dat_overall_sex, a))
    cif <- get_cif_at_age(fit, age_target = max_default, cause=1)
    tibble::tibble(section="Overall", group="Overall", index_age=a,
                   risk_ci=fmt_pct_ci(cif$est, cif$lcl, cif$ucl),
                   n_at_risk=nrow(make_index_subset(dat_overall_sex, a)),
                   n_events_obs=sum(make_index_subset(dat_overall_sex, a)$status_num==1, na.rm=TRUE))
  }) |> bind_rows()

  # sex
  sex_levels <- c("Female","Male")
  sex_block <- lapply(sex_levels, function(sx){
    dat_sx <- dat_overall_sex |> filter(sex == sx)
    lapply(index_ages, function(a){
      max_use <- if (sx == "Male") max_male else max_default
      fit <- fit_prodlim_lefttrunc(make_index_subset(dat_sx, a))
      cif <- get_cif_at_age(fit, age_target = max_use, cause=1)
      tibble::tibble(section="Sex", group=sx, index_age=a,
                     risk_ci=fmt_pct_ci(cif$est, cif$lcl, cif$ucl),
                     n_at_risk=nrow(make_index_subset(dat_sx, a)),
                     n_events_obs=sum(make_index_subset(dat_sx, a)$status_num==1, na.rm=TRUE))
    }) |> bind_rows()
  }) |> bind_rows()

  # race (expects dat_race has race factor with "White"/"Black")
  race_block <- lapply(c("White","Black"), function(rc){
    dat_rc <- dat_race |> filter(race == rc)
    lapply(index_ages, function(a){
      max_use <- if (rc == "Black") max_black else max_default
      fit <- fit_prodlim_lefttrunc(make_index_subset(dat_rc, a))
      cif <- get_cif_at_age(fit, age_target = max_use, cause=1)
      tibble::tibble(section="Race", group=rc, index_age=a,
                     risk_ci=fmt_pct_ci(cif$est, cif$lcl, cif$ucl),
                     n_at_risk=nrow(make_index_subset(dat_rc, a)),
                     n_events_obs=sum(make_index_subset(dat_rc, a)$status_num==1, na.rm=TRUE))
    }) |> bind_rows()
  }) |> bind_rows()

  bind_rows(overall, sex_block, race_block) |>
    mutate(outcome = outcome_name,
           index_age = factor(paste0(index_age," years"), levels=paste0(index_ages," years")))
}

tab2_dem <- block_for(
  dat_overall_sex = dementia_set,
  dat_race        = dementia_set1,
  outcome_name    = "Dementia",
  max_default     = max_age_dem_default,
  max_black       = max_age_dem_black,
  max_male        = max_age_dem_default
)

tab2_mci <- block_for(
  dat_overall_sex = mci_set,
  dat_race        = mci_set1,
  outcome_name    = "MCI",
  max_default     = max_age_mci_default,
  max_black       = max_age_mci_black,
  max_male        = max_age_mci_male
)

tab2_long <- bind_rows(tab2_dem, tab2_mci)

# build display (like your earlier table)
tab2_disp <- tab2_long |>
  mutate(row_label = case_when(
    outcome == "Dementia" & section=="Overall" & group=="Overall" ~ "Overall",
    outcome == "Dementia" & section=="Sex" ~ group,
    outcome == "Dementia" & section=="Race" ~ group,
    outcome == "MCI" & section=="Overall" & group=="Overall" ~ "Overall",
    outcome == "MCI" & section=="Sex" ~ group,
    outcome == "MCI" & section=="Race" ~ group,
    TRUE ~ group
  )) |>
  select(outcome, section, row_label, index_age, risk_ci) |>
  tidyr::pivot_wider(names_from = index_age, values_from = risk_ci)

gt_table2 <- tab2_disp |>
  gt(groupname_col = "outcome", rowname_col = "row_label") |>
  tab_header(title = "Table 2. Lifetime risk of dementia and MCI after select index ages to advanced age") |>
  tab_source_note("Estimates are cumulative incidence (%) with 95% confidence intervals, accounting for competing risk of death and left truncation.") |>
  tab_source_note("For Black participants, estimates are reported through the maximum age supported by follow-up (dementia: age 102; MCI: age 100). For male participants in MCI analyses, estimates are reported through age 104.") |>
  opt_table_outline()

saveRDS(list(long = tab2_long, gt = gt_table2), "output/tables/table2_gt.rds")
gt_table2
