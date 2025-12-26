# tables/table1_demographics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(gt)
})

# EXPECTS: dementia_set (dementia analytic), mci_set (MCI analytic)
stopifnot(exists("dementia_set"), exists("mci_set"))

summ_table1 <- function(dat, outcome = c("Dementia","MCI")) {
  outcome <- match.arg(outcome)

  dat <- dat |> mutate(
    female = ifelse(as.character(msex) == "0" | as.character(sex) == "Female", 1, 0),
    white  = ifelse(race7 == 1, 1, 0),
    latino = ifelse(diversity_grp == 3, 1, 0),
    apoe4  = ifelse(anye4 == 1, 1, 0),
    followup_years = fu_year_last
  )

  # status is 0/1/2 for both datasets
  alive_free <- sum(dat$status == 0, na.rm = TRUE)
  inc_cases  <- sum(dat$status == 1, na.rm = TRUE)
  deaths_wo  <- sum(dat$status == 2, na.rm = TRUE)

  tibble::tibble(
    Characteristic = c(
      "Age at baseline, years",
      "Age at last visit, years",
      "Follow up years",
      "Female, n (%)",
      "White, n (%)",
      "Latino, n (%)",
      "Education, years",
      "APOE ε4 carriers, n (%)",
      "Study, n (%)",
      "ROS",
      "MAP",
      "MARS",
      "AA",
      "LATC",
      paste0("Alive free of ", outcome, ", n (%)"),
      paste0("Incident cases of ", outcome, ", n (%)"),
      paste0("Deaths without ", outcome, ", n (%)"),
      "Age at death, years",
      "Interval between age at last and at death (median [IQR])"
    ),
    value = c(
      sprintf("%.2f (%.2f)", mean(dat$age_bl, na.rm=TRUE), sd(dat$age_bl, na.rm=TRUE)),
      sprintf("%.2f (%.2f)", mean(dat$age_last, na.rm=TRUE), sd(dat$age_last, na.rm=TRUE)),
      sprintf("%.2f (%.2f)", mean(dat$followup_years, na.rm=TRUE), sd(dat$followup_years, na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$female==1, na.rm=TRUE), 100*mean(dat$female==1, na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$white==1, na.rm=TRUE), 100*mean(dat$white==1, na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$latino==1, na.rm=TRUE), 100*mean(dat$latino==1, na.rm=TRUE)),
      sprintf("%.2f (%.2f)", mean(dat$educ, na.rm=TRUE), sd(dat$educ, na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$apoe4==1, na.rm=TRUE), 100*mean(dat$apoe4==1, na.rm=TRUE)),
      "",
      sprintf("%d (%.1f)", sum(dat$study=="ROS", na.rm=TRUE), 100*mean(dat$study=="ROS", na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$study=="MAP", na.rm=TRUE), 100*mean(dat$study=="MAP", na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$study=="MARS", na.rm=TRUE), 100*mean(dat$study=="MARS", na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$study=="AA", na.rm=TRUE), 100*mean(dat$study=="AA", na.rm=TRUE)),
      sprintf("%d (%.1f)", sum(dat$study=="LATC", na.rm=TRUE), 100*mean(dat$study=="LATC", na.rm=TRUE)),
      sprintf("%d (%.1f)", alive_free, 100*alive_free/nrow(dat)),
      sprintf("%d (%.1f)", inc_cases,  100*inc_cases/nrow(dat)),
      sprintf("%d (%.1f)", deaths_wo,  100*deaths_wo/nrow(dat)),
      sprintf("%.2f (%.2f)", mean(dat$age_death, na.rm=TRUE), sd(dat$age_death, na.rm=TRUE)),
      {
        dif <- dat$diff_death_last
        sprintf("%.2f [%.2f, %.2f]", median(dif, na.rm=TRUE),
                quantile(dif, 0.25, na.rm=TRUE),
                quantile(dif, 0.75, na.rm=TRUE))
      }
    )
  )
}

tab_dem <- summ_table1(dementia_set, outcome="Dementia") |> rename(Dementia = value)
tab_mci <- summ_table1(mci_set, outcome="MCI") |> rename(MCI = value)

tab1 <- tab_dem |> select(Characteristic, Dementia) |>
  left_join(tab_mci |> select(Characteristic, MCI), by="Characteristic")

gt_tab1 <- tab1 |>
  gt(rowname_col = "Characteristic") |>
  tab_header(title = "Table 1. Characteristics of study participants in estimating lifetime risk of dementia and MCI") |>
  cols_label(Dementia = sprintf("Dementia (N = %d)", nrow(dementia_set)),
             MCI     = sprintf("MCI (N = %d)", nrow(mci_set))) |>
  tab_source_note("Values are mean (SD) unless otherwise indicated. Follow-up years is based on fu_year_last. Percentages are column percentages. APOE ε4 carrier percentage is calculated among participants with non-missing APOE ε4 status.") |>
  opt_table_outline()

saveRDS(gt_tab1, "output/tables/table1_gt.rds")
gt_tab1
