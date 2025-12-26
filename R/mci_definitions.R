# R/mci_definitions.R
# Define MCI onset ages under:
#   - primary: your "confirmed at next available evaluation or death" logic (approx)
#   - first MCI: first observed MCI
#   - confirmed within 1y/2y: first MCI that is later confirmed (MCI/dementia) within window

suppressPackageStartupMessages({
  library(dplyr)
})

as_dx_int <- function(x) {
  # dcfdx_3gp: 1=NCI, 2=MCI, 3=Dementia
  suppressWarnings(as.integer(x))
}

first_mci_onset_age <- function(df_id) {
  df_id <- df_id |> arrange(age_at_visit)
  dx <- as_dx_int(df_id$dcfdx_3gp)
  idx <- which(dx == 2)
  if (length(idx) == 0) return(NA_real_)
  df_id$age_at_visit[min(idx)]
}

# Confirmed MCI onset age:
# Find ANY MCI visit i such that there exists a later visit j (j>i) within window_years
# where dx is MCI or Dementia (dx>=2).
# If the next visit is NCI, we continue searching forward.
confirmed_mci_onset_age <- function(df_id, window_years = Inf) {
  df_id <- df_id |> arrange(age_at_visit)
  age <- df_id$age_at_visit
  dx  <- as_dx_int(df_id$dcfdx_3gp)

  # candidate MCI positions
  mci_pos <- which(dx == 2)
  if (length(mci_pos) == 0) return(NA_real_)

  for (i in mci_pos) {
    ai <- age[i]
    # search forward within window
    jset <- which(age > ai & (age - ai) <= window_years)
    if (length(jset) == 0) next
    # look for first confirming dx>=2 (MCI or dementia)
    if (any(dx[jset] >= 2, na.rm = TRUE)) {
      return(ai)  # onset age is the *first* MCI that later confirms
    }
  }
  NA_real_
}

# Primary onset:
# Use your already derived event time (incmci_iii/time2mci_iii) whenever possible,
# otherwise fall back to confirmed_mci_onset_age(window=Inf).
primary_mci_onset_age <- function(df_id) {
  df_id <- df_id |> arrange(age_at_visit)
  # if you already computed incident + time2 variable at baseline row
  if (all(c("incmci_iii","time2mci_iii","age_bl") %in% names(df_id))) {
    # take first row info
    inc <- df_id$incmci_iii[1]
    if (!is.na(inc) && inc == 1 && !is.na(df_id$time2mci_iii[1])) {
      return(df_id$age_bl[1] + df_id$time2mci_iii[1])
    }
  }
  confirmed_mci_onset_age(df_id, window_years = Inf)
}
