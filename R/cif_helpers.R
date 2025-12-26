# R/cif_helpers.R
# Helpers for Aalen–Johansen cumulative incidence with left truncation (age as time)
# Assumes:
#   - entry age: new_age_bl = age_bl - origin_age
#   - exit time: persontime = exit_age - origin_age
#   - status: 0=censor, 1=event of interest, 2=death without event

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(prodlim)
})

origin_age_default <- 55
eps_default <- 1e-6

ensure_sex <- function(dat, sex_var = "sex", msex_var = "msex") {
  dat <- as.data.frame(dat)
  if (!sex_var %in% names(dat)) {
    if (msex_var %in% names(dat)) {
      dat[[sex_var]] <- factor(dat[[msex_var]]) |>
        forcats::fct_recode("Male" = "1", "Female" = "0")
    }
  }
  dat
}

ensure_status_num <- function(dat, status_var = "status") {
  dat <- as.data.frame(dat)
  if (!status_var %in% names(dat)) stop("Missing status column: ", status_var)
  # status may be factor/character/numeric
  s <- dat[[status_var]]
  if (is.factor(s)) s <- as.character(s)
  if (is.character(s)) s <- suppressWarnings(as.integer(s))
  dat$status_num <- as.integer(s)
  dat
}

# Make an index-age subset: keep those alive + event-free at index age (using entry age).
# With left truncation, simplest consistent filter is: age_bl >= index_age.
make_index_subset <- function(dat, index_age, origin_age = origin_age_default) {
  dat <- as.data.frame(dat)
  if (!("age_bl" %in% names(dat))) stop("Need age_bl in dat")
  dat <- dat |> filter(age_bl >= index_age)
  # if new_age_bl/persontime not present, compute from exit_age if possible
  if (!("new_age_bl" %in% names(dat))) dat$new_age_bl <- dat$age_bl - origin_age
  if (!("persontime" %in% names(dat))) {
    if ("exit_age" %in% names(dat)) dat$persontime <- dat$exit_age - origin_age
  }
  dat
}

# Fit Aalen–Johansen via prodlim with left truncation
fit_prodlim_lefttrunc <- function(dat, origin_age = origin_age_default) {
  dat <- ensure_status_num(dat)
  if (!all(c("new_age_bl", "persontime", "status_num") %in% names(dat))) {
    stop("dat must contain new_age_bl, persontime, status/status_num")
  }
  dat <- dat |> mutate(
    new_age_bl = pmax(new_age_bl, 0),
    # avoid entry==exit edge cases
    persontime = ifelse(persontime <= new_age_bl, new_age_bl + eps_default, persontime)
  )
  prodlim::prodlim(
    Hist(time = persontime, event = status_num, entry = new_age_bl) ~ 1,
    data = dat
  )
}

# Get cuminc estimate and pointwise CI at a target absolute age
# NOTE: prodlim v2025.x does not expose SE/CI from predict() in many cases.
# We return est always; CI may be NA depending on prodlim build.
get_cif_at_age <- function(fit, age_target, origin_age = origin_age_default, cause = 1) {
  t_target <- age_target - origin_age
  est <- as.numeric(predict(fit, type = "cuminc", cause = cause, times = t_target))
  tibble::tibble(est = est, lcl = NA_real_, ucl = NA_real_)
}

# Extract a full curve over a grid of ages (absolute age)
get_curve <- function(fit, ages, origin_age = origin_age_default, cause = 1) {
  times <- ages - origin_age
  est <- as.numeric(predict(fit, type = "cuminc", cause = cause, times = times))
  tibble::tibble(age = ages, est = est)
}

fmt_pct_ci <- function(est, lcl, ucl, digits = 0) {
  if (is.na(est)) return("—")
  estp <- round(est * 100, digits)
  if (is.na(lcl) || is.na(ucl)) return(paste0(estp))
  lclp <- round(lcl * 100, digits)
  uclp <- round(ucl * 100, digits)
  paste0(estp, " (", lclp, ", ", uclp, ")")
}

# Convenience: estimate lifetime risk at max_age for a dataset + index age
estimate_lifetime_risk <- function(dat, index_age, max_age, origin_age = origin_age_default, cause = 1) {
  dsub <- make_index_subset(dat, index_age = index_age, origin_age = origin_age)
  fit  <- fit_prodlim_lefttrunc(dsub, origin_age = origin_age)
  cif  <- get_cif_at_age(fit, age_target = max_age, origin_age = origin_age, cause = cause)
  tibble::tibble(
    index_age = index_age,
    max_age = max_age,
    est = cif$est,
    lcl = cif$lcl,
    ucl = cif$ucl,
    n_at_risk = nrow(dsub)
  )
}

# Count observed events after index age (not projected)
count_events_after_index <- function(dat, index_age) {
  dsub <- make_index_subset(dat, index_age = index_age)
  tibble::tibble(
    index_age = index_age,
    n_at_risk = nrow(dsub),
    n_events_obs = sum(dsub$status_num == 1, na.rm = TRUE)
  )
}

# Plot overlay of multiple curves (expects list of tibbles with age, est)
plot_cif_overlay <- function(curves_named, title = NULL, legend_pos = c(0.02, 0.98), primary_name = NULL) {
  df <- dplyr::bind_rows(lapply(names(curves_named), function(nm) {
    curves_named[[nm]] |> mutate(def_label = nm)
  }))
  df$def_label <- factor(df$def_label, levels = names(curves_named))

  # primary solid, others dashed
  df <- df |> mutate(linetype = ifelse(def_label == primary_name, "solid", "dashed"))

  ggplot(df, aes(x = age, y = est, linetype = linetype, group = def_label)) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(labels = function(x) round(100 * x)) +
    scale_linetype_identity() +
    labs(x = "Age (years)", y = "Cumulative incidence (%)", title = title) +
    theme_bw() +
    theme(
      legend.position = legend_pos,
      legend.justification = c(0, 1),
      legend.title = element_blank()
    )
}
