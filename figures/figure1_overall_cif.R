# figures/figure1_overall_cif.R
# Figure 1: overall cumulative incidence of dementia and MCI + competing death (optional)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  source("R/cif_helpers.R")
})

stopifnot(exists("dementia_set"), exists("mci_set"))

# Dementia overall curve
fit_dem <- fit_prodlim_lefttrunc(make_index_subset(dementia_set, 55))
ages <- 55:105
dem_curve <- get_curve(fit_dem, ages, cause=1) |> mutate(outcome="Dementia")
dem_death <- get_curve(fit_dem, ages, cause=2) |> mutate(outcome="Death without dementia")

# MCI overall curve
fit_mci <- fit_prodlim_lefttrunc(make_index_subset(mci_set, 55))
mci_curve <- get_curve(fit_mci, ages, cause=1) |> mutate(outcome="MCI")
mci_death <- get_curve(fit_mci, ages, cause=2) |> mutate(outcome="Death without MCI")

p1 <- bind_rows(dem_curve, dem_death) |>
  ggplot(aes(x=age, y=est, linetype=outcome)) +
  geom_line(linewidth=0.9) +
  scale_y_continuous(labels = function(x) round(100*x)) +
  labs(x="Age (years)", y="Cumulative incidence (%)", title="Figure 1A. Dementia and competing death") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.02,0.98), legend.justification=c(0,1))

p2 <- bind_rows(mci_curve, mci_death) |>
  ggplot(aes(x=age, y=est, linetype=outcome)) +
  geom_line(linewidth=0.9) +
  scale_y_continuous(labels = function(x) round(100*x)) +
  labs(x="Age (years)", y="Cumulative incidence (%)", title="Figure 1B. MCI and competing death") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.02,0.98), legend.justification=c(0,1))

saveRDS(list(fig1A = p1, fig1B = p2), "output/figures/figure1_overall_cif.rds")

p1
p2
