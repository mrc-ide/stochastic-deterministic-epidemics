orderly2::orderly_parameters(assumptions = "central", deterministic = TRUE)

orderly2::orderly_shared_resource(util.R = "util.R")

orderly2::orderly_resource(c("pars/alpha_ve_high/deterministic/info.csv", "pars/alpha_ve_high/deterministic/proposal.csv", "pars/alpha_ve_high/stochastic/info.csv", "pars/alpha_ve_high/stochastic/proposal.csv", "pars/alpha_ve_low/deterministic/info.csv", "pars/alpha_ve_low/deterministic/proposal.csv", "pars/alpha_ve_low/stochastic/info.csv", "pars/alpha_ve_low/stochastic/proposal.csv", "pars/booster_ve_high/deterministic/info.csv", "pars/booster_ve_high/deterministic/proposal.csv", "pars/booster_ve_high/stochastic/info.csv", "pars/booster_ve_high/stochastic/proposal.csv", "pars/booster_ve_low/deterministic/info.csv", "pars/booster_ve_low/deterministic/proposal.csv", "pars/booster_ve_low/stochastic/info.csv", "pars/booster_ve_low/stochastic/proposal.csv", "pars/central/deterministic/info.csv", "pars/central/deterministic/proposal.csv", "pars/central/stochastic/info.csv", "pars/central/stochastic/proposal.csv", "pars/crim_death_high/deterministic/info.csv", "pars/crim_death_high/deterministic/proposal.csv", "pars/crim_death_high/stochastic/info.csv", "pars/crim_death_high/stochastic/proposal.csv", "pars/crim_death_low/deterministic/info.csv", "pars/crim_death_low/deterministic/proposal.csv", "pars/crim_death_low/stochastic/info.csv", "pars/crim_death_low/stochastic/proposal.csv", "pars/crim_hospi_high/deterministic/info.csv", "pars/crim_hospi_high/deterministic/proposal.csv", "pars/crim_hospi_high/stochastic/info.csv", "pars/crim_hospi_high/stochastic/proposal.csv", "pars/crim_hospi_low/deterministic/info.csv", "pars/crim_hospi_low/deterministic/proposal.csv", "pars/crim_hospi_low/stochastic/info.csv", "pars/crim_hospi_low/stochastic/proposal.csv", "pars/crim_infect_high/deterministic/info.csv", "pars/crim_infect_high/deterministic/proposal.csv", "pars/crim_infect_high/stochastic/info.csv", "pars/crim_infect_high/stochastic/proposal.csv", "pars/crim_infect_low/deterministic/info.csv", "pars/crim_infect_low/deterministic/proposal.csv", "pars/crim_infect_low/stochastic/info.csv", "pars/crim_infect_low/stochastic/proposal.csv", "pars/delta_ve_high/deterministic/info.csv", "pars/delta_ve_high/deterministic/proposal.csv", "pars/delta_ve_high/stochastic/info.csv", "pars/delta_ve_high/stochastic/proposal.csv", "pars/delta_ve_low/deterministic/info.csv", "pars/delta_ve_low/deterministic/proposal.csv", "pars/delta_ve_low/stochastic/info.csv", "pars/delta_ve_low/stochastic/proposal.csv", "pars/fixed_si_high/deterministic/info.csv", "pars/fixed_si_high/deterministic/proposal.csv", "pars/fixed_si_high/stochastic/info.csv", "pars/fixed_si_high/stochastic/proposal.csv", "pars/fixed_si_low/deterministic/info.csv", "pars/fixed_si_low/deterministic/proposal.csv", "pars/fixed_si_low/stochastic/info.csv", "pars/fixed_si_low/stochastic/proposal.csv", "pars/mu_d_summer/deterministic/info.csv", "pars/mu_d_summer/deterministic/proposal.csv", "pars/mu_d_summer/stochastic/info.csv", "pars/mu_d_summer/stochastic/proposal.csv", "pars/mu_d_winter/deterministic/info.csv", "pars/mu_d_winter/deterministic/proposal.csv", "pars/mu_d_winter/stochastic/info.csv", "pars/mu_d_winter/stochastic/proposal.csv", "data/vaccine_efficacy_alpha.csv", "data/vaccine_efficacy_delta.csv", "data/vaccine_efficacy_omicron.csv", "data/vaccine_uptake.csv", "data/support_severity.csv"))

orderly2::orderly_dependency(
  "sarscov2_data",
  "latest",
  c(data_vaccination.csv = "outputs/data_vaccination.csv",
    weighted_prior_ranges.csv = "outputs/weighted_prior_ranges.csv"))

orderly2::orderly_artefact(
  "fitted hyperparameters for priors",
  c("parameters_base.rds", "parameters_info.csv", "parameters_prior.csv", "parameters_proposal.csv", "parameters_transform.R"))
orderly2::orderly_artefact("supplementary figures", "fig_sup_vacc_age.png")

library(sircovid)
library(spimalot)
library(tidyr)
library(dplyr)
library(forcats)
library(magrittr)
library(ggplot2)
library(scales)

orderly2::orderly_resource("R/support.R")
orderly2::orderly_resource("R/priors.R")
orderly2::orderly_resource("R/baseline.R")
orderly2::orderly_resource("R/transform.R")
orderly2::orderly_resource("R/vaccine.R")
source("R/support.R")
source("R/priors.R")
source("R/baseline.R")
source("R/transform.R")
source("R/vaccine.R")

source("util.R")

version_check("sircovid", "0.15.1")
version_check("spimalot", "0.8.30")

## Define date at which the data is capped for analysis
date <- "2022-02-24"

## Five epochs after starting with a single strain model (without vaccination)
## * mid August 2020: Alpha appears, expand strains
## * early December 2020: vaccination starts, expand vaccine classes 
##                        but without boosters
## * early March 2021: delta appears, expand strains
## * mid September 2021: booster programme starts, expand vaccine classes
## * early November 2021: omicron appears, rotate strains
epoch_dates <- c("2020-09-17", "2020-12-07", "2021-03-08", "2021-09-14", "2021-11-01")

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(assumptions, deterministic)

## The baselines are always region-specific
regions <- sircovid::regions("england")

baseline <- lapply(regions, create_baseline,
                   date, NULL, # setting restart_date to NULL
                   epoch_dates, pars$info, assumptions)
names(baseline) <- regions

message("Writing parameters_info.csv")
write_csv(pars$info, "parameters_info.csv")
message("Writing parameters_proposal.csv")
write_csv(pars$proposal, "parameters_proposal.csv")
message("Writing parameters_prior.csv")
write_csv(pars$prior, "parameters_prior.csv")

message("Writing parameters_base.rds")
saveRDS(baseline, "parameters_base.rds")

message("Writing parameters_transform.R")
fs::file_copy("R/transform.R",
              "parameters_transform.R", overwrite = TRUE)

message("Printing supplementary figures")
png("fig_sup_vacc_age.png", units = "in", width = 6, height = 6, res = 300)
supl_fig_vac_age()
dev.off()
