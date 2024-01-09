orderly2::orderly_parameters(date = NULL,
                             pars_set = NULL,
                             assumptions = "central",
                             deterministic = FALSE,
                             adaptive_proposal = FALSE)

orderly2::orderly_resource(c("pars/central/bb/2020/deterministic/info.csv",
                             "pars/central/bb/2020/deterministic/proposal.csv",
                             "pars/central/bb/2020/deterministic_adaptive/info.csv",
                             "pars/central/bb/2020/deterministic_adaptive/proposal.csv",
                             "pars/central/bb/2020/stochastic/info.csv",
                             "pars/central/bb/2020/stochastic/proposal.csv",
                             "pars/central/bb/2021/deterministic/info.csv",
                             "pars/central/bb/2021/deterministic/proposal.csv",
                             "pars/central/bb/2021/deterministic_adaptive/info.csv",
                             "pars/central/bb/2021/deterministic_adaptive/proposal.csv",
                             "pars/central/bb/2021/stochastic/info.csv",
                             "pars/central/bb/2021/stochastic/proposal.csv",
                             "vaccine_data/vaccine_efficacy_alpha.csv",
                             "vaccine_data/vaccine_efficacy_delta.csv",
                             "vaccine_data/vaccine_efficacy_omicron.csv",
                             "vaccine_data/vaccine_uptake.csv"))

orderly2::orderly_dependency(
  "sarscov2_data",
  "latest",
  c(data_vaccination.csv = "data/data_vaccination.csv",
    weighted_prior_ranges.csv = "data/weighted_prior_ranges.csv"))

orderly2::orderly_artefact(
  "fitted hyperparameters for priors",
  c("parameters_base.rds",
    "parameters_info.csv",
    "parameters_prior.csv",
    "parameters_proposal.csv",
    "parameters_transform.R"))

library(sircovid)
library(spimalot)
library(tidyr)
library(dplyr)
library(forcats)
library(magrittr)

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("R/support.R")
orderly2::orderly_resource("R/priors.R")
orderly2::orderly_resource("R/baseline.R")
orderly2::orderly_resource("R/transform.R")
orderly2::orderly_resource("R/vaccine.R")
orderly2::orderly_resource("R/multiregion.R")
source("R/support.R")
source("R/priors.R")
source("R/baseline.R")
source("R/transform.R")
source("R/vaccine.R")
source("R/multiregion.R")
source("util.R")

version_check("sircovid", "0.14.13")
version_check("spimalot", "0.8.24")

epoch_dates <- c("2020-12-07", "2021-03-08")
restart_date <- NULL

multiregion <- FALSE
model_type <- "BB"

if (!adaptive_proposal) adaptive_proposal <- NULL

if (!deterministic && !is.null(adaptive_proposal)) {
    stop("adaptive pMCMC method not yet implemented, please ensure deterministic = TRUE if adaptive_proposal != NULL")
}


## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(model_type, pars_set, adaptive_proposal, assumptions, deterministic,
                             multiregion)

## The baselines are always region-specific
regions <- sircovid::regions("england")
baseline <- lapply(regions, create_baseline,
                   date, model_type, restart_date,
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
