orderly2::orderly_parameters(region = "london", deterministic = TRUE, short_run = TRUE, assumptions = "central")

orderly2::orderly_shared_resource(util.R = "util.R")

orderly2::orderly_dependency(
  "sarscov2_data",
  "latest",
  c("data/england_region_data.csv" = "outputs/england_region_data.csv",
    "data/serology.csv" = "outputs/serology_for_inference.csv"))
orderly2::orderly_dependency(
  "sarscov2_parameters",
  "latest(parameter:assumptions == this:assumptions && parameter:deterministic == this:deterministic)",
  c("parameters/base.rds" = "parameters_base.rds",
    "parameters/info.csv" = "parameters_info.csv",
    "parameters/prior.csv" = "parameters_prior.csv",
    "parameters/proposal.csv" = "parameters_proposal.csv",
    "parameters/transform.R" = "parameters_transform.R"))

orderly2::orderly_artefact("pMCMC trace plots", "outputs/pmcmc_traceplots.pdf")
orderly2::orderly_artefact("pMCMC trace plots multipage", "outputs/pmcmc_traceplots_separate.pdf")
orderly2::orderly_artefact("PMCMC results for combined task", "outputs/fit.rds")

library(sircovid)
library(spimalot)
library(tidyr)

orderly2::orderly_resource("data.R")
orderly2::orderly_resource("support.R")
source("data.R")
source("support.R")

source("util.R")

version_check("sircovid", "0.15.1")
version_check("spimalot", "0.8.30")

date <- "2022-02-24"

## We're effectively NOT trimming any data stream as backfill is not an issue here
trim_deaths <- 0
trim_pillar2 <- 0
trim_pillar2_date <- FALSE
adm_backfill_date <- date

## MCMC control (only applies if short_run = FALSE)
if (deterministic) {
  burnin <- 5000
  n_mcmc <- 30000
  n_sample <- 1000
  chains <- 8
  kernel_scaling <- 0.1
} else {
  burnin <- 1000
  n_mcmc <- 5000
  n_sample <- 1000
  chains <- 4
  kernel_scaling <- 0.2
}

region <- spimalot::spim_check_region(region, multiregion = FALSE)

pars <- spimalot::spim_fit_pars_load("parameters", region, assumptions,
                                     kernel_scaling)
pars <- simplify_transform(pars, "parameters", date)

## Fix all unused parameters
## (parameters not impacting fitting before the date parameter)
pars <- fix_unused_parameters(pars, date)

restart_date <- readRDS("parameters/base.rds")[[region[[1]]]]$restart_date

## NOTE: only currently using compiled compare for the deterministic
## model as it has a bigger increase in speed here.
control <- spimalot::spim_control(
  short_run, chains, deterministic = deterministic,
  multiregion = FALSE, severity = TRUE, demography = TRUE,
  date_restart = restart_date, adaptive_proposal = deterministic,
  n_mcmc = n_mcmc, burnin = burnin, n_sample = n_sample,
  compiled_compare = deterministic)

data_rtm <- read_csv("data/england_region_data.csv")
data_serology <- read_csv("data/serology.csv")

data <- spim_data(
  date, region, data_rtm, data_serology, trim_deaths, trim_pillar2,
  adm_backfill_date, trim_pillar2_date, full_data = FALSE)

filter <- spimalot::spim_particle_filter(data, pars$mcmc,
                                         control$particle_filter,
                                         deterministic)

## To run the model at this point, we just need to run:
##
## > filter$run(pars$mcmc$model(pars$mcmc$initial()))

## This bit takes ages, of course
samples <- spimalot::spim_fit_run(pars, filter, control$pmcmc)

## This is the data set including series that we do not fit to, and
## with the full series of carehomes deaths.
data_full <- spim_data(
  date, region, data_rtm, data_serology,
  trim_deaths, trim_pillar2,
  adm_backfill_date, trim_pillar2_date, full_data = TRUE)

## This is new, and used only in sorting out the final outputs. Some
## explanation would be useful.
data_inputs <- list(rtm = data_rtm,
                    full = data_full,
                    fitted = data)

dat <- spimalot::spim_fit_process(samples, pars, data_inputs,
                                  control$particle_filter)
dat <- add_full_proposal(dat, pars)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat$fit, "outputs/fit.rds")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$fit$samples),
  width = 16, height = 9)

write_pdf(
  "outputs/pmcmc_traceplots_separate.pdf",
  spimalot::spim_plot_fit_traces_separate(dat$fit$samples),
  width = 16, height = 9)
