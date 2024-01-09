orderly2::orderly_parameters(region = "london",
                             date = NULL,
                             pars_set = NULL,
                             deterministic = FALSE,
                             adaptive_proposal = FALSE,
                             short_run = FALSE,
                             assumptions = "central")

orderly2::orderly_artefact(
  "Fitting outputs",
  "outputs/fit.rds",
  "parameters_base.rds",
  "parameters_info.csv",
  "parameters_prior.csv",
  "parameters_proposal.csv",
  "parameters_transform.R"
)
orderly2::orderly_artefact(
  "Traceplots",
  "outputs/pmcmc_traceplots.pdf"
)

orderly2::orderly_dependency(
  "sarscov2_data",
  "latest",
  c("data/rtm.csv" = "data/uk_rtm.csv",
    "data/serology.csv" = "data/serology.csv"))
orderly2::orderly_dependency(
  "sarscov2_parameters",
  "latest(parameter:date == this:date && parameter:pars_set == this:pars_set && parameter:assumptions == this:assumptions && parameter:deterministic == this:deterministic && parameter:adaptive_proposal == this:adaptive_proposal)",
  c("parameters/base.rds" = "parameters_base.rds",
    "parameters/info.csv" = "parameters_info.csv",
    "parameters/prior.csv" = "parameters_prior.csv",
    "parameters/proposal.csv" = "parameters_proposal.csv",
    "parameters/transform.R" = "parameters_transform.R"))

library(sircovid)
library(spimalot)
library(dplyr)
library(tidyr)

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("data.R")
source("data.R")
source("util.R")

version_check("sircovid", "0.14.13")
version_check("spimalot", "0.8.24")

if (!deterministic && adaptive_proposal) {
    stop("adaptive pMCMC method not yet implemented, please ensure deterministic = TRUE if adaptive_proposal != NULL")
}

model_type <- "BB"

trim_deaths <- 0
trim_pillar2 <- 0

## MCMC control (only applies if short_run = FALSE)
burnin <- 5000
n_mcmc <- 15000
chains <- 4

total_pars <- read.csv("parameters/info.csv")

kernel_scaling <- 2.38^2/sum(total_pars$initial > 1e-08)

region <- spimalot::spim_check_region(region, multiregion = FALSE)

pars <- spimalot::spim_fit_pars_load("parameters", region, assumptions,
                                     kernel_scaling)

restart_date <- readRDS("parameters/base.rds")[[region[[1]]]]$restart_date

## This will probably want much more control, if we are to support
## using rrq etc to create a multinode job; some of that will depend a
## bit on the combination of multiregion and deterministic I think

#adaptive_proposal is set to default proposal which needs to be edited by calling
#the pmcmc object inside the control i.e control$pmcmc$adaptive_proposal(...)

control <- spimalot::spim_control(
  short_run, chains, deterministic = deterministic, date_restart = restart_date,
  compiled_compare = deterministic, n_mcmc = n_mcmc, burnin = burnin,
  n_particles = 192, adaptive_proposal = adaptive_proposal)

if (adaptive_proposal) {
# c(initial_scaling, scaling_increment, initial_weight)
control$pmcmc$adaptive_proposal[c(1,2,4)] <- c(
  kernel_scaling, kernel_scaling*.02, 100) 
}

data_rtm <- read_csv("data/rtm.csv")
data_serology <- read_csv("data/serology.csv")

data <- spim_data(
  date, region, model_type, data_rtm, data_serology,
  trim_deaths, trim_pillar2,
  full_data = FALSE)


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
  date, region, model_type, data_rtm,
  data_serology, trim_deaths, trim_pillar2,
  full_data = TRUE)

## This is new, and used only in sorting out the final outputs. Some
## explanation would be useful.
data_inputs <- list(rtm = data_rtm,
                    full = data_full,
                    fitted = data)

dat <- spimalot::spim_fit_process(samples, pars, data_inputs,
    control$particle_filter)


dir.create("outputs", FALSE, TRUE)
saveRDS(dat$fit, "outputs/fit.rds")
#saveRDS(dat$restart, "outputs/restart.rds")

#dir.create("outputs/parameters", FALSE, TRUE)
#spimalot::spim_pars_pmcmc_save(dat$fit$parameters, "outputs/parameters")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$fit$samples),
  width = 16, height = 9)
