source("global_util.R")
source("global_vaccine.R")

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
