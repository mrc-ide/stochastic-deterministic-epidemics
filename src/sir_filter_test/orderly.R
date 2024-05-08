
orderly2::orderly_parameters(short_run = TRUE, skip_filter_test = FALSE)

orderly2::orderly_dependency("sir_fits",
                             'latest(parameter:short_run == this:short_run && parameter:deterministic == FALSE && parameter:adaptive_proposal == FALSE)',
                             c("stochastic_fit.rds" = "outputs/fit.rds",
                               "sir.R" = "sir.R"))

orderly2::orderly_dependency("sir_fits",
                             "latest(parameter:short_run == this:short_run && parameter:deterministic == TRUE && parameter:adaptive_proposal == FALSE)",
                             c("deterministic_fit.rds" = "outputs/fit.rds"))

orderly2::orderly_dependency("sir_fits",
                             "latest(parameter:short_run == this:short_run && parameter:deterministic == TRUE && parameter:adaptive_proposal == TRUE)",
                             c("deterministic_adaptive_fit.rds" = "outputs/fit.rds"))

orderly2::orderly_artefact("Fit objects for downstream usage",
                           c("stochastic_fit.rds",
                             "deterministic_fit.rds",
                             "deterministic_adaptive_fit.rds"))

orderly2::orderly_artefact("Filter outputs",
                           c("outputs/filter_data.rds",
                             "outputs/filter_samples.rds"))

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("support.R")
source("util.R")
source("support.R")

library(mcstate)
library(odin.dust)

stochastic_fit <- readRDS("stochastic_fit.rds")
deterministic_adaptive_fit <- readRDS("deterministic_adaptive_fit.rds")

det_adap_sir_fit <- deterministic_adaptive_fit$samples
stoch_sir_fit <- stochastic_fit$samples

det_adap_sir_pars <- deterministic_adaptive_fit$pars
stoch_sir_pars <- stochastic_fit$pars

det_adap_pars <- det_adap_sir_fit$pars_full
stoch_pars <- stoch_sir_fit$pars_full

if (short_run) {
  n_filter_iter <- 10
} else {
  n_filter_iter <- 1000
}

dir.create("outputs", FALSE, TRUE)

if (skip_filter_test) {
  saveRDS(NULL, "outputs/filter_samples.rds")
  saveRDS(NULL, "outputs/filter_data.rds")
} else {
  ## -----------------------------------------------------------------------------
  ## Part 6 Varying single parameter (sample) value across fixed particle filter size
  ## -----------------------------------------------------------------------------
  
  ## Build stochastic filter
  set.seed(1234)
  
  data <- stoch_sir_fit$predict$filter$data
  sir <- odin.dust::odin_dust("sir.R")
  compare <- stoch_sir_fit$predict$filter$compare
  index <- stoch_sir_fit$predict$filter$index
  n_particles <- 192
  
  n_threads <- control_cores()
  
  filter_stochastic <-
    mcstate::particle_filter$new(data, sir, n_particles = n_particles,
                                 compare = NULL, index = index, seed = 1L,
                                 n_threads = n_threads)
  
  #sample_pars_index <- sample(nrow(stoch_pars), 1)
  sample_pars_index <- which.max(stoch_pars[, "gamma"])
  ## parameters to index and test 
  det_adap_pars <- det_adap_pars[sample_pars_index, ]
  stoch_pars <- stoch_pars[sample_pars_index, ]
  param_test <- "beta"
  
  ## sample values near true value (deterministic and stochastic outputs)
  ## number of values to try
  n_try <- 50
  
  ## parameter values at either side of the fitted parameter value
  param_values <- seq(from = 0,
                      to = det_adap_pars[param_test] * 3,
                      length.out = n_try)
  
  ## run filter on single parameter
  ## (sample around "true" deterministic parameter value)
  ## calculate the log-posterior (likelihood + prior)
  det_filtered_samples <- vapply(seq_along(param_values),
                                 function(i) {
                                   det_adap_pars[param_test] <- param_values[i]
                                   filter_stochastic$run(stoch_sir_pars$model(det_adap_pars)) +
                                     stoch_sir_pars$prior(det_adap_pars)
                                 },
                                 numeric(1)
  )
  
  
  
  ## run filter on single parameter
  ## (sample around "true" stochastic parameter value)
  ## calculate the log-posterior (likelihood + prior)
  stoch_filtered_sample <- vapply(seq_along(param_values),
                                  function(i) {
                                    stoch_pars[param_test] <- param_values[i]
                                    filter_stochastic$run(stoch_sir_pars$model(stoch_pars)) + 
                                      stoch_sir_pars$prior(stoch_pars)
                                  },
                                  numeric(1)
  )
  
  stoch_filtered_sample <- unlist(stoch_filtered_sample)
  
  ## save filter outputs
  filter_samples <- list(det_filtered_samples, stoch_filtered_sample)
  names(filter_samples) <- c("det_filtered_samples", "stoch_filtered_sample")
  saveRDS(filter_samples, "outputs/filter_samples.rds")
  
  ## -----------------------------------------------------------------------------
  ## Part 7 Varying particle filter size across parameter sample
  ## -----------------------------------------------------------------------------
  
  ## define the number of particle sizes to test
  ## the number of particle filter iterations
  det_adap_pars <- det_adap_sir_fit$pars
  stoch_pars <- stoch_sir_fit$pars
  n_particles <- c(1,192,1024) #2 ^ c(6:10)
  filter_data <- list(NULL)
  
  run_n_particles <- function(n, pars) {
    filter <- mcstate::particle_filter$new(data, sir,
                                           n_particles = n,
                                           compare = compare, index = index,
                                           seed = 1L, n_threads = n_threads)
    
    vapply(seq_len(nrow(pars)),
           function(i) {
             filter$run(
               stoch_sir_pars$model(
                 pars[i, ]
               )
             ) + stoch_sir_pars$prior(pars[i, ])
           },
           numeric(1))
  }
  
  deterministic_filtered <-
    vapply(n_particles, function (i) run_n_particles(i, det_adap_pars),
           numeric(nrow(det_adap_pars)))
  
  stochastic_filtered <-
    vapply(n_particles, function (i) run_n_particles(i, stoch_pars),
           numeric(nrow(stoch_pars)))
  
  filter_data <- list(deterministic_filtered = deterministic_filtered,
                      stochastic_filtered = stochastic_filtered,
                      n_particles = n_particles)
  
  saveRDS(filter_data, "outputs/filter_data.rds")  
}
