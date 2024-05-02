orderly2::orderly_parameters(short_run = TRUE,
                            deterministic = TRUE,
                            adaptive_proposal = TRUE)

orderly2::orderly_dependency("sir_data", "latest()",
                            c("data.rds" = "outputs/data.rds",
                              "true_history.rds" = "outputs/true_history.rds",
                              "sir.R" = "sir.R"))

orderly2::orderly_artefact("Traceplots",
                           "figs/traceplots.png")
orderly2::orderly_artefact("Fit object",
                           "outputs/fit.rds")

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("plot.R")
orderly2::orderly_resource("support.R")
source("util.R")
source("plot.R")
source("support.R")

library(mcstate)
library(odin.dust)

if (!deterministic && adaptive_proposal) {
    stop("adaptive pMCMC method not yet implemented, please ensure deterministic = TRUE if adaptive_proposal != NULL")
}

sir <- odin.dust::odin_dust("sir.R")

incidence <- readRDS("data.rds")
true_history <- readRDS("true_history.rds")

dt <- 0.25
data <- mcstate::particle_filter_data(incidence,
  initial_time = 0,
  time = "day",
  rate = 1 / dt)

proposal_kernel <- rbind(c(0.00057, 0.00034), c(0.00034, 0.00026))

pars <- mcstate::pmcmc_parameters$new(
  list(mcstate::pmcmc_parameter("beta", 0.2, min = 0, max = 1,
                                prior = function(p) log(1e-10)),
       mcstate::pmcmc_parameter("gamma", 0.1, min = 0, max = 1,
                                prior = function(p) log(1e-10))),
  proposal = proposal_kernel)

if (short_run) {
  burnin <- 1
  n_steps <- 20
  n_sample <- 10
  n_chains <- 4
} else {
  burnin <- 5000
  n_steps <- 15000
  n_sample <- 10000
  n_chains <- 4
}
n_steps_retain <- ceiling(n_sample / n_chains)

n_threads_total <- control_cores()
if (deterministic) {
  n_workers <- min(n_chains, n_threads_total)
  n_threads <- n_threads_total / n_workers
  p <- mcstate::particle_deterministic$new(data, sir, compare = NULL, index = index,
                                            n_threads = n_threads)
} else {
  max_workers <- 4
  pos <- seq_len(max_workers)
  n_workers <- max(pos[n_threads_total %% pos == 0 & pos <= n_chains])
  n_threads <- n_threads_total / n_workers
  
  n_particles <- 192
  p <- mcstate::particle_filter$new(data, sir, n_particles = n_particles,
                                    compare = NULL, index = index,
                                    n_threads = n_threads)
}
n_threads <- n_threads_total / n_workers

control <- 
  mcstate::pmcmc_control(n_steps = n_steps, n_burnin = burnin,
                         n_threads_total = n_threads_total,
                         n_workers = n_workers, n_chains = n_chains,
                         n_steps_retain = n_steps_retain, save_state = TRUE,
                         adaptive_proposal = adaptive_proposal,
                         save_trajectories = TRUE, progress = TRUE)

samples <-  mcstate::pmcmc(pars, p, control = control)

fit <- list(samples = samples,
            data = data,
            true_history = true_history,
            pars = pars)

dir.create("outputs", FALSE, TRUE)
saveRDS(fit, "outputs/fit.rds")

dir.create("figs", FALSE, TRUE)
write_png("figs/traceplots.png", width = 3000, height = 1800, res = 200,
          plot_traceplots(fit$samples))
