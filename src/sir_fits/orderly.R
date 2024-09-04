orderly2::orderly_parameters(short_run = TRUE,
                             deterministic = TRUE,
                             data_seed = 1L,
                             fit_lambda = FALSE,
                             recoveries_data = FALSE)

orderly2::orderly_dependency("sir_data", 
                             "latest(parameter:data_seed == this:data_seed)",
                            c("data.rds" = "outputs/data.rds",
                              "true_history.rds" = "outputs/true_history.rds",
                              "sir.R" = "sir.R"))

orderly2::orderly_artefact("Traceplots",
                           "figs/traceplots.png")
orderly2::orderly_artefact("Fit object",
                           "outputs/fit.rds")
orderly2::orderly_artefact(
  "Model code for downstream compilation",
  "sir.R"
)

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("plot.R")
orderly2::orderly_resource("support.R")
source("util.R")
source("plot.R")
source("support.R")

library(odin2)
library(dust2)
library(monty)

sir <- odin2::odin("sir.R")

data <- readRDS("data.rds")
true_history <- readRDS("true_history.rds")

if (!recoveries_data) {
  data$recoveries <- NA
}

dt <- 0.25

if (fit_lambda) {
  prior <- monty::monty_dsl({
    beta ~ Uniform(0, 1)
    gamma ~ Uniform(0, 1)
    lambda ~ Uniform(0, 1000)
  })
  
  sir_packer <- monty::monty_packer(c("beta", "gamma", "lambda"))
  
  proposal_vcv <- rbind(c(0.00057, 0.00034, 0), 
                           c(0.00034, 0.00026, 0),
                           c(0, 0, 5))
  
} else {
  prior <- monty::monty_dsl({
    beta ~ Uniform(0, 1)
    gamma ~ Uniform(0, 1)
  })
  
  sir_packer <- monty::monty_packer(c("beta", "gamma"), 
                                    fixed = list(lambda = 10))
  
  proposal_vcv <- rbind(c(0.00057, 0.00034), c(0.00034, 0.00026))
}

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
  filter <- dust2::dust_unfilter_create(sir(), 0, data, n_particles = 1)
} else {
  max_workers <- 4
  pos <- seq_len(max_workers)
  n_workers <- max(pos[n_threads_total %% pos == 0 & pos <= n_chains])
  n_threads <- n_threads_total / n_workers
  
  n_particles <- 1024
  # rerun_every <- 100
  filter <- dust2::dust_filter_create(sir(), 0, data, n_particles = n_particles,
                                      n_threads = n_threads)
  
}

likelihood <- dust2::dust_filter_monty(filter, sir_packer)

posterior <- prior + likelihood

sampler <- monty::monty_sampler_random_walk(proposal_vcv)

runner <- monty::monty_runner_parallel(n_workers)

samples <- monty::monty_sample(posterior, sampler, n_steps, 
                               n_chains = n_chains,
                               initial = sir_packer$pack(list(beta = 0.2,
                                                              gamma = 0.1,
                                                              lambda = 10)))

fit <- list(samples = samples,
            data = data,
            true_history = true_history)

dir.create("outputs", FALSE, TRUE)
saveRDS(fit, "outputs/fit.rds")

dir.create("figs", FALSE, TRUE)
write_png("figs/traceplots.png", width = 3000, height = 1800, res = 200,
          plot_traceplots(fit$samples))
