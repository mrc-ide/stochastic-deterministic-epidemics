orderly2::orderly_parameters(short_run = TRUE,
                             n_seeds = 1L,
                             skip_filter_test = FALSE,
                             fit_lambda = FALSE,
                             recoveries_data = FALSE)

if (fit_lambda) {
  fit_lambda_name <- "fitted lambda"
} else {
  fit_lambda_name <- "fixed lambda"
}

if (recoveries_data) {
  recoveries_data_name <- "with recoveries data"
} else {
  recoveries_data_name <- "no recoveries data"
}

for (seed in seq_len(n_seeds)) {
  orderly2::orderly_dependency(
    "sir_filter_test",
    quote(latest(parameter:short_run == this:short_run && parameter:data_seed == environment:seed && parameter:skip_filter_test == this:skip_filter_test  && parameter:fit_lambda == this:fit_lambda && parameter:recoveries_data == this:recoveries_data)),
    c("inputs/${seed}/deterministic_fit.rds" = "deterministic_fit.rds",
      "inputs/${seed}/stochastic_fit.rds" = "stochastic_fit.rds",
      "inputs/${seed}/filter_data.rds" = "outputs/filter_data.rds",
      "inputs/${seed}/filter_samples.rds" = "outputs/filter_samples.rds",
      "figs/${seed}. traceplots_deterministic ${fit_lambda_name} ${recoveries_data_name}.png" = "figs/traceplots_deterministic.png",
      "figs/${seed}. traceplots_stochastic ${fit_lambda_name} ${recoveries_data_name}.png" = "figs/traceplots_stochastic.png"))
}

orderly2::orderly_artefact("Plots",
                           c(paste0("figs/", seq_len(n_seeds), ". ", "parameter_correlation", " ", fit_lambda_name, " ", recoveries_data_name, ".png"),
                             paste0("figs/", seq_len(n_seeds), ". ", "model_fit_prevalence", " ", fit_lambda_name, " ", recoveries_data_name, ".png"),
                             paste0("figs/", seq_len(n_seeds), ". ", "model_fit_incidence", " ", fit_lambda_name, " ", recoveries_data_name, ".png"),
                             paste0("figs/", seq_len(n_seeds), ". ", "model_fit_recoveries", " ", fit_lambda_name, " ", recoveries_data_name, ".png")))

if (!skip_filter_test) {
  orderly2::orderly_artefact("Filter test plots",
                             c(paste0("figs/particle_filter_", seq_len(n_seeds), ".png"),
                               paste0("figs/filter_samples_", seq_len(n_seeds), ".png")))  
}

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("plot.R")
orderly2::orderly_resource("support.R")
source("util.R")
source("plot.R")
source("support.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)

dat <- lapply(seq_len(n_seeds), read_inputs)

# Create directory for figures
dir.create("figs", showWarnings = FALSE)

for (data_seed in seq_len(n_seeds)) {
  ggsave(filename = paste0("figs/", data_seed, ". ", "parameter_correlation", " ", fit_lambda_name, " ", recoveries_data_name, ".png"), 
         plot = plot_parameter_correlation_ggplot(dat[[data_seed]]), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/", data_seed, ". ", "model_fit_prevalence", " ", fit_lambda_name, " ", recoveries_data_name, ".png"), 
         plot = plot_sir_model(dat[[data_seed]], "I"), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/", data_seed, ". ", "model_fit_incidence", " ", fit_lambda_name, " ", recoveries_data_name, ".png"), 
         plot = plot_sir_model(dat[[data_seed]], "cases_inc"), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/", data_seed, ". ", "model_fit_recoveries", " ", fit_lambda_name, " ", recoveries_data_name, ".png"), 
         plot = plot_sir_model(dat[[data_seed]], "recoveries_inc"), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  if (!skip_filter_test) {
    ggsave(filename = paste0("figs/filter_samples_", data_seed, ".png"), 
           plot = plot_filter_samples(dat[[data_seed]]), 
           bg = "white", width = 15, height = 9, dpi = 200)
    
    ggsave(filename = paste0("figs/particle_filter_", data_seed, ".png"), 
           plot = create_particle_filter_plot(dat[[data_seed]]$filter_data), 
           bg = "white", width = 15, height = 9, dpi = 200)
  }  
}
