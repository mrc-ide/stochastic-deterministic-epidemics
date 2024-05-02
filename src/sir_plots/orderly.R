orderly2::orderly_parameters(short_run = TRUE)

orderly2::orderly_dependency("sir_filter_test",
                             "latest(parameter:short_run == this:short_run)",
                            c("deterministic_fit.rds" = "deterministic_fit.rds",
                              "deterministic_adaptive_fit.rds" = "deterministic_adaptive_fit.rds",
                              "stochastic_fit.rds" = "stochastic_fit.rds",
                              "filter_data.rds" = "outputs/filter_data.rds",
                              "filter_samples.rds" = "outputs/filter_samples.rds"))

orderly2::orderly_artefact("Plots",
                           c("figs/create_particle_filter_plot.png",
                             "figs/plot_combined_parameter_correlation_heatmap.png",
                             "figs/plot_filter_samples.png",
                             "figs/plot_parameter_correlation_ggplot.png",
                             "figs/plot_sir_model.png"))

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("support.R")
source("util.R")
source("support.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)

deterministic_fit <- readRDS("deterministic_fit.rds")
deterministic_adaptive_fit <- readRDS("deterministic_adaptive_fit.rds")
stochastic_fit <- readRDS("stochastic_fit.rds")

det_sir_fit <- deterministic_fit$samples
det_adap_sir_fit <- deterministic_adaptive_fit$samples
stoch_sir_fit <- stochastic_fit$samples

filter_data <- readRDS("filter_data.rds")
filter_samples <- readRDS("filter_samples.rds")

## THIS IS A TERRIBLE FIX. ENSURE THAT stoch_filtered_sample in the filter_test is renamed
names(filter_samples)[2] <- "stoch_filtered_samples"

true_history <- stochastic_fit$true_history

incidence <- true_history["cases_inc", 1, ]

sample_pars_index <- which.max(stoch_sir_fit$pars_full[, "gamma"])
param_values <- seq(from = 0,
                    to = det_adap_sir_fit$pars_full[sample_pars_index, ]["beta"] * 3,
                    length.out = length(filter_samples$det_filtered_samples))

det_df <- create_filter_df(filter_samples, "det", param_values)
stoch_df <- create_filter_df(filter_samples, "stoch", param_values)

samples_df <- rbind(det_df, stoch_df)

# Create directory for figures
dir.create("figs", showWarnings = FALSE)

# List of plotting functions and their arguments
plot_functions <- list(
  plot_parameter_correlation_ggplot = list(stoch_sir_fit, det_sir_fit, det_adap_sir_fit),
  plot_combined_parameter_correlation_heatmap = list(stoch_sir_fit, det_sir_fit, det_adap_sir_fit),
  plot_sir_model = list(det_sir_fit, det_adap_sir_fit, stoch_sir_fit, true_history, incidence, "I"),
  plot_filter_samples = list(samples_df),
  create_particle_filter_plot = list(filter_data)
)

# Generate and save each plot
lapply(names(plot_functions), function(func_name) {
  file_name <- paste0("figs/", func_name, ".png")
  generate_and_save_plots(func_name, plot_functions[[func_name]], file_name)
})
