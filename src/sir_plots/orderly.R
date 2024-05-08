orderly2::orderly_parameters(short_run = TRUE, skip_filter_test = FALSE)

orderly2::orderly_dependency("sir_filter_test",
                             "latest(parameter:short_run == this:short_run)",
                            c("deterministic_fit.rds" = "deterministic_fit.rds",
                              "deterministic_adaptive_fit.rds" = "deterministic_adaptive_fit.rds",
                              "stochastic_fit.rds" = "stochastic_fit.rds",
                              "filter_data.rds" = "outputs/filter_data.rds",
                              "filter_samples.rds" = "outputs/filter_samples.rds"))

orderly2::orderly_artefact("Plots",
                           c("figs/combined_parameter_correlation_heatmap.png",
                             "figs/parameter_correlation.png",
                             "figs/model_fit_prevalence.png",
                             "figs/model_fit_incidence.png"))

if (!skip_filter_test) {
  orderly2::orderly_artefact("Filter test plots",
                             c("figs/particle_filter.png",
                               "figs/filter_samples.png"))  
}

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

true_history <- stochastic_fit$true_history

incidence <- deterministic_fit$data$cases

# Create directory for figures
dir.create("figs", showWarnings = FALSE)

ggsave(filename = "figs/parameter_correlation.png", 
       plot = plot_parameter_correlation_ggplot(stoch_sir_fit, det_sir_fit, 
                                                det_adap_sir_fit), 
       bg = "white", width = 15, height = 9, dpi = 200)

ggsave(filename = "figs/combined_parameter_correlation_heatmap.png", 
       plot = plot_combined_parameter_correlation_heatmap(
         stoch_sir_fit, det_sir_fit, det_adap_sir_fit), 
       bg = "white", width = 15, height = 9, dpi = 200)

ggsave(filename = "figs/model_fit_prevalence.png", 
       plot = plot_sir_model(det_sir_fit, det_adap_sir_fit, stoch_sir_fit, 
                             true_history, incidence, "I"), 
       bg = "white", width = 15, height = 9, dpi = 200)

ggsave(filename = "figs/model_fit_incidence.png", 
       plot = plot_sir_model(det_sir_fit, det_adap_sir_fit, stoch_sir_fit, 
                             true_history, incidence, "cases_inc"), 
       bg = "white", width = 15, height = 9, dpi = 200)

if (!skip_filter_test) {
  ## THIS IS A TERRIBLE FIX. ENSURE THAT stoch_filtered_sample in the filter_test is renamed
  names(filter_samples)[2] <- "stoch_filtered_samples"
  
  sample_pars_index <- which.max(stoch_sir_fit$pars_full[, "gamma"])
  param_values <- seq(from = 0,
                      to = det_adap_sir_fit$pars_full[sample_pars_index, ]["beta"] * 3,
                      length.out = length(filter_samples$det_filtered_samples))
  
  det_df <- create_filter_df(filter_samples, "det", param_values)
  stoch_df <- create_filter_df(filter_samples, "stoch", param_values)
  
  samples_df <- rbind(det_df, stoch_df)
  
  ggsave(filename = "figs/filter_samples.png", 
         plot = plot_filter_samples(samples_df), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = "figs/particle_filter.png", 
         plot = create_particle_filter_plot(filter_data), 
         bg = "white", width = 15, height = 9, dpi = 200)
}


