orderly2::orderly_parameters(short_run = TRUE,
                             n_seeds = 1L,
                             skip_filter_test = FALSE)

for (seed in seq_len(n_seeds)) {
  orderly2::orderly_dependency(
    "sir_filter_test",
    quote(latest(parameter:short_run == this:short_run && parameter:data_seed == environment:seed && parameter:skip_filter_test == this:skip_filter_test)),
    c("inputs/${seed}/deterministic_fit.rds" = "deterministic_fit.rds",
      "inputs/${seed}/stochastic_fit.rds" = "stochastic_fit.rds",
      "inputs/${seed}/filter_data.rds" = "outputs/filter_data.rds",
      "inputs/${seed}/filter_samples.rds" = "outputs/filter_samples.rds"))  
}

orderly2::orderly_artefact("Plots",
                           c(paste0("figs/combined_parameter_correlation_heatmap_", seq_len(n_seeds), ".png"),
                             paste0("figs/parameter_correlation_", seq_len(n_seeds), ".png"),
                             paste0("figs/model_fit_prevalence_", seq_len(n_seeds), ".png"),
                             paste0("figs/model_fit_incidence_", seq_len(n_seeds), ".png")))

if (!skip_filter_test) {
  orderly2::orderly_artefact("Filter test plots",
                             c(paste0("figs/particle_filter_", seq_len(n_seeds), ".png"),
                               paste0("figs/filter_samples.png_", seq_len(n_seeds), ".png")))  
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
  ggsave(filename = paste0("figs/parameter_correlation_", data_seed, ".png"), 
         plot = plot_parameter_correlation_ggplot(dat[[data_seed]]), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/combined_parameter_correlation_heatmap_", data_seed, ".png"), 
         plot = plot_combined_parameter_correlation_heatmap(dat[[data_seed]]), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/model_fit_prevalence_", data_seed, ".png"), 
         plot = plot_sir_model(dat[[data_seed]], "I"), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  ggsave(filename = paste0("figs/model_fit_incidence_", data_seed, ".png"), 
         plot = plot_sir_model(dat[[data_seed]], "cases_inc"), 
         bg = "white", width = 15, height = 9, dpi = 200)
  
  if (!skip_filter_test) {
    ## THIS IS A TERRIBLE FIX. ENSURE THAT stoch_filtered_sample in the filter_test is renamed
    names(filter_samples)[2] <- "stoch_filtered_samples"
    
    sample_pars_index <- which.max(stoch_sir_fit$pars_full[, "gamma"])
    param_values <- seq(from = 0,
                        to = det_sir_fit$pars_full[sample_pars_index, ]["beta"] * 3,
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
}



