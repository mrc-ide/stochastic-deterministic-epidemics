# Function to plot parameter correlations
plot_parameter_correlation_ggplot <- function(stoch_sir_fit, det_sir_fit, det_adap_sir_fit) {  
  # Combine the data into a single data frame
  df_stoch <- data.frame(beta = stoch_sir_fit$pars[, 1], gamma = stoch_sir_fit$pars[, 2], Type = "Stochastic")
  df_det <- data.frame(beta = det_sir_fit$pars[, 1], gamma = det_sir_fit$pars[, 2], Type = "Deterministic")
  df_adap <- data.frame(beta = det_adap_sir_fit$pars[, 1], gamma = det_adap_sir_fit$pars[, 2], Type = "Adaptive")
  df <- rbind(df_stoch, df_det, df_adap)
  
  # Plot using ggplot2
  ggplot(df, aes(x = beta, y = gamma, color = Type)) +
    geom_point(alpha = 0.5) + # Translucent points
    scale_color_viridis_d() +  # Colorblind-friendly palette
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Parameter Correlation Plot",
      color = "Model Type",
      x = expression(beta),
      y = expression(gamma)
    ) +
    xlim(0.15, 0.35) +
    ylim(0.05, 0.2)
}

# Function to create a density data frame for heatmap
create_density_df <- function(df, x_bins, y_bins) {
  df %>% 
    mutate(
      x_bin = cut_width(beta, width = x_bins, boundary = 0),
      y_bin = cut_width(gamma, width = y_bins, boundary = 0)
    ) %>% 
    group_by(x_bin, y_bin, Type) %>% 
    summarise(count = n(), .groups = 'drop')
}

# Function to plot combined parameter correlation heatmap
plot_combined_parameter_correlation_heatmap <- function(stoch_sir_fit, det_sir_fit, det_adap_sir_fit) {
  # Combine the data into a single data frame
  df_stoch <- data.frame(beta = stoch_sir_fit$pars[, 1], gamma = stoch_sir_fit$pars[, 2], Type = "Stochastic")
  df_det <- data.frame(beta = det_sir_fit$pars[, 1], gamma = det_sir_fit$pars[, 2], Type = "Deterministic")
  df_adap <- data.frame(beta = det_adap_sir_fit$pars[, 1], gamma = det_adap_sir_fit$pars[, 2], Type = "Adaptive")
  df <- rbind(df_stoch, df_det, df_adap)

  # Plot using ggplot2 with density heatmap
  ggplot(df, aes(x = beta, y = gamma, fill = ..density..)) +
    stat_density_2d(geom = "raster", contour = FALSE) +
    scale_fill_viridis_c(begin = 0, end = 1) +  # Adjust color scale to increase brightness at low end
    facet_wrap(~ Type) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.position = "right",
      panel.grid.major = element_blank(),  # Remove grid lines
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Combined Parameter Space Exploration Heatmap",
      fill = "Density",
      x = expression(beta),
      y = expression(gamma)
    ) +
    xlim(0.15, 0.35) +
    ylim(0.05, 0.2)
}

# Function to process simulation results for state plot
process_sim_results <- function(sim_data, name_prefix, state) {
  res <- sim_data$trajectories$state[state, , ]
  y <- cbind(colMeans(res),
             t(apply(res, 2, quantile, probs = c(0.025, 0.975))))
  colnames(y) <- c(paste0(name_prefix, "_mean"),
                   paste0(name_prefix, "_lb"),
                   paste0(name_prefix, "_ub"))
  return(y)
}

# Function to plot SIR states
plot_sir_model <- function(det_sir_fit, det_adap_sir_fit, stoch_sir_fit, true_history, incidence, state) {
  # Process simulation results
  y1 <- process_sim_results(det_sir_fit, "det", state)
  y2 <- process_sim_results(det_adap_sir_fit, "adap", state)
  y3 <- process_sim_results(stoch_sir_fit, "stoch", state)
  
  # Add true history data to the stochastic results
  y3 <- cbind(y3, data = true_history[state, , ])
  
  # Combine all results
  y <- cbind(y1, y2, y3)
  y <- as.data.frame(y)
  
  # Prepare date for x-axis
  date <- true_history["time", , ]
    # Use viridis palette for colorblind-friendly colors
  viridis_colors <- viridis(4)
  fit_cols <- setNames(viridis_colors, c("det_mean", "adap_mean", "stoch_mean", "data"))

  if (state == "I") {
    ylab <- "Infection prevalence"
  } else if (state == "cases_inc") {
    ylab <- "Infection incidence"
  }
  
  # Create the plot with consistent aesthetics
  ggplot(y, aes(x = date)) +
    geom_line(aes(y = det_mean, color = "det_mean"), size = 1) +
    geom_ribbon(aes(ymin = det_lb, ymax = det_ub, fill = "det_mean"), alpha = 0.25) +
    geom_line(aes(y = adap_mean, color = "adap_mean"), size = 1) +
    geom_ribbon(aes(ymin = adap_lb, ymax = adap_ub, fill = "adap_mean"), alpha = 0.25) +
    geom_line(aes(y = stoch_mean, color = "stoch_mean"), size = 1) +
    geom_ribbon(aes(ymin = stoch_lb, ymax = stoch_ub, fill = "stoch_mean"), alpha = 0.25) +
    geom_point(aes(y = data, color = "truth"), size = 1.5) +
    scale_color_manual(values = fit_cols) +
    scale_fill_manual(values = fit_cols) +
    labs(
      y = ylab,
      x = "Time (Days)",
      title = "SIR Model States",
      color = "Type",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank(),
      axis.line = element_blank()  # Remove the solid axis lines
    )
}

# Function to create particle filter scatter plot with boxplot for publication
create_particle_filter_plot <- function(filter_data) {
  # Prepare data
  colnames(filter_data$deterministic_filtered) <- filter_data$n_particles
  df_det <- as.data.frame(filter_data$deterministic_filtered) %>%
    mutate(variable = "MLP_deterministic_sample") %>%
    pivot_longer(!variable, names_to = "n_particles")

  colnames(filter_data$stochastic_filtered) <- filter_data$n_particles
  df_stoch <- as.data.frame(filter_data$stochastic_filtered) %>%
    mutate(variable = "MLP_stochastic_sample") %>%
    pivot_longer(!variable, names_to = "n_particles")

  df_combined <- rbind(df_det, df_stoch) %>%
    mutate(n_particles = as.numeric(n_particles), n_particles = factor(n_particles))

  # Plot
  ggplot(df_combined, aes(y = value, x = n_particles, color = variable, fill = variable)) +
    geom_point(position = position_jitterdodge(jitter.width = .25, dodge.width = 0.6), size = 1, pch = 19) +
    stat_boxplot(geom = "errorbar", width = .25, position = position_dodge(0.6), color = "#00000055") +
    geom_boxplot(width = .25, position = position_dodge(0.6), outlier.shape = NA, alpha = 0.5, color = "black") +
    scale_color_manual(values = c("MLP_deterministic_sample" = "#0000ff55", "MLP_stochastic_sample" = "#ff000055")) +
    scale_fill_manual(values = c("MLP_deterministic_sample" = "#0000ff55", "MLP_stochastic_sample" = "#ff000055")) +
    theme_classic() +
    theme(legend.position = c(.8, .1), panel.grid.major = element_line(), panel.grid.minor = element_line()) +
    labs(fill = "Fitted pipeline", color = "Fitted pipeline", title = "Particle-Filter Attributable Stochasticity") +
    ylab("MLP")
}

# Function to create a data frame for each filter type
create_filter_df <- function(filter_samples, filter_type, param_values) {
  data.frame(
    Beta = param_values,
    LogPosterior = filter_samples[[paste0(filter_type, "_filtered_samples")]],
    Filter = filter_type
  )
}

# Function to plot the data, excluding the first value for each filter type
plot_filter_samples <- function(samples_df) {
  # Filter out the first value for each filter type
  filtered_df <- samples_df %>% 
                 group_by(Filter) %>% 
                 slice(-1)  # Remove the first row for each group

  ggplot(filtered_df, aes(x = Beta, y = LogPosterior, group = Filter, color = Filter)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c("det" = "blue", "stoch" = "red")) +
    labs(
      x = expression(beta),
      y = "Log-posterior",
      title = "Parameter Log-posterior Value Across Fixed Particle Size (Particles = 192)",
      color = "Filter Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_blank(),
      axis.line = element_blank()
    )
}
