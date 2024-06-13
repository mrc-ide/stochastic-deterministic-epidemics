
short_run <- FALSE
skip_filter_test <- TRUE
n_seeds <- 5L
fit_lambda <- FALSE
recoveries_data <- FALSE

## Create SIR model and SIR data simulation
for (data_seed in seq_len(n_seeds)) {
  orderly2::orderly_run("sir_data",
                        parameters = list(data_seed = data_seed))
}

## Fit SIR model

# Deterministic
for (data_seed in seq_len(n_seeds)) {
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = TRUE,
                                          data_seed = data_seed,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data))
}
  

# Stochastic
for (data_seed in seq_len(n_seeds)) {
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = FALSE,
                                          data_seed = data_seed,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data))
}


## Run particle filter on SIR samples
for (data_seed in seq_len(n_seeds)) {
  orderly2::orderly_run('sir_filter_test',
                        parameters = list(short_run = short_run,
                                          data_seed = data_seed,
                                          skip_filter_test = skip_filter_test,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data))  
}


## Create plots from SIR fits and particle filter samples
orderly2::orderly_run("sir_plots",
                      list(short_run = short_run,
                           n_seeds = n_seeds,
                           skip_filter_test = skip_filter_test,
                           fit_lambda = fit_lambda,
                           recoveries_data = recoveries_data))
