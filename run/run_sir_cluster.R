
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

## ---------------------------
## Run in the cluster
## ---------------------------

## 1. Basic cluster setup
hipercow::hipercow_init(driver = "windows")
hipercow::hipercow_provision()

## Fit SIR model

## 2. Run fits ----

det_fit <- hipercow::task_create_bulk_expr(
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = TRUE,
                                          data_seed = data_seed,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data)),
  data.frame(data_seed = seq_len(n_seeds)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 4)
)
det_result <- hipercow::hipercow_bundle_result(det_fit$name)

stoch_fit <- hipercow::task_create_bulk_expr(
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = FALSE,
                                          data_seed = data_seed,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data)),
  data.frame(data_seed = seq_len(n_seeds)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)
stoch_result <- hipercow::hipercow_bundle_result(stoch_fit$name)


## 3. Run particle filter on SIR samples

filter_test <- hipercow::task_create_bulk_expr(
  orderly2::orderly_run('sir_filter_test',
                        parameters = list(short_run = short_run,
                                          data_seed = data_seed,
                                          skip_filter_test = skip_filter_test,
                                          fit_lambda = fit_lambda,
                                          recoveries_data = recoveries_data)),
  data.frame(data_seed = seq_len(n_seeds)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)

filter_test_result <- hipercow::hipercow_bundle_result(filter_test$name)


## 4. Create plots from SIR fits and particle filter samples
orderly2::orderly_run("sir_plots",
                      list(short_run = short_run,
                           n_seeds = n_seeds,
                           skip_filter_test = skip_filter_test,
                           fit_lambda = fit_lambda,
                           recoveries_data = recoveries_data))

