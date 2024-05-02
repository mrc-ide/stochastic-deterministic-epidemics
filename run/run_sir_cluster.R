
short_run <- TRUE

## Create SIR model and SIR data simulation
orderly2::orderly_run("sir_data")

## ---------------------------
## Run in the cluster
## ---------------------------

## 1. Basic cluster setup
hipercow::hipercow_init(driver = "windows")
hipercow::hipercow_provision()

## Fit SIR model

## 2. Run fits ----

det_fit <- hipercow::task_create_expr(
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = TRUE,
                                          adaptive_proposal = FALSE)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 4)
)
det_result <- hipercow::task_result(det_fit)


det_adap_fit <- hipercow::task_create_expr(
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = TRUE,
                                          adaptive_proposal = TRUE)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 4)
)
det_adap_result <- hipercow::task_result(det_adap_fit)


stoch_fit <- hipercow::task_create_expr(
  orderly2::orderly_run('sir_fits',
                        parameters = list(short_run = short_run,
                                          deterministic = FALSE,
                                          adaptive_proposal = FALSE)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)
stoch_result <- hipercow::task_result(stoch_fit)


## 3. Run particle filter on SIR samples

filter_test <- hipercow::task_create_expr(
  orderly2::orderly_run('sir_filter_test',
                        parameters = list(short_run = short_run)),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)
filter_test_result <- hipercow::task_result(filter_test)


## 4. Create plots from SIR fits and particle filter samples
orderly2::orderly_run("sir_plots",
                      list(short_run = short_run))
