## 1. sarscov2_data
orderly2::orderly_run("sarscov2_data")

## 2. sarscov2_parameters 
orderly2::orderly_run("sarscov2_parameters", 
                      parameters = list(deterministic = TRUE,
                                        assumptions = "central"))

orderly2::orderly_run("sarscov2_parameters", 
                      parameters = list(deterministic = FALSE,
                                        assumptions = "central"))

## ---------------------------
## Run in the cluster
## ---------------------------

## Basic cluster setup
hipercow::hipercow_init(driver = "windows")
hipercow::hipercow_provision()

regions <- sircovid::regions("england")

short_run <- FALSE

#----

## Deterministic
det_fits <- 
  hipercow::task_create_bulk_expr(
    orderly2::orderly_run('sarscov2_fits',
                          parameters = list(region = region,
                                            short_run = short_run,
                                            deterministic = TRUE,
                                            assumptions = "central")),
    data.frame(region = regions),
    resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                             cores = 8))

res <- hipercow::hipercow_bundle_result(det_fits$name)

det_combined <- hipercow::task_create_expr(
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(short_run = short_run,
                                          deterministic = TRUE,
                                          assumptions = "central")),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)
res <- hipercow::task_result(det_combined)


## Stochastic
stoch_fits <- 
  hipercow::task_create_bulk_expr(
    orderly2::orderly_run('sarscov2_fits',
                          parameters = list(region = region,
                                            short_run = short_run,
                                            deterministic = FALSE,
                                            assumptions = "central")),
    data.frame(region = regions),
    resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                             cores = 32))

res <- hipercow::hipercow_bundle_result(stoch_fits$name)

stoch_combined <- hipercow::task_create_expr(
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(short_run = short_run,
                                          deterministic = FALSE,
                                          assumptions = "central")),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 32)
)
res <- hipercow::task_result(stoch_combined)
