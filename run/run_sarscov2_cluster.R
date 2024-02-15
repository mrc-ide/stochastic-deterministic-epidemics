## 1. sarscov2_data
orderly2::orderly_run("sarscov2_data")

## 2. sarscov2_parameters 
orderly2::orderly_run("sarscov2_parameters", 
                      parameters = list(deterministic = TRUE,
                                        assumptions = "central"))

## ---------------------------
## Run in the cluster
## ---------------------------

## 1. Basic cluster setup
hipercow::hipercow_init(driver = "windows")
hipercow::hipercow_provision()

regions <- sircovid::regions("england")

#----

## 2. Short runs ----
fits <- 
  hipercow::task_create_bulk_call(
    function(x) {
      orderly2::orderly_run('sarscov2_fits',
                            parameters = list(region = x,
                                              short_run = TRUE,
                                              deterministic = TRUE,
                                              assumptions = "central"))},
    regions,
    resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                             cores = 8))
batch <- fits$name

## Collect results
res <- hipercow::hipercow_bundle_result(batch)

# Combined
combined <- hipercow::task_create_expr(
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(short_run = TRUE,
                                          deterministic = TRUE,
                                          assumptions = "central")),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 8)
)
combined_result <- hipercow::task_result(combined)

#----

## 3. Long runs ----
fits <- 
  hipercow::task_create_bulk_call(
    function(x) {
      orderly2::orderly_run('sarscov2_fits',
                            parameters = list(region = x,
                                              short_run = FALSE,
                                              deterministic = TRUE,
                                              assumptions = "central"))},
    regions,
    resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                             cores = 8))
batch <- fits$name

## Collect results
res <- hipercow::hipercow_bundle_result(batch)

# Combined
# Combined
combined <- hipercow::task_create_expr(
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(short_run = FALSE,
                                          deterministic = TRUE,
                                          assumptions = "central")),
  resources = hipercow::hipercow_resources(queue = 'AllNodes',
                                           cores = 8)
)
combined_result <- hipercow::task_result(combined)
