
short_run <- TRUE

## Create SIR model and SIR data simulation
orderly2::orderly_run("data")

## Fit SIR model
runs <- list(
  list(short_run = short_run, deterministic = TRUE, adaptive_proposal = FALSE),
  list(short_run = short_run, deterministic = TRUE, adaptive_proposal = TRUE),
  list(short_run = short_run, deterministic = FALSE, adaptive_proposal = FALSE)
)

lapply(runs, function(params) {
  orderly2::orderly_run("fits", params)
})

## Run particle filter on SIR samples
orderly2::orderly_run("filter_test", 
                        list(short_run = short_run))

## Create plots from SIR fits and particle filter samples
orderly2::orderly_run("plots",
                      list(short_run = short_run))
