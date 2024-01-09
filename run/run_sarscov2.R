date <- "2021-09-13"
pars_set <- "2021"
short_run <- TRUE

regions <- sircovid::regions("england")


## Run data task
orderly2::orderly_run("sarscov2_data")

## Run parameters task
pars_deterministic <- orderly2::orderly_run(
  "sarscov2_parameters",
  parameters = list(date = date,
                    pars_set = pars_set,
                    deterministic = TRUE,
                    adaptive_proposal = FALSE,
                    assumptions = "central"))

pars_deterministic_adaptive <- orderly2::orderly_run(
  "sarscov2_parameters",
  parameters = list(date = date,
                    pars_set = pars_set,
                    deterministic = TRUE,
                    adaptive_proposal = TRUE,
                    assumptions = "central"))

pars_stochastic <- orderly2::orderly_run(
  "sarscov2_parameters",
  parameters = list(date = date,
                    pars_set = pars_set,
                    deterministic = FALSE,
                    adaptive_proposal = FALSE,
                    assumptions = "central"))

## Run deterministic fits
fits_deterministic <-
  lapply(X = regions,
         FUN = function(x) {
           orderly2::orderly_run('sarscov2_fits',
                                 parameters = list(region = x,
                                                   date = date,
                                                   pars_set = pars_set,
                                                   short_run = short_run,
                                                   deterministic = TRUE,
                                                   adaptive_proposal = FALSE,
                                                   assumptions = "central"))})

combined_deterministic <-
  orderly2::orderly_run('sarscov2_fits_combined',
                       parameters = list(date = date,
                                         pars_set = pars_set,
                                         short_run = short_run,
                                         deterministic = TRUE,
                                         adaptive_proposal = FALSE,
                                         assumptions = "central"))

## Run deterministic adaptive fits
fits_deterministic_adaptive <-
  lapply(X = regions,
         FUN = function(x) {
           orderly2::orderly_run('sarscov2_fits',
                                 parameters = list(region = x,
                                                   date = date,
                                                   pars_set = pars_set,
                                                   short_run = short_run,
                                                   deterministic = TRUE,
                                                   adaptive_proposal = TRUE,
                                                   assumptions = "central"))})

combined_deterministic_adaptive <-
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(date = date,
                                          pars_set = pars_set,
                                          short_run = short_run,
                                          deterministic = TRUE,
                                          adaptive_proposal = TRUE,
                                          assumptions = "central"))

## Run deterministic adaptive fits
fits_stochastic <-
  lapply(X = regions,
         FUN = function(x) {
           orderly2::orderly_run('sarscov2_fits',
                                 parameters = list(region = x,
                                                   date = date,
                                                   pars_set = pars_set,
                                                   short_run = short_run,
                                                   deterministic = FALSE,
                                                   adaptive_proposal = FALSE,
                                                   assumptions = "central"))})

combined_stochastic <-
  orderly2::orderly_run('sarscov2_fits_combined',
                        parameters = list(date = date,
                                          pars_set = pars_set,
                                          short_run = short_run,
                                          deterministic = FALSE,
                                          adaptive_proposal = FALSE,
                                          assumptions = "central"))
