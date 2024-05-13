read_inputs <- function(data_seed) {
  list(
    deterministic_fit = 
      readRDS(paste0("inputs/", data_seed, "/deterministic_fit.rds")),
    stochastic_fit =
      readRDS(paste0("inputs/", data_seed, "/stochastic_fit.rds")),
    filter_data = 
      readRDS(paste0("inputs/", data_seed, "/filter_data.rds")),
    filter_samples = 
      readRDS(paste0("inputs/", data_seed, "/filter_samples.rds"))
  )
}  
  