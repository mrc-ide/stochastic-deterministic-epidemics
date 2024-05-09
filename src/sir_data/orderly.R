orderly2::orderly_artefact(
  "Model code for downstream compilation",
  "sir.R"
)
orderly2::orderly_artefact(
  "Simulated data and true epidemic history",
  c("outputs/data.rds",
    "outputs/true_history.rds")
)

orderly2::orderly_resource("sir.R")
orderly2::orderly_resource("support.R")
source("support.R")

library(odin.dust)

sir <- odin.dust::odin_dust("sir.R")

dir.create("outputs", FALSE, TRUE)
saveRDS(sir, "outputs/model.rds")

pars <- list(beta = 0.2,
             gamma = 0.1)

mod <- sir$new(pars, 0, 1, seed = 1L)

n_days <- 130

y <- mod$simulate(seq(0, 4 * n_days, by = 4))
rownames(y) <- names(mod$info()$index)
saveRDS(y, "outputs/true_history.rds")

set.seed(1)
cases_model <- y["cases_inc", , seq_len(n_days) + 1]
cases_data <- rpois(n_days, lambda = cases_model + 1e-6)

data <- data.frame(cases = cases_data,
                   day = seq_len(n_days))
saveRDS(data, "outputs/data.rds")


