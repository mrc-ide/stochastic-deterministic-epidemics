
orderly2::orderly_parameters(data_seed = 1L)

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

pars <- list(beta = 0.2,
             gamma = 0.1)

mod <- sir$new(pars, 0, 1, seed = data_seed)

## run until number of infectives goes to 0
t <- 1
y <- mod$simulate(c(0, 4))
index_I <- which(names(mod$info()$index) == "I")
inf_zero <- y[index_I, , t + 1] == 0
while (!inf_zero) {
  t <- t + 1
  y <- abind::abind(y, mod$simulate(4 * t), along = 3)
  inf_zero <- y[index_I, , t + 1] == 0
}
rownames(y) <- names(mod$info()$index)
saveRDS(y, "outputs/true_history.rds")

n_days <- t
set.seed(data_seed)
cases_model <- y["cases_inc", , seq_len(n_days) + 1]
cases_data <- rpois(n_days, lambda = cases_model + rexp(n_days, 1e6))

data <- data.frame(cases = cases_data,
                   day = seq_len(n_days))
saveRDS(data, "outputs/data.rds")
