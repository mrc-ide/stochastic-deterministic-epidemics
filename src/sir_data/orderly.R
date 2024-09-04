
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

library(odin2)

sir <- odin2::odin("sir.R")

dir.create("outputs", FALSE, TRUE)

pars <- list(beta = 0.2,
             gamma = 0.1,
             lambda = 10)

sys <- dust2::dust_system_create(sir(), pars, n_particles = 1,
                                 seed = data_seed, dt = 0.25)
dust2::dust_system_set_state_initial(sys)

## run until number of infectives goes to 0
t <- 1
y <- dust2::dust_system_simulate(sys, c(0, t))
index_I <- 3
inf_zero <- y[index_I, t + 1] == 0
while (!inf_zero) {
  t <- t + 1
  y <- cbind(y, dust2::dust_system_simulate(sys, t))
  inf_zero <- y[index_I, t + 1] == 0
}
rownames(y) <- c("S", "R", "I", "cases_cumul", "cases_inc", "recoveries_inc")
saveRDS(y, "outputs/true_history.rds")

n_days <- t
set.seed(data_seed)
model_cases <- y["cases_inc", seq_len(n_days) + 1]
data_cases <- rpois(n_days, lambda = model_cases + rexp(n_days, 1e6))

model_recoveries <- y["recoveries_inc", seq_len(n_days) + 1]
data_recoveries <- rpois(n_days, lambda = model_recoveries + rexp(n_days, 1e6))

data <- data.frame(cases = data_cases,
                   recoveries = data_recoveries,
                   time = seq_len(n_days))
saveRDS(data, "outputs/data.rds")
