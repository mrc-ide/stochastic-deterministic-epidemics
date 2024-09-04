p_SI <- 1 - exp(-beta * dt * I / N)
p_IR <- 1 - exp(-gamma * dt)
n_IR <- Binomial(I, p_IR)
n_SI <- Binomial(S, p_SI)

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(cases_cumul) <- cases_cumul + n_SI
update(cases_inc) <- if (time %% 1 == 0) n_SI else cases_inc + n_SI
update(recoveries_inc) <- if (time %% 1 == 0) n_IR else recoveries_inc + n_IR
# update(cases_inc) <- cases_inc + n_SI
# update(recoveries_inc) <- recoveries_inc + n_IR

initial(S) <- N - min(I0, N)
initial(R) <- 0
initial(I) <- min(I0, N)
initial(cases_cumul) <- 0
initial(cases_inc) <- 0
initial(recoveries_inc) <- 0
#initial(cases_inc, zero_every = 1) <- 0
#initial(recoveries_inc, zero_every = 1) <- 0

I0 <- Poisson(lambda)

beta <- parameter(0.2)
gamma <- parameter(0.1)
lambda<- parameter(10)
N <- parameter(1000)

exp_noise <- parameter(1e6)
cases <- data()
recoveries <- data()


## Likelihood
modelled_cases <- cases_inc + Exponential(rate = exp_noise)
cases ~ Poisson(modelled_cases)

modelled_recoveries <- recoveries_inc + Exponential(rate = exp_noise)
recoveries ~ Poisson(modelled_recoveries)
