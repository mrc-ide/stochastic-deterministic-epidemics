p_SI <- 1 - exp(-beta * dt * I / N)
p_IR <- 1 - exp(-gamma * dt)
n_IR <- rbinom(I, p_IR)
n_SI <- rbinom(S, p_SI)

update(time) <- (step + 1) * dt
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(cases_cumul) <- cases_cumul + n_SI
update(cases_inc) <- if (step %% freq == 0) n_SI else cases_inc + n_SI
update(recoveries_inc) <- if (step %% freq == 0) n_IR else recoveries_inc + n_IR

initial(time) <- 0
initial(S) <- N - min(I0, N)
initial(R) <- 0
initial(I) <- min(I0, N)
initial(cases_cumul) <- 0
initial(cases_inc) <- 0
initial(recoveries_inc) <- 0

I0 <- rpois(lambda)

beta <- user(0.2)
gamma <- user(0.1)
lambda<- user(10)
N <- user(1000)

freq <- user(4)
dt <- 1.0 / freq

exp_noise <- user(1e6)
cases <- data()
recoveries <- data()

modelled_cases <- cases_inc + rexp(exp_noise)
compare(cases) ~ poisson(modelled_cases)

modelled_recoveries <- recoveries_inc + rexp(exp_noise)
compare(recoveries) ~ poisson(modelled_recoveries)
