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

initial(time) <- 0
initial(S) <- N - I0
initial(R) <- 0
initial(I) <- I0
initial(cases_cumul) <- 0
initial(cases_inc) <- 0

beta <- user(0.2)
gamma <- user(0.1)
I0 <- user(10)
N <- user(1000)

freq <- user(4)
dt <- 1.0 / freq

exp_noise <- user(1e6)
cases <- data()

modelled_cases <- cases_inc + rexp(exp_noise)
compare(cases) ~ poisson(modelled_cases)
