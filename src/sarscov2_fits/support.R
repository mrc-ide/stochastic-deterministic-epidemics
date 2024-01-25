simplify_transform <- function(pars, path, date) {
  
  e <- new.env()
  sys.source(file.path(path, "transform.R"), e)
  
  make_transform <- e$make_transform
  pars$transform <- make_transform(pars$base, date)
  
  pars$mcmc <- spimalot:::spim_pars_mcmc_single(pars$info, pars$prior, 
                                                pars$proposal, pars$transform)
  
  
  pars$base$epoch_dates <-
    pars$base$epoch_dates[pars$base$epoch_dates <= sircovid_date(date)]
  keep_strain_epochs <- pars$base$strain_epochs <= length(pars$base$epoch_dates)
  pars$base$strain_epochs <- pars$base$strain_epochs[keep_strain_epochs]
  pars$base$date <- date
  
  pars
}


fix_unused_parameters <- function(pars, date) {
  
  ## Automatically detect which betas to fix
  ## We need to keep all betas up to the first one with date greater than or 
  ## equal to the date parameter
  i <- max(which(pars$base$beta_date < sircovid::sircovid_date(date)))
  beta_fixed <- setdiff(pars$base$beta_names, sprintf("beta%d", seq_len(i + 1)))
  
  ## Now we will fix other parameters that have no impact before the date
  ## parameter
  
  ## Note firstly that the following parameters are required whatever
  ## the date parameter is:
  ## "alpha_admission", "alpha_D", "alpha_death_hosp", "alpha_H", "p_G_D","p_H",
  ## "p_H_D", "p_ICU", "p_ICU_D", "p_W_D", "start_date"
  
  ## Now we declare the date from which these parameters have an impact
  pars_dates <- list(
    ## Various changepoint parameters. We must include them the day after the
    ## previous changepoint
    mu_D = "2020-04-02",
    mu_D_2 = "2020-09-16",
    mu_D_3 = "2020-12-02",
    mu_D_4 = "2021-02-05",
    mu_D_5 = "2021-11-05",
    mu_gamma_H = "2020-12-02",
    mu_gamma_H_2 = "2021-01-02",
    mu_gamma_H_3 = "2021-03-02",
    mu_gamma_H_4 = "2021-06-02",
    p_G_D_2 = "2020-05-02",
    p_H_2 = "2021-11-05",
    p_ICU_2 = "2020-04-02",
    
    ## Pillar 2 parameters - we start fitting pillar 2 from 2020-06-18
    ## This is a Thursday, so first weekend day is 2020-06-20
    p_NC_15_24 = "2020-06-18",
    p_NC_25_49 = "2020-06-18",
    p_NC_50_64 = "2020-06-18",
    p_NC_65_79 = "2020-06-18",
    p_NC_80_plus = "2020-06-18",
    p_NC_weekend_15_24 = "2020-06-20",  
    p_NC_weekend_25_49 = "2020-06-20", 
    p_NC_weekend_50_64 = "2020-06-20", 
    p_NC_weekend_65_79 = "2020-06-20", 
    p_NC_weekend_80_plus = "2020-06-20", 
    rho_pillar2_tests = "2020-06-18",
    
    ## alpha parameters
    rel_p_D_alpha = "2020-09-17",
    rel_p_H_alpha = "2020-09-17",
    rel_p_ICU_alpha = "2020-09-17",
    seed_date_alpha = "2020-09-17",
    ta_alpha = "2020-09-17",
    
    ## delta parameters
    rel_p_D_delta = "2021-03-08",
    rel_p_H_delta = "2021-03-08",
    rel_p_ICU_delta = "2021-03-08",
    seed_date_delta = "2021-03-08",
    ta_delta = "2021-03-08",
    
    ## omicron parameters
    rel_p_D_omicron = "2021-11-01",
    rel_p_H_omicron = "2021-11-01",
    rel_p_ICU_omicron = "2021-11-01",
    seed_date_omicron = "2021-11-01",
    ta_omicron = "2021-11-01"
  )

  ## Fix parameters
  fixed <- c(beta_fixed, names(pars_dates)[pars_dates > date])
  pars$mcmc <- pars$mcmc$fix(pars$mcmc$initial()[fixed])
  
  pars
}

add_full_proposal <- function(dat, pars) {
  
  new_prop <- dat$fit$parameters$proposal
  old_prop <- pars$proposal
  
  new_prop_pars <- new_prop$name
  
  full_prop <- data.frame(0 * old_prop)
  full_prop[new_prop_pars, new_prop_pars] <- new_prop[, -c(1, 2)]
  full_prop <- cbind(data.frame(region = dat$fit$samples$info$region,
                                name = rownames(full_prop)),
                     full_prop)
  rownames(full_prop) <- NULL
  dat$fit$parameters$proposal <- full_prop
  
  dat
}
