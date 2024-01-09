compute_severity <- function(pars, severity_data) {
  dt <- 0.25 # TODO: tidy this up
  expected <- c("mu_D", "mu_D_2", "p_G_D", "p_G_D_CHR", "p_H",
                "p_H_2", "p_H_3", "p_H_D", "p_ICU", "p_ICU_2",
                "p_ICU_D", "p_W_D")
  stopifnot(all(expected %in% names(pars)))
  ## WARNING: this is a hack
  list2env(pars[expected], environment())

  # a. Probability of death in different hospital compartments
  p_D_date <- sircovid::sircovid_date(c("2020-04-01", "2020-06-01",
                                        "2020-10-01", "2020-12-15",
                                        "2021-01-15", "2021-02-01"))
  mu_D_vec <- c(1, mu_D, mu_D, mu_D_2, mu_D_2, mu_D)
  p_ICU_D_value <- p_ICU_D * mu_D_vec # ICU
  p_H_D_value <- p_H_D * mu_D_vec # General ward (never triaged for ICU)
  p_W_D_value <- p_W_D * mu_D_vec # Stepdown care (gen ward, after ICU)

  # b. Probability of ICU admission
  p_ICU_date <- sircovid::sircovid_date(c("2020-04-01", "2020-06-01"))
  p_ICU_value <- c(p_ICU, p_ICU_2)

  # c. Probablity of hospitalisation
  p_H_date <- sircovid::sircovid_date(
    c("2020-10-01", "2020-12-15", "2021-02-15"))
  p_H_value <- c(p_H, p_H_2, p_H_3)

  # d. Probability of being admitted with positive PCR
  p_star_date <- sircovid::sircovid_date(c("2020-03-15", "2020-07-01",
                                           "2020-09-20", "2021-06-27"))
  p_star_value <- c(0.1, 0.42, 0.2, 0.45)

  # e. Simplifying assumptions
  severity_data[severity_data$Name == "p_sero_pos_1", 2:18] <- 0.85
  severity_data[severity_data$Name == "p_sero_pos_2", 2:18] <- 0.85
  severity_data[severity_data$Name == "p_G_D", 2:18] <- 0.05

  severity <- sircovid::lancelot_parameters_severity(
    dt,
    severity_data,
    p_H = list(value = p_H_value, date = p_H_date),
    p_ICU = list(value = p_ICU_value, date = p_ICU_date),
    p_ICU_D = list(value = p_ICU_D_value, date = p_D_date),
    p_H_D = list(value = p_H_D_value, date = p_D_date),
    p_W_D = list(value = p_W_D_value, date = p_D_date),
    p_G_D = list(value = p_G_D),
    p_G_D_CHR = list(value = p_G_D_CHR),
    p_star = list(value = p_star_value, date = p_star_date))

  severity
}


compute_progression <- function(pars, progression_data) {
  dt <- 0.25 # TODO: make this flexible
  expected <- c("mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4")
  stopifnot(all(expected %in% names(pars)))
  ## WARNING: this is a hack
  list2env(pars[expected], environment())

  k_parameters <-
    progression_data[grep("^k_", progression_data$parameter), ]
  gammas <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "value"]
  names(gammas) <-
    progression_data[grep("^gamma_", progression_data$parameter),
                     "parameter"]
  gammas <- as.list(gammas)

  # Reduce length of stay; same dates apply
  mu_gamma_H_date <- sircovid::sircovid_date(c("2020-12-01",
                                               "2021-01-01",
                                               "2021-03-01",
                                               "2021-06-01",
                                               "2021-09-01"))
  mu_gamma_H_value <- c(1, 1 / mu_gamma_H, 1 / mu_gamma_H_2, 1 / mu_gamma_H_3,
                        1 / mu_gamma_H_4)

  gamma_E <- gammas$gamma_E
  gamma_ICU_pre <- gammas$gamma_ICU_pre
  gamma_ICU_D <- gammas$gamma_ICU_D
  gamma_ICU_W_D <- gammas$gamma_ICU_W_D
  gamma_ICU_W_R <- gammas$gamma_ICU_W_R
  gamma_H_R_value <- gammas$gamma_H_R * mu_gamma_H_value
  gamma_H_D_value <- gammas$gamma_H_D * mu_gamma_H_value
  gamma_W_R_value <- gammas$gamma_W_R * mu_gamma_H_value
  gamma_W_D_value <- gammas$gamma_W_D * mu_gamma_H_value
  
  # Time to diagnosis if admitted without test
  gamma_U_value <- 1 / 3

  gamma_PCR_pre_value <- 0.1922243
  gamma_PCR_pos_value <- 0.083
  
  progression <- sircovid::lancelot_parameters_progression(
    dt,
    gamma_E = list(value = gamma_E),
    gamma_ICU_pre = list(value = gamma_ICU_pre),
    gamma_H_D = list(value = gamma_H_D_value, date = mu_gamma_H_date),
    gamma_H_R = list(value = gamma_H_R_value, date = mu_gamma_H_date),
    gamma_ICU_D = list(value = gamma_ICU_D),
    gamma_ICU_W_D = list(value = gamma_ICU_W_D),
    gamma_ICU_W_R = list(value = gamma_ICU_W_R),
    gamma_W_D = list(value = gamma_W_D_value, date = mu_gamma_H_date),
    gamma_W_R = list(value = gamma_W_R_value, date = mu_gamma_H_date),
    gamma_U = list(value = gamma_U_value),
    gamma_PCR_pre = list(value = gamma_PCR_pre_value),
    gamma_PCR_pos = list(value = gamma_PCR_pos_value)
    )
  progression[k_parameters$parameter] <- k_parameters$value

  ## These could possibly be moved to the sircovid as defaults
  progression$k_sero_pre_1 <- 1
  progression$gamma_sero_pre_1 <- 1 / 13
  progression$k_sero_pre_2 <- 1
  progression$gamma_sero_pre_2 <- 1 / 13
  progression$k_PCR_pre <- 1
  progression$k_PCR_pos <- 1
  progression$k_sero_pos_1 <- 1
  progression$gamma_sero_pos_1 <- 1 / 200
  progression$k_sero_pos_2 <- 1
  progression$gamma_sero_pos_2 <- 1 / 400

  progression
}


compute_observation <- function(pars, model_type, pillar2_age_bands, region) {

  if (model_type == "BB") {
    expected_model_specific <- c("p_NC", "p_NC_weekend",
                                 "rho_pillar2_tests")
  } else {
    expected_model_specific <- c("phi_pillar2_cases",
                                 "phi_pillar2_cases_weekend",
                                 paste0("phi_pillar2_cases_", pillar2_age_bands),
                                 paste0("phi_pillar2_cases_weekend_", pillar2_age_bands),
                                 "alpha_pillar2_cases")
  }
  expected <- c("alpha_D", "alpha_H", "alpha_death_hosp", expected_model_specific)
  stopifnot(all(expected %in% names(pars)))
  ## WARNING: this is a hack
  list2env(pars[expected], environment())

  observation <- sircovid::lancelot_parameters_observation()

  if (model_type == "BB") {
    observation[paste0("p_NC_", pillar2_age_bands)] <- p_NC
    observation[paste0("p_NC_weekend_", pillar2_age_bands)] <- p_NC_weekend
    observation$rho_pillar2_tests <- rho_pillar2_tests
  } else {
    if (region %in% sircovid::regions("england")) {
      ## Not fitting to under 15s, just give them the same value as 15 to 24
      observation$phi_pillar2_cases_under15 <- phi_pillar2_cases_15_24
      observation$phi_pillar2_cases_15_24 <- phi_pillar2_cases_15_24
      observation$phi_pillar2_cases_25_49 <- phi_pillar2_cases_25_49
      observation$phi_pillar2_cases_50_64 <- phi_pillar2_cases_50_64
      observation$phi_pillar2_cases_65_79 <- phi_pillar2_cases_65_79
      observation$phi_pillar2_cases_80_plus <- phi_pillar2_cases_80_plus
      observation$phi_pillar2_cases_weekend_under15 <- phi_pillar2_cases_weekend_15_24
      observation$phi_pillar2_cases_weekend_15_24 <- phi_pillar2_cases_weekend_15_24
      observation$phi_pillar2_cases_weekend_25_49 <- phi_pillar2_cases_weekend_25_49
      observation$phi_pillar2_cases_weekend_50_64 <- phi_pillar2_cases_weekend_50_64
      observation$phi_pillar2_cases_weekend_65_79 <- phi_pillar2_cases_weekend_65_79
      observation$phi_pillar2_cases_weekend_80_plus <- phi_pillar2_cases_weekend_80_plus
    } else {
      observation[paste0("phi_pillar2_cases_", pillar2_age_bands)] <- phi_pillar2_cases
      observation[paste0("phi_pillar2_cases_weekend_", pillar2_age_bands)] <- phi_pillar2_cases_weekend
    }
    observation$kappa_pillar2_cases <- 1 / pars[["alpha_pillar2_cases"]]
  }

  ## kappa for hospital data streams (not all will actually be used)
  observation$kappa_ICU <- 1 / alpha_H
  observation$kappa_general <- 1 / alpha_H
  observation$kappa_hosp <- 1 / alpha_H
  observation$kappa_admitted <- 1 / alpha_H
  observation$kappa_diagnoses <- 1 / alpha_H
  observation$kappa_all_admission <- 1 / alpha_H
  observation$kappa_death_hosp <- 1 / alpha_death_hosp

  ## kappa for death data streams (not all will actually be used)
  observation$kappa_death_carehomes <- 1 / alpha_D
  observation$kappa_death_comm <- 1 / alpha_D
  observation$kappa_death_non_hosp <- 1 / alpha_D
  observation$kappa_death <- 1 / alpha_D

  observation
}


apply_assumptions <- function(baseline, assumptions) {
  ## TODO: validation here that we are using a valid set of
  ## assumptions, more important as we add more things that depend on
  ## this.
  
  stopifnot(assumptions %in% names(baseline$rel_severity_delta_omicron))
  baseline$rel_severity_delta_omicron <-
    baseline$rel_severity_delta_omicron[[assumptions]]

  stopifnot(assumptions %in% names(baseline$cross_immunity_delta_omicron))
  baseline$cross_immunity_delta_omicron <-
    baseline$cross_immunity_delta_omicron[[assumptions]]
  
  stopifnot(assumptions %in% names(baseline$severity_cross_multiplier_omicron))
  baseline$severity_cross_multiplier_omicron <-
    baseline$severity_cross_multiplier_omicron[[assumptions]]
  
  stopifnot(assumptions %in% names(baseline$natural_waning_rate))
  baseline$natural_waning_rate <- baseline$natural_waning_rate[[assumptions]]

  baseline
}


## This will get simplified considerably once we drop the previous
## two-stage fitting; that will be needed to bring in 3 and 4 stage
## fitting really.
make_transform <- function(baseline) {
  
  expected <- c("date", "model_type", "region", "restart_date", "epoch_dates",
                "beta_date", "beta_names", "pillar2_age_bands",
                "severity_data", "progression_data",
                "sens_and_spec", "initial_seed_size", "initial_seed_pattern",
                "natural_waning_rate", "cross_immunity_wildtype_alpha",
                "cross_immunity_alpha_delta", "cross_immunity_delta_omicron", 
                "rel_gamma_alpha", "rel_gamma_alpha_delta",
                "rel_gamma_delta_omicron",
                ## Lots of vaccination things
                "rel_severity_wildtype_alpha",
                "rel_severity_alpha_delta",
                "rel_severity_delta_omicron",
                "severity_cross_multiplier_delta",
                "severity_cross_multiplier_omicron",
                "vaccine_eligibility_min_age",
                "vaccine_progression_rate",
                "vaccine_schedule",
                "vaccine_schedule_effect",
                "vaccine_uptake",
                "vaccine_mean_days_between_doses",
                "vaccine_index_dose2",
                "vaccine_index_booster",
                "vaccine_days_to_effect",
                #Thom adds
                "prop_pfizer",
                "average_vacc_efficacy_alpha",
                "average_vacc_efficacy_delta",
                #
                "n_doses",
                "strain_seed_size",
                "strain_seed_pattern")
  stopifnot(setequal(expected, names(baseline)))

  epoch_dates <- baseline$epoch_dates

  ## WARNING: vaccine_eligibility_min_age and vaccine_uptake (probably
  ## akong others) are not actually used here because we use a
  ## schedule that has been built.  So these exist only so that they
  ## can be used in onward simulations that do not update these
  ## parameters (such as the MTPs).  This will be tidied up once we
  ## support partial parameter updating.

  expected <- c("eps", "m_CHR", "m_CHW", "start_date", baseline$beta_names,
                ## severity
                "mu_D", "mu_D_2", "p_G_D", "p_G_D_CHR", "p_H",
                "p_H_2", "p_H_3", "p_H_D", "p_ICU", "p_ICU_2",
                "p_ICU_D", "p_W_D",
                ## progression
                "mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4",
                # multistrain
                "ta_delta", "seed_date_delta",
                "rel_p_H_delta", "rel_p_ICU_delta", "rel_p_D_delta",
                ## observation
                c("alpha_D", "alpha_H", "alpha_death_hosp"),
                if (baseline$model_type == "BB")
                  c("p_NC", "p_NC_weekend",
                      "rho_pillar2_tests") else
                if (baseline$model_type == "NB")
                    c("phi_pillar2_cases", "phi_pillar2_cases_weekend",
                      paste0("phi_pillar2_cases_", baseline$pillar2_age_bands),
                      paste0("phi_pillar2_cases_weekend_",
                             baseline$pillar2_age_bands),
                      "alpha_pillar2_cases"))

  `%||%` <- function(a, b) if (is.null(a)) b else a

  function(pars) {
    stopifnot(setequal(expected, names(pars)))
    beta_value <- unname(pars[baseline$beta_names])
    pars <- as.list(pars) # using list access below

    progression <- compute_progression(pars, baseline$progression_data)
    severity <- compute_severity(pars, baseline$severity_data)
    observation <- compute_observation(pars, baseline$model_type,
                                       baseline$pillar2_age_bands,
                                       baseline$region)

    vaccine_schedule <- baseline$vaccine_schedule_effect
    
    stage_parameters <- function(strains, vaccine_doses) {
      
      if (strains == "Alpha") {
        
        cross_immunity <- 1
        strain_seed_date <- NULL
        strain_seed_size <- NULL
        strain_seed_pattern <- NULL
        strain_transmission <- 1
        rel_severity <- lapply(
          baseline$rel_severity_alpha_delta, function(x) x[, 1, , drop = FALSE])
        strain_rel_gamma <- baseline$rel_gamma_alpha
        
        strain_rel_p_hosp_if_sympt <- 1
        strain_rel_p_icu <- 1
        strain_rel_p_death <- 1
        
      } else if (strains == "Alpha_Delta") {
        
        cross_immunity <- baseline$cross_immunity_alpha_delta
        strain_seed_date <- pars$seed_date_delta
        strain_seed_size <- baseline$strain_seed_size
        strain_seed_pattern <- baseline$strain_seed_pattern
        strain_transmission <- c(1, pars$ta_delta)
        rel_severity <- baseline$rel_severity_alpha_delta
        strain_rel_gamma <- baseline$rel_gamma_alpha_delta
        
        ## Values for: Alpha, Delta, Alpha -> Delta, Delta -> Alpha
        strain_rel_p_hosp_if_sympt <-
          c(1, pars$rel_p_H_delta, 
            pars$rel_p_H_delta * baseline$severity_cross_multiplier_delta, 1)
        strain_rel_p_icu <-
          c(1, pars$rel_p_ICU_delta, pars$rel_p_ICU_delta, 1)
        strain_rel_p_death <-
          c(1, pars$rel_p_D_delta, pars$rel_p_D_delta, 1)
        
      } else {
        stop("'strains' input not supported")
      }
      
      if (vaccine_doses == 0) {
        
        rel_severity <- lapply(
          rel_severity, function(x) 1)
        vaccine_progression_rate <- 0
        vaccine_schedule <- NULL
        vaccine_index_dose2 <- NULL
        vaccine_index_booster <- NULL
        ## this is just the default value
        n_doses <- 2L
        
      } else if (vaccine_doses == 2) {
        
        rel_severity <- lapply(
          rel_severity, function(x) x[, , 1:4, drop = FALSE])
        vaccine_progression_rate <- baseline$vaccine_progression_rate[1:4]
        vaccine_schedule$doses <- vaccine_schedule$doses[, 1:2, , drop = FALSE]
        vaccine_schedule$n_doses <- 2L
        vaccine_index_dose2 <- baseline$vaccine_index_dose2
        vaccine_index_booster <- NULL
        n_doses <- 2L
        
      } else {
        stop("vaccine_doses must be 0 or 2")
      }
      
      sircovid::lancelot_parameters(
        start_date = pars$start_date,
        region = baseline$region,
        beta_date = baseline$beta_date,
        beta_value = beta_value,
        
        severity = severity,
        progression = progression,
        observation = observation,
        sens_and_spec = baseline$sens_and_spec,
        
        initial_seed_size = baseline$initial_seed_size,
        initial_seed_pattern = baseline$initial_seed_pattern,
        
        eps = pars$eps,
        m_CHW = pars$m_CHW,
        m_CHR = pars$m_CHR,
        
        strain_transmission = strain_transmission, 
        strain_seed_date = strain_seed_date, 
        strain_seed_size = strain_seed_size,
        strain_seed_pattern = strain_seed_pattern,
        
        strain_rel_p_hosp_if_sympt = strain_rel_p_hosp_if_sympt,
        strain_rel_p_icu = strain_rel_p_icu,
        strain_rel_p_death = strain_rel_p_death,
        strain_rel_p_G_D = strain_rel_p_death,
        rel_susceptibility = rel_severity$rel_susceptibility,
        rel_p_sympt = rel_severity$rel_p_sympt,
        rel_p_hosp_if_sympt = rel_severity$rel_p_hosp_if_sympt,
        rel_p_death = rel_severity$rel_p_death,
        rel_infectivity = rel_severity$rel_infectivity,
        
        strain_rel_gamma_E = strain_rel_gamma$E,
        strain_rel_gamma_A = strain_rel_gamma$A,
        strain_rel_gamma_P = strain_rel_gamma$P,
        strain_rel_gamma_C_1 = strain_rel_gamma$C_1,
        strain_rel_gamma_C_2 = strain_rel_gamma$C_2,
        
        vaccine_progression_rate = vaccine_progression_rate,
        vaccine_schedule = vaccine_schedule,
        vaccine_index_dose2 = vaccine_index_dose2,
        vaccine_index_booster = vaccine_index_booster,
        n_doses = n_doses,
        
        waning_rate = baseline$natural_waning_rate,
        cross_immunity = cross_immunity)
      
    }
    
    p1 <- stage_parameters("Alpha", 0)
    p2 <- stage_parameters("Alpha", 2)
    p3 <- stage_parameters("Alpha_Delta", 2)
    
    epochs <- list(
      mcstate::multistage_epoch(
        epoch_dates[1], p2, sircovid::inflate_state_vacc_classes),
      mcstate::multistage_epoch(
        epoch_dates[2], p3, sircovid::inflate_state_strains))
    mcstate::multistage_parameters(p1, epochs = epochs)
  }
}
