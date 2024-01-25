compute_severity <- function(pars, severity_data, assumptions) {
  dt <- 0.25 # TODO: tidy this up
  expected <- c("mu_D", "mu_D_2", "mu_D_3", "mu_D_4", "mu_D_5",
                "p_G_D",
                "p_G_D_2", "p_H", "p_H_2", "p_H_D", "p_ICU", "p_ICU_2",
                "p_ICU_D", "p_W_D")
  if (assumptions == "mu_d_summer") {
    expected <- c(expected, "mu_D_6")
  }
  stopifnot(all(expected %in% names(pars)))
  ## WARNING: this is a hack
  list2env(pars[expected], environment())

  # a. Probability of death in different hospital compartments
  p_D_date <- sircovid::sircovid_date(c("2020-04-01",
                                        "2020-07-01", "2020-09-15",
                                        "2020-10-15", "2020-12-01",
                                        "2021-02-04",
                                        "2021-04-01", "2021-11-04",
                                        "2021-12-31"))
  mu_D_vec <- c(1, mu_D, mu_D, mu_D_2, mu_D_2, mu_D_3, mu_D_4, mu_D_4, mu_D_5)
  
  if (assumptions == "mu_d_winter") {
    
    p_D_date <- sircovid::sircovid_date(
      c("2020-04-01",
        "2020-07-01", "2020-09-15",
        # Linear increase from  09 to 12
        "2020-12-01",
        "2021-02-04",
        "2021-04-01", "2021-11-04",
        "2021-12-31"))
    
    mu_D_vec <- c(1, mu_D, mu_D, mu_D_2, mu_D_3, mu_D_4, mu_D_4, mu_D_5)
    
  } else if (assumptions == "mu_d_summer") {
    
    p_D_date <- sircovid::sircovid_date(
      c("2020-04-01",
        "2020-07-01", "2020-09-15",
        "2020-10-15", "2020-12-01",
        "2021-02-04",
        # Flat from 04 to 06 then linear from 06 to 11
        "2021-04-01", "2021-06-01",
        "2021-11-04",
        "2021-12-31"))
    
    mu_D_vec <- c(1, mu_D, mu_D, mu_D_2, mu_D_2, mu_D_3, mu_D_4, mu_D_4,
                  mu_D_5, mu_D_6)
  }
  
  p_ICU_D_value <- p_ICU_D * mu_D_vec # ICU
  p_H_D_value <- p_H_D * mu_D_vec # General ward (never triaged for ICU)
  p_W_D_value <- p_W_D * mu_D_vec # Stepdown care (gen ward, after ICU)
  
  # b. Probability of ICU admission
  p_ICU_date <- sircovid::sircovid_date(c("2020-04-01", "2020-06-01"))
  p_ICU_value <- c(p_ICU, p_ICU_2)
  
  # c. Probablity of hospitalisation
  p_H_date <- sircovid::sircovid_date(c("2021-11-04", "2021-12-31"))
  p_H_value <- c(p_H, p_H_2)
  
  # d. Probability of dying outside hospital
  p_G_D_date <- sircovid::sircovid_date(c("2020-05-01", "2020-07-01"))
  p_G_D_value <- c(p_G_D, p_G_D_2)

  # d. Probability of being admitted with positive PCR
  p_star_date <- sircovid::sircovid_date(c("2020-03-15", "2020-07-01",
                                           "2020-09-20", "2021-06-27",
                                           "2021-12-01", "2022-01-01"))
  p_star_value <- c(0.1, 0.42, 0.2, 0.45, 0.45, 0.33)

  # e. Simplifying assumptions
  severity_data[severity_data$Name == "p_sero_pos_1", 2:18] <- 0.85
  severity_data[severity_data$Name == "p_sero_pos_2", 2:18] <- 0.85
  
  severity <- sircovid::lancelot_parameters_severity(
    dt,
    severity_data,
    p_H = list(value = p_H_value, date = p_H_date),
    p_ICU = list(value = p_ICU_value, date = p_ICU_date),
    p_ICU_D = list(value = p_ICU_D_value, date = p_D_date),
    p_H_D = list(value = p_H_D_value, date = p_D_date),
    p_W_D = list(value = p_W_D_value, date = p_D_date),
    p_G_D = list(value = p_G_D_value, date = p_G_D_date),
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
  mu_gamma_H_value <- c(1, 1 / mu_gamma_H, 1 / mu_gamma_H_2,
                        1 / mu_gamma_H_3, 1 / mu_gamma_H_4)

  gamma_E <- gammas$gamma_E
  gamma_ICU_pre <- gammas$gamma_ICU_pre
  gamma_ICU_D <- gammas$gamma_ICU_D
  gamma_ICU_W_D <- gammas$gamma_ICU_W_D
  gamma_ICU_W_R <- gammas$gamma_ICU_W_R
  gamma_H_R_value <- gammas$gamma_H_R * mu_gamma_H_value
  gamma_H_D_value <- gammas$gamma_H_D * mu_gamma_H_value
  gamma_W_R_value <- gammas$gamma_W_R * mu_gamma_H_value
  gamma_W_D_value <- gammas$gamma_W_D * mu_gamma_H_value
  gamma_PCR_pre_value <- 1 / 2
  gamma_PCR_pos_value <- 0.083
  
  # Time to diagnosis if admitted without test
  gamma_U_date <- sircovid::sircovid_date(c("2020-03-15", "2020-04-05",
                                            "2020-09-15", "2020-11-15",
                                            "2020-12-10", "2021-01-20",
                                            "2021-08-01", "2022-03-01"))
  gamma_U_value <- 1 / c(2.1, 1.2, 1.2, 0.7, 1.3, 0.7, 0.4, 0.9)

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
    gamma_U = list(value = gamma_U_value, date = gamma_U_date),
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
  progression$gamma_sero_pos_1 <- 1 / 400
  progression$k_sero_pos_2 <- 1
  progression$gamma_sero_pos_2 <- 1 / 1000

  progression
}


compute_observation <- function(pars, pillar2_age_bands, region) {

  expected <- c("alpha_D", "alpha_H", "alpha_admission", "alpha_death_hosp",
                paste0("p_NC_", pillar2_age_bands),
                paste0("p_NC_weekend_", pillar2_age_bands),
                "rho_pillar2_tests")
  stopifnot(all(expected %in% names(pars)))
  ## WARNING: this is a hack
  list2env(pars[expected], environment())

  observation <- sircovid::lancelot_parameters_observation()

  observation$p_NC_under15 <- p_NC_15_24
  observation$p_NC_15_24 <- p_NC_15_24
  observation$p_NC_25_49 <- p_NC_25_49
  observation$p_NC_50_64 <- p_NC_50_64
  observation$p_NC_65_79 <- p_NC_65_79
  observation$p_NC_80_plus <- p_NC_80_plus
  observation$p_NC_weekend_under15 <- p_NC_weekend_15_24
  observation$p_NC_weekend_15_24 <- p_NC_weekend_15_24
  observation$p_NC_weekend_25_49 <- p_NC_weekend_25_49
  observation$p_NC_weekend_50_64 <- p_NC_weekend_50_64
  observation$p_NC_weekend_65_79 <- p_NC_weekend_65_79
  observation$p_NC_weekend_80_plus <- p_NC_weekend_80_plus
  observation$rho_pillar2_tests <- rho_pillar2_tests

  ## kappa for hospital data streams (not all will actually be used)
  observation$kappa_ICU <- 1 / alpha_H
  observation$kappa_general <- 1 / alpha_H
  observation$kappa_hosp <- 1 / alpha_H
  observation$kappa_admitted <- 1 / alpha_H
  observation$kappa_diagnoses <- 1 / alpha_H
  observation$kappa_all_admission <- 1 / alpha_admission
  observation$kappa_death_hosp <- 1 / alpha_death_hosp

  ## kappa for death data streams (not all will actually be used)
  observation$kappa_death_carehomes <- 1 / alpha_D
  observation$kappa_death_comm <- 1 / alpha_D
  observation$kappa_death_non_hosp <- 1 / alpha_D
  observation$kappa_death <- 1 / alpha_D

  observation
}


apply_assumptions <- function(baseline, assumptions) {
  ## TODO: will use this function again when setting up other sensitivities
  ##       so keeping it as placeholder for now

  baseline
}


## This will get simplified considerably once we drop the previous
## two-stage fitting; that will be needed to bring in 3 and 4 stage
## fitting really.
make_transform <- function(baseline, date = NULL) {
  
  expected <- c("date", "model_type", "region", "assumptions",
                "restart_date", "epoch_dates",
                "intrinsic_severity_dates", "strain_epochs",
                "beta_date", "beta_names", "pillar2_age_bands",
                "severity_data", "progression_data",
                "sens_and_spec", "initial_seed_size", "initial_seed_pattern",
                "natural_waning_rate", "cross_immunity_wildtype",
                "cross_immunity_alpha", "cross_immunity_delta",
                "cross_immunity_omicron",
                "rel_gamma_wildtype",
                "rel_gamma_wildtype_alpha", "rel_gamma_alpha_delta",
                "rel_gamma_delta_omicron",
                ## Lots of vaccination things
                "rel_severity_wildtype_alpha",
                "rel_severity_alpha_delta",
                "rel_severity_delta_omicron",
                "severity_cross_multiplier_alpha",
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
                "vacc_skip_to",
                "vacc_skip_progression_rate",
                "vacc_skip_weight",
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

  expected <- c("start_date", baseline$beta_names,
                ## severity
                "mu_D", "mu_D_2", "mu_D_3", "mu_D_4", "mu_D_5",
                "p_G_D", "p_G_D_2", "p_H", "p_H_2", "p_H_D", "p_ICU", "p_ICU_2",
                "p_ICU_D", "p_W_D",
                ## progression
                "mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4",
                # multistrain
                "ta_alpha", "seed_date_alpha",
                "ta_delta", "seed_date_delta",
                "ta_omicron", "seed_date_omicron",
                "rel_p_H_alpha", "rel_p_H_delta", "rel_p_H_omicron",
                "rel_p_ICU_alpha", "rel_p_ICU_delta", "rel_p_ICU_omicron",
                "rel_p_D_alpha", "rel_p_D_delta", "rel_p_D_omicron",
                ## observation
                "alpha_D", "alpha_H", "alpha_admission", "alpha_death_hosp",
                paste0("p_NC_", baseline$pillar2_age_bands),
                paste0("p_NC_weekend_", baseline$pillar2_age_bands),
                "rho_pillar2_tests")
  if (baseline$assumptions == "mu_d_summer") {
    expected <- c(expected, "mu_D_6")
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a

  function(pars) {
    stopifnot(setequal(expected, names(pars)))
    beta_value <- unname(pars[baseline$beta_names])
    pars <- as.list(pars) # using list access below
    
    if (is.null(date)) {
      date <- baseline$date
    }

    progression <- compute_progression(pars, baseline$progression_data)
    severity <- compute_severity(pars, baseline$severity_data,
                                 baseline$assumptions)
    observation <- compute_observation(pars, baseline$pillar2_age_bands,
                                       baseline$region)

    vaccine_schedule <- baseline$vaccine_schedule_effect
    
    stage_parameters <- function(strains, vaccine_doses) {
      
      if (strains == "Wildtype") {
        
        cross_immunity <- baseline$cross_immunity_wildtype
        strain_seed_date <- NULL
        strain_seed_size <- NULL
        strain_seed_pattern <- NULL
        strain_transmission <- 1
        rel_severity <- lapply(
          baseline$rel_severity_wildtype_alpha, function(x) x[, 1, , drop = FALSE])
        strain_rel_gamma <- baseline$rel_gamma_wildtype
        
        strain_rel_p_hosp_if_sympt <- 1
        strain_rel_p_icu <- 1
        strain_rel_p_death <- 1
        
      } else if (strains == "Wildtype_Alpha") {

        # Note cross immunity is c(against strain 2, against strain 1)
        cross_immunity <-
          c(baseline$cross_immunity_alpha, baseline$cross_immunity_wildtype)
        strain_seed_date <- pars$seed_date_alpha # new param
        strain_seed_size <- baseline$strain_seed_size
        strain_seed_pattern <- baseline$strain_seed_pattern
        strain_transmission <- c(1, pars$ta_alpha) # new param
        
        ## Assume the same rel_severity for Alpha as Wildtype 
        rel_severity <- baseline$rel_severity_wildtype_alpha
        strain_rel_gamma <- baseline$rel_gamma_wildtype_alpha
        
        ## Values for: Wildtype, Alpha, Wildtype -> Alpha, Alpha -> Wildtype
        strain_rel_p_hosp_if_sympt <-
          c(1, pars$rel_p_H_alpha,
            pars$rel_p_H_alpha * baseline$severity_cross_multiplier_alpha, 1)
        strain_rel_p_icu <-
          c(1, pars$rel_p_ICU_alpha, pars$rel_p_ICU_alpha, 1)
        strain_rel_p_death <-
          c(1, pars$rel_p_D_alpha, pars$rel_p_D_alpha, 1)
        
      } else if (strains == "Alpha_Delta") {
        
        # Note cross immunity is c(against strain 2, against strain 1)
        cross_immunity <-
          c(baseline$cross_immunity_delta, baseline$cross_immunity_alpha)
        strain_seed_date <- pars$seed_date_delta
        strain_seed_size <- baseline$strain_seed_size
        strain_seed_pattern <- baseline$strain_seed_pattern
        strain_transmission <- pars$ta_alpha * c(1, pars$ta_delta)
        rel_severity <- baseline$rel_severity_alpha_delta
        strain_rel_gamma <- baseline$rel_gamma_alpha_delta
        
        ## Values for: Alpha, Delta, Alpha -> Delta, Delta -> Alpha
        strain_rel_p_hosp_if_sympt <-
          pars$rel_p_H_alpha *
          c(1, pars$rel_p_H_delta, 
            pars$rel_p_H_delta * baseline$severity_cross_multiplier_delta,
            baseline$severity_cross_multiplier_alpha)
        strain_rel_p_icu <-
          pars$rel_p_ICU_alpha *
          c(1, pars$rel_p_ICU_delta, pars$rel_p_ICU_delta, 1)
        strain_rel_p_death <-
          pars$rel_p_D_alpha *
          c(1, pars$rel_p_D_delta, pars$rel_p_D_delta, 1)

      } else if (strains == "Delta_Omicron") {
        
        # Note cross immunity is c(against strain 2, against strain 1)
        cross_immunity <-
          c(baseline$cross_immunity_omicron, baseline$cross_immunity_delta)
        strain_seed_date <- pars$seed_date_omicron
        strain_seed_size <- baseline$strain_seed_size
        strain_seed_pattern <- baseline$strain_seed_pattern
        strain_transmission <-
          pars$ta_alpha * pars$ta_delta * c(1, pars$ta_omicron)
        rel_severity <- baseline$rel_severity_delta_omicron
        strain_rel_gamma <- baseline$rel_gamma_delta_omicron
        
        ## Values for: Delta, Omicron, Delta -> Omicron, Omicron -> Delta
        strain_rel_p_hosp_if_sympt <-
          pars$rel_p_H_alpha * pars$rel_p_H_delta *
          c(1, pars$rel_p_H_omicron, 
            pars$rel_p_H_omicron *
              baseline$severity_cross_multiplier_omicron$rel_p_hosp_if_sympt,
            baseline$severity_cross_multiplier_delta)
        strain_rel_p_icu <-
          pars$rel_p_ICU_alpha * pars$rel_p_ICU_delta *
          c(1, pars$rel_p_ICU_omicron, pars$rel_p_ICU_omicron, 1)
        strain_rel_p_death <-
          pars$rel_p_D_alpha * pars$rel_p_D_delta *
          c(1, pars$rel_p_D_omicron, 
            pars$rel_p_D_omicron *
              baseline$severity_cross_multiplier_omicron$rel_p_death,
            1)
        
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
        
        vacc_skip_progression_rate <- NULL
        vacc_skip_to <- NULL
        vacc_skip_weight <- NULL
        
      } else if (vaccine_doses == 2) {
        
        rel_severity <- lapply(
          rel_severity, function(x) x[, , 1:5, drop = FALSE])
        vaccine_progression_rate <- baseline$vaccine_progression_rate[1:5]
        vaccine_schedule$doses <- vaccine_schedule$doses[, 1:2, , drop = FALSE]
        vaccine_schedule$n_doses <- 2L
        vaccine_index_dose2 <- baseline$vaccine_index_dose2
        vaccine_index_booster <- NULL
        n_doses <- 2L
        
        vacc_skip_progression_rate <- NULL
        vacc_skip_to <- NULL
        vacc_skip_weight <- NULL
        
      } else if (vaccine_doses == 3) {
        
        vaccine_progression_rate <- baseline$vaccine_progression_rate[1:7]
        vaccine_schedule$doses <- vaccine_schedule$doses[, 1:3, , drop = FALSE]
        vaccine_schedule$n_doses <- 3L
        vaccine_index_dose2 <- baseline$vaccine_index_dose2
        vaccine_index_booster <- baseline$vaccine_index_booster
        n_doses <- 3L
        
        vacc_skip_progression_rate <- baseline$vacc_skip_progression_rate[1:7]
        vacc_skip_to <- baseline$vacc_skip_to
        vacc_skip_weight <- baseline$vacc_skip_weight
        
      } else {
        stop("vaccine_doses must be 0, 2 or 3")
      }
      
      sircovid::lancelot_parameters(
        start_date = pars$start_date,
        region = baseline$region,
        beta_date = baseline$beta_date,
        beta_value = beta_value,
        carehome_beds = 0,
        
        severity = severity,
        progression = progression,
        observation = observation,
        sens_and_spec = baseline$sens_and_spec,
        
        initial_seed_size = baseline$initial_seed_size,
        initial_seed_pattern = baseline$initial_seed_pattern,
        
        eps = 0,
        m_CHW = 0,
        m_CHR = 0,
        
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
        
        vacc_skip_progression_rate = vacc_skip_progression_rate,
        vacc_skip_to = vacc_skip_to,
        vacc_skip_weight = vacc_skip_weight,
        
        waning_rate = baseline$natural_waning_rate,
        cross_immunity = cross_immunity)
      
    }
    
    date <- sircovid::sircovid_date(date)
    p1 <- stage_parameters("Wildtype", 0)
    if (date < epoch_dates[1]) {
      ret <- p1
    } else {
      p2 <- stage_parameters("Wildtype_Alpha", 0)
      epochs <- list(mcstate::multistage_epoch(
        epoch_dates[1], p2, sircovid::inflate_state_strains))
      if (date >= epoch_dates[2]) {
        p3 <- stage_parameters("Wildtype_Alpha", 2)
        epochs <- append(epochs, list(mcstate::multistage_epoch(
          epoch_dates[2], p3, sircovid::inflate_state_vacc_classes)))
      }
      if (date >= epoch_dates[3]) {
        p4 <- stage_parameters("Alpha_Delta", 2)
        epochs <- append(epochs, list(mcstate::multistage_epoch(
          epoch_dates[3], p4, sircovid::rotate_strains)))
                    
      }
      if (date >= epoch_dates[4]) {
        p5 <- stage_parameters("Alpha_Delta", 3)
        epochs <- append(epochs, list(mcstate::multistage_epoch(
          epoch_dates[4], p5, sircovid::inflate_state_vacc_classes)))
      }
      if (date >= epoch_dates[5]) {
        p6 <- stage_parameters("Delta_Omicron", 3)
        epochs <- append(epochs, list(mcstate::multistage_epoch(
          epoch_dates[5], p6, sircovid::rotate_strains)))
      }
      ret <- mcstate::multistage_parameters(p1, epochs = epochs)
    }

    ret
  }
}
