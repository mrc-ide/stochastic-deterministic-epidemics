
## WARNING: not everything in here actually has an effect;
## vaccine_eligibility is a notable example.  Look in transform to see
## which things are actually used, and a further warning about this!

create_baseline <- function(region, date, restart_date,
                            epoch_dates, pars_info, assumptions) {
  message(sprintf("Creating baseline for '%s'", region))
  pars_info <- pars_info[pars_info$region == region | is.na(pars_info$region), ]
  pars_info <- setNames(pars_info$initial, pars_info$name)
  
  model_type <- "BB"
  
  ## 1. Set up basic variables ----
  date <- as.Date(date)
  dt <- 0.25
  vaccine_days_to_effect_dose1 <- 0
  vaccine_days_to_effect_dose2 <- 7
  vaccine_days_to_effect_booster <- 7
  vaccine_days_to_effect_booster2 <- 7
  vaccine_days_to_effect <- c(vaccine_days_to_effect_dose1,
                              vaccine_days_to_effect_dose2,
                              vaccine_days_to_effect_booster,
                              vaccine_days_to_effect_booster2)
  
  ## 2. Read in data ----
  
  # a. Local resources
  
  # Vaccines efficacy
  vaccine_efficacy_alpha <- read_csv("data/vaccine_efficacy_alpha.csv")
  vaccine_efficacy_delta <- read_csv("data/vaccine_efficacy_delta.csv")
  vaccine_efficacy_omicron <- read_csv("data/vaccine_efficacy_omicron.csv")

  # Vaccine uptake
  uptake_by_age <- read.csv("data/vaccine_uptake.csv", row.names = "group")
  
  # b. Dependencies
  data_vaccination <- read_csv("data_vaccination.csv")
  
  # c. External data (default parameters)
  # TODO: might make more sense to put these back into this task at some point
  progression_data <- read_csv(spimalot_file("extdata/support_progression.csv"))
  severity_data <- read_csv("data/support_severity.csv")
  
  
  ## 3. Set-up basic model parameters ----
  # Beta change points - A vector of date (strings) for the beta parameters.
  beta_date <-
    c("2020-03-16", ##  1. PM advises WFH, against non-essential travel etc
      "2020-03-23", ##  2. PM announces full lockdown
      "2020-03-25", ##  3. lockdown into full effect
      "2020-05-11", ##  4. initial easing of lockdown
      "2020-06-15", ##  5. non-essential shops can open
      "2020-07-04", ##  6. restaurants, pubs etc can open
      "2020-08-01", ##  7. "Eat out to help out" scheme starts
      "2020-09-01", ##  8. Schools reopen
      "2020-09-14", ##  9. "Rule of six" introduced
      "2020-10-14", ## 10. Tiered system introduced
      "2020-10-31", ## 11. lockdown announced
      "2020-11-05", ## 12. lockdown 2 starts
      "2020-12-02", ## 13. lockdown 2 ends
      "2020-12-18", ## 14. school Christmas holidays
      "2020-12-25", ## 15. last day of holidays season relaxation
      "2021-01-05", ## 16. Lockdown 3 starts
      "2021-03-08", ## 17. Step 1 of roadmap: schools reopen
      "2021-04-01", ## 18. Semi-arbitrary - school holidays / restart date
      "2021-04-19", ## 19. Step 2 of roadmap: outdoors hospitality (04-12)
      ##     and schools return (04-19)
      "2021-05-17", ## 20. Step 3 of roadmap: indoors hospitality
      "2021-06-21", ## 21. Step 3.5 - "freedom day" delayed/Euros last group match
      "2021-07-03", ## 22. Euros quarter final
      "2021-07-11", ## 23. Euros 2020 final - peak in transmission
      "2021-07-19", ## 24. Step 4
      "2021-08-15", ## 25. Summer festivals / holidays
      "2021-09-01", ## 26. Schools return
      "2021-09-22", ## 27. Mid-point between school start and half term
      ##    (help the model stabilise a long period of time)
      "2021-10-01", ## 28. Point before sharp increase in cases/hospitalisations
      "2021-10-22", ## 29. School half-term - Point before recent plateau in cases
      "2021-11-01", ## 30. Schools return
      "2021-12-08", ## 31. Announcement of move to Plan B
      "2021-12-23", ## 32. Starting of X-mas holidays
      "2022-01-04", ## 33. Schools return from X-mas holidays
      "2022-01-19", ## 34. Announcement of end of Plan B
      "2022-01-27", ## 35. End of Plan B
      "2022-02-24") ## 36. End date
  ## End of self-isolation in England
  
  ## Validate beta_date
  beta_date <- spimalot::spim_pars_check_beta_date(beta_date)
  beta_names <- sprintf("beta%d", seq_along(beta_date))
  
  pillar2_age_bands <- c("15_24", "25_49", "50_64", "65_79", "80_plus")
  
  # Set of parameters that will be fitted for each model type
  to_fit_all <- c(
    # observation
    "alpha_D", "alpha_H", "alpha_admission", "alpha_death_hosp",
    # direct
    "start_date", beta_names,
    # severity
    "mu_D", "mu_D_2", "mu_D_3", "mu_D_4", "mu_D_5",
    "p_G_D", "p_G_D_2", "p_H",
    "p_H_2", "p_H_D", "p_ICU", "p_ICU_2",
    "p_ICU_D", "p_W_D",
    # progression
    "mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4",
    # pillar 2 parameters
    paste0("p_NC_", pillar2_age_bands),
    paste0("p_NC_weekend_", pillar2_age_bands),
    "rho_pillar2_tests",
    # multistrain, direct
    "ta_alpha", "seed_date_alpha",
    "ta_delta", "seed_date_delta",
    "ta_omicron", "seed_date_omicron",
    "rel_p_H_alpha", "rel_p_ICU_alpha", "rel_p_D_alpha",
    "rel_p_H_delta", "rel_p_ICU_delta", "rel_p_D_delta",
    "rel_p_H_omicron", "rel_p_ICU_omicron", "rel_p_D_omicron"
    )
  if (assumptions == "mu_d_summer") {
    to_fit_all <- c(to_fit_all, "mu_D_6")
  }
  stopifnot(setequal(to_fit_all, names(pars_info)))
  
  # Waning of natural immunity rate (rate of exponential distribution)
  natural_waning_rate <- 1 / (3 * 365) # 3 years
  
  
  ## 4. Set-up fixed parameter assumptions and sensitivity analysis ----
  
  # a. Cross-immunity against VOCs
  
  # Note we assume perfect immunity against any given variant once recovered
  
  # Infection
  cross_immunity_wildtype <- 1
  cross_immunity_alpha <- 0.96 # 95% protection
  cross_immunity_delta <- 0.85 # 85% protection
  cross_immunity_omicron <- 0.25 # 25% protection
  
  # Hospitalisation by Alpha or Delta
  severity_cross_multiplier_alpha <- 1 - 0.85
  severity_cross_multiplier_delta <- 1 - 0.85
  
  # Hospitalisation & death by Omicron
  # Nyberg et al. Lancet 2022; 399: 1303–12, reports below values as HR
  # conditional on a positive test, unvaccinated status and documented prior
  # infection with historic variant. Respective SA use 95%CI from same paper.
  severity_cross_multiplier_omicron <- list(
    rel_p_hosp_if_sympt = 0.55,
    rel_p_death = 0.18 / 0.55
  )
  
  # Relative duration of SI and PCR positivity (compared to Wildtype)
  # We assume Alpha is 6%, Delta 13% and Omicron 25% shorter than Wildtype
  rel_si_wildtype <- 1
  rel_si_alpha <- 0.94
  rel_si_delta <- 0.87
  rel_si_omicron <- 0.75
  
  
  # Note that cross_immunity against infection parameters are actual cross-
  # immunity, whereas death and hospitalisation are relative severity parameters.
  if (assumptions == "crim_infect_high") {
    
    cross_immunity_omicron <- 0.35
  
  } else if (assumptions == "crim_infect_low") {
    
    cross_immunity_omicron <- 0.20
  
  } else if (assumptions == "crim_death_low") {
    
    severity_cross_multiplier_omicron$rel_p_death <- 0.57 / 0.63
    
  } else if (assumptions == "crim_death_high") {
    
    severity_cross_multiplier_omicron$rel_p_death <- 0.06 / 0.48
    
  } else if (assumptions == "crim_hospi_low") {
    
    severity_cross_multiplier_omicron$rel_p_hosp_if_sympt <- 0.63
    
  } else if (assumptions == "crim_hospi_high") {
    
    severity_cross_multiplier_omicron$rel_p_hosp_if_sympt <- 0.48
    
  } else if (assumptions == "fixed_si_high") {
    
    rel_si_alpha <- 1
    rel_si_delta <- 1
    rel_si_omicron <- 1
    
  } else if (assumptions == "fixed_si_low") {
    
    rel_si_wildtype <- 0.75
    rel_si_alpha <- 0.75
    rel_si_delta <- 0.75
    rel_si_omicron <- 0.75
    
  } else {
    
    ve_assumption_names <- unique(c(colnames(vaccine_efficacy_alpha), 
                                    colnames(vaccine_efficacy_delta),
                                    colnames(vaccine_efficacy_omicron)))
    other_assumptions <- c("central", "booster_ve_high", "booster_ve_low",
                           "mu_d_summer", "mu_d_winter", "alpha_ve_low",
                           "alpha_ve_high", "delta_ve_low", "delta_ve_high")
    
    stopifnot(assumptions %in% c(ve_assumption_names, other_assumptions))
    
  }
  
  strain_seed_size <- 7 * 2e-6 * sum(sircovid:::sircovid_population(region))
  strain_seed_pattern <- rep(1, 7 * 4)
  
  
  # Serial interval of variants
  # Central mean SI for Wildtype is 5.2 days, based on calculation
  # adapted from Svensson (https://www2.math.su.se/matstat/reports/seriea/2005/rep14/report.pdf)
  # so that SI = E + (P^2 + P·C_1 + C_1^2) / (P + C_1)
  ret_si <- data.frame(
    parameter = c("k_E", "gamma_E"),
    value = c(2, 0.865))
  progression_data <- rbind(progression_data, ret_si)
  
  # We ensure that mean[T_I_C_1 + T_I_C_2] is unchanged
  mean_E <- 2 / 0.865
  mean_P <- 1.68
  mean_C_1 <- 2.14
  mean_C_2 <- 1.86
  
  rel_C_2_wildtype <- (mean_C_2 + (1 - rel_si_wildtype) * mean_C_1) / mean_C_2
  rel_C_2_alpha <- (mean_C_2 + (1 - rel_si_alpha) * mean_C_1) / mean_C_2
  rel_C_2_delta <- (mean_C_2 + (1 - rel_si_delta) * mean_C_1) / mean_C_2
  rel_C_2_omicron <- (mean_C_2 + (1 - rel_si_omicron) * mean_C_1) / mean_C_2
  
  rel_gamma_wildtype <- list(E = 1 / rel_si_wildtype,
                             A = 1 / rel_si_wildtype,
                             P = 1 / rel_si_wildtype,
                             C_1 = 1 / rel_si_wildtype,
                             C_2 = 1 / rel_C_2_wildtype)
  
  rel_gamma_wildtype_alpha <- list(E = c(1 / rel_si_wildtype, 1 / rel_si_alpha),
                                   A = c(1 / rel_si_wildtype, 1 / rel_si_alpha),
                                   P = c(1 / rel_si_wildtype, 1 / rel_si_alpha),
                                   C_1 = c(1 / rel_si_wildtype, 1 / rel_si_alpha),
                                   C_2 = c(1 / rel_C_2_wildtype, 1 / rel_C_2_alpha))
  
  rel_gamma_alpha_delta <- list(E = c(1 / rel_si_alpha, 1 / rel_si_delta),
                                A = c(1 / rel_si_alpha, 1 / rel_si_delta),
                                P = c(1 / rel_si_alpha, 1 / rel_si_delta),
                                C_1 = c(1 / rel_si_alpha, 1 / rel_si_delta),
                                C_2 = c(1 / rel_C_2_alpha, 1 / rel_C_2_delta))
  
  rel_gamma_delta_omicron <- list(E = c(1 / rel_si_delta, 1 / rel_si_omicron),
                                  A = c(1 / rel_si_delta, 1 / rel_si_omicron),
                                  P = c(1 / rel_si_delta, 1 / rel_si_omicron),
                                  C_1 = c(1 / rel_si_delta, 1 / rel_si_omicron),
                                  C_2 = c(1 / rel_C_2_delta, 1 / rel_C_2_omicron))
  
  
  ## 5. Set-up vaccination parameters and assumptions ----
  
  # Vaccine eligibility by age; min_18 only here, though we never
  # actually use it in the fits which are driven by data.
  ## TODO: remove from here; best placed in parameters simulation task/will remove form here soon.  
  vaccine_eligibility_min_age <- 5
  
  # Average duration of stay in each vaccinated compartment
  mean_days_between_doses <- 7 * 11 # second dose starts 12 weeks after first
  
  stopifnot(length(unique(vaccine_efficacy_delta$week_wane)) == 1,
            length(unique(vaccine_efficacy_omicron$week_wane)) == 1,
            vaccine_efficacy_delta$week_wane[1] == vaccine_efficacy_omicron$week_wane[1])
  
  mean_time_to_waned <- (vaccine_efficacy_delta$week_wane[1]) * 7
  
  time_to_dose_1_effect <- 7 * 3
  
  vaccine_efficacy_delta$week_wane <- NULL
  vaccine_efficacy_omicron$week_wane <- NULL
  
  # Compartment progression rates
  # These are progression rates OUT of the specified compartment, where
  # zeroes mean movement is controlled by vaccination schedule rather than
  # a rate parameter and booster_waned is an absorbing state
  vaccine_progression_rate <- c(0,                         # unvaccinated 
                                1 / time_to_dose_1_effect, # first dose no effect
                                0,                         # first dose full effect
                                1/ mean_time_to_waned,     # second dose
                                0,                         # waned
                                1 / mean_time_to_waned,    # booster
                                0)                         # booster_waned
  
  ## parameters controlling the "vaccine skip moves":
  ## there is one vaccine skip move from vaccine strata 4 (two doses, not waned)
  ## to 6 (boosted)
  ## base progression rate is 0 (progression will be controlled by doses anyway)
  ## weighting of non-waned relative to waned for doses is 1 (i.e. equal weight)
  vacc_skip_to <- c(0L, 0L, 0L, 6L, 0L, 0L, 0L)
  vacc_skip_progression_rate <- rep(0, 7)
  vacc_skip_weight <- c(0, 0, 0, 1, 0, 0, 0)
  
  # Proportion eligible for booster - all over 50s and CEV proportion from ONS
  vaccine_booster_proportion <- c(0, 0, 0, 0.004, # under 20s
                                  0.004, 0.004, # 20 - 29
                                  0.005, 0.005, # 30 - 39
                                  0.008, 0.008, # 40 - 49
                                  rep(1L, 9)) # 50 - 80+ general pop, CHW, CHR)
  
  
  ## 6. Set-up vaccine efficacy from report table ----
  
  # Data on vaccine doses by age, date and type of vaccine
  vacc_doses_by_age_date_vaccine <-
    data_vaccination %>%
    dplyr::mutate(age_band_min = replace(age_band_min, age_band_min == 16, 15)) %>%
    dplyr::filter(!is.na(age_band_min)) %>%
    dplyr::group_by(vaccine, age_band_min) %>%
    dplyr::summarise(n = sum(first_dose, na.rm = TRUE)) %>%
    dplyr::group_by(age_band_min) %>%
    dplyr::mutate(freq = n / sum(n, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = c(age_band_min),
                       values_from = freq, names_from = vaccine)
  
  ## check all proportions sum to 1
  stopifnot(
    all(abs(rowSums(
      vacc_doses_by_age_date_vaccine[3:17, -1], na.rm = TRUE) - 1) < 1e-8))
  
  prop_pfizer <- apply(cbind(vacc_doses_by_age_date_vaccine$Pfizer,
    vacc_doses_by_age_date_vaccine$Moderna), 1, sum, na.rm = TRUE)
  
  ## set %PF in CHW / CHR <-  80+
  # prop_pfizer <- c(prop_pfizer, rep(prop_pfizer[17], 2))
  # 100% PF in 0-10
  prop_pfizer[sircovid:::sircovid_age_bins()$end < 10] <- 1
  
  
  ## Add assumed booster_waned VE for Alpha and Delta - we have values for Omicron
  vaccine_efficacy_alpha <- 
    add_vaccine_efficacy_booster_waned(vaccine_efficacy_alpha)
  vaccine_efficacy_delta <- 
    add_vaccine_efficacy_booster_waned(vaccine_efficacy_delta)
  
  
  ## Calculate average vaccine efficacy and rel_severity
  average_vacc_efficacy_alpha <-
    calculate_average_vacc_efficacy(vaccine_efficacy_alpha, prop_pfizer)
  rel_severity_alpha <- lapply(average_vacc_efficacy_alpha, function(e)
    get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                 e$infection, e$transmission))
  
  average_vacc_efficacy_delta <-
    calculate_average_vacc_efficacy(vaccine_efficacy_delta, prop_pfizer)
  rel_severity_delta <- lapply(
    average_vacc_efficacy_delta, function(e)
      get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                   e$infection, e$transmission))
  
  average_vacc_efficacy_omicron <-
    calculate_average_vacc_efficacy(vaccine_efficacy_omicron, prop_pfizer)
  rel_severity_omicron <- lapply(
    average_vacc_efficacy_omicron, function(e)
      get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                   e$infection, e$transmission))
  
  
  ## Calculate rel_severity across strain transitions - one-way re_severity
  ## sensitivity analyses only, so all other strains keep central values
  rel_severity_wildtype_alpha <-
    lapply(names(rel_severity_alpha),
           function (x)
             rel_severity_strains(rel_severity_alpha[["central"]],
                                  rel_severity_alpha[[x]]))
  names(rel_severity_wildtype_alpha) <- names(rel_severity_alpha)
  
  rel_severity_alpha_delta <-
    lapply(names(rel_severity_delta),
           function(x)
             rel_severity_strains(rel_severity_alpha[["central"]],
                                  rel_severity_delta[[x]]))
  names(rel_severity_alpha_delta) <- names(rel_severity_delta)
  
  rel_severity_delta_omicron <- 
    lapply(names(rel_severity_omicron),
           function(x) 
             rel_severity_strains(rel_severity_delta[["central"]],
                                  rel_severity_omicron[[x]]))
  names(rel_severity_delta_omicron) <- names(rel_severity_omicron)
  
  
  ## Set right VE assumption values across sensitivity analyses
  if (assumptions == "booster_ve_high") {
    
    rel_severity_delta_omicron <- rel_severity_delta_omicron$booster_ve_high
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$central
    rel_severity_alpha_delta <- rel_severity_alpha_delta$central
    
  } else if (assumptions == "booster_ve_low") {
    
    rel_severity_delta_omicron <- rel_severity_delta_omicron$booster_ve_low
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$central
    rel_severity_alpha_delta <- rel_severity_alpha_delta$central
    
  } else if (assumptions == "alpha_ve_high") {
    
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$alpha_ve_high
    rel_severity_alpha_delta <- 
      rel_severity_strains(rel_severity_alpha[[assumptions]],
                           rel_severity_delta[["central"]])
    rel_severity_delta_omicron <- rel_severity_delta_omicron$central
    
  } else if (assumptions == "alpha_ve_low") {
    
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$alpha_ve_low
    rel_severity_alpha_delta <- 
      rel_severity_strains(rel_severity_alpha[[assumptions]],
                           rel_severity_delta[["central"]])
    rel_severity_delta_omicron <- rel_severity_delta_omicron$central
    
  } else if (assumptions == "delta_ve_low") {
    
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$central
    rel_severity_alpha_delta <- rel_severity_alpha_delta$delta_ve_low
    rel_severity_delta_omicron <-
      rel_severity_strains(rel_severity_delta[[assumptions]],
                           rel_severity_omicron[["central"]])
    
  } else if (assumptions == "delta_ve_high") {
    
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$central
    rel_severity_alpha_delta <- rel_severity_alpha_delta$delta_ve_high
    rel_severity_delta_omicron <-
      rel_severity_strains(rel_severity_delta[[assumptions]],
                           rel_severity_omicron[["central"]])
    
  } else {
    
    rel_severity_wildtype_alpha <- rel_severity_wildtype_alpha$central
    rel_severity_delta_omicron <- rel_severity_delta_omicron$central
    rel_severity_alpha_delta <- rel_severity_alpha_delta$central

  }
  
  sens_and_spec <-
    sircovid::lancelot_parameters_sens_and_spec(sero_sensitivity_1 = 1,
                                                sero_specificity_1 = 0.99,
                                                sero_sensitivity_2 = 1,
                                                sero_specificity_2 = 0.99,
                                                pillar2_sensitivity = 1,
                                                pillar2_specificity = 1,
                                                ons_sensitivity = 1,
                                                ons_specificity = 1,
                                                react_sensitivity = 1,
                                                react_specificity = 1)
  
  ## Note that vaccine_uptake[i, j] is proportional uptake of dose j for group i 
  vaccine_uptake <- 
    array(uptake_by_age$central, c(length(uptake_by_age$central), 2))
  vaccine_uptake <- cbind(vaccine_uptake,
                          vaccine_uptake[, 1] * vaccine_booster_proportion)
  
  data_vaccination <- data_vaccination %>%
    dplyr::rename(dose1 = first_dose,
                  dose2 = second_dose,
                  dose3 = third_dose,
                  dose4 = fourth_dose,
                  dose5 = fifth_dose)
  n_doses <- 3
  dose_start_dates <- c("2020-12-08",
                        "2020-12-08",
                        "2021-09-15")
  vaccination <- 
    spimalot::spim_vaccination_data(as.Date("2022-02-24"), region,
                                    vaccine_uptake, vaccine_days_to_effect,
                                    data_vaccination, n_doses, dose_start_dates,
                                    carehomes = FALSE)
  vaccination$schedule$doses <- 
    vaccination$schedule$doses[1:17, , , drop = FALSE]
  
  vaccine_schedule_real <- vaccination$schedule
  
  ## shift doses to account for time between dose and effect of dose
  vaccine_schedule_effect <- shift_doses(vaccine_schedule_real,
                                         vaccine_days_to_effect)
  
  n_doses <- vaccination$schedule$n_doses
  vaccine_index_booster <- 5L
  vaccine_index_dose2 <- 3L
  
  ## Initial seeding: seed 10 per million over 1 day (4 steps)
  initial_seed_size <- 10e-6 * sum(sircovid:::sircovid_population(region))
  initial_seed_pattern <- c(1, 1, 1, 1)
  
  
  ## Dates for which we will calculate intrinsic severity
  ## (i.e. excluding immunity)
  intrinsic_severity_dates <- c(
    # Have a range for evaluating severity at Emergence
    "Emergence1" = sircovid_date(as.Date("2020-03-16")),
    "Emergence2" = sircovid_date(as.Date("2020-04-16")),
    "Emergence3" = sircovid_date(as.Date("2020-05-16")),
    "Emergence4" = sircovid_date(as.Date("2020-06-16")),
    # Monthly intervals during the winter 2020/21 peak
    "Sep 2020" = sircovid_date(as.Date("2020-09-01")),
    "Oct 2020" = sircovid_date(as.Date("2020-10-01")),
    "Nov 2020" = sircovid_date(as.Date("2020-11-01")),
    "Dec 2020" = sircovid_date(as.Date("2020-12-01")),
    "Jan 2021" = sircovid_date(as.Date("2021-01-01")),
    "Feb 2021" = sircovid_date(as.Date("2021-02-01")),
    "Mar 2021" = sircovid_date(as.Date("2021-03-01")),
    # Omicron era
    "Omicron1" = sircovid_date(as.Date("2021-11-01")),
    "Omicron2" = sircovid_date(as.Date("2021-12-01")),
    "Omicron3" = sircovid_date(as.Date("2022-01-01")),
    "Omicron4" = sircovid_date(as.Date("2022-02-01"))
  )
  ## The epochs where each strain is introduced
  strain_epochs <- c(
    "Wildtype" = 0,
    "Alpha" = 1,
    "Delta" = 3,
    "Omicron" = 5
  )
  
  ## TODO: We save epoch date as a sircovid date, but restart date as
  ## a string; we should pick one, but this does have some
  ## implications for onward functions.

  ## TODO: Some of these need to be present for the parameter construction
  ## to work but we will replace later.
  ret <- list(
    date = date,
    model_type = model_type,
    region = region,
    assumptions = assumptions,
    restart_date = restart_date,
    epoch_dates = sircovid::sircovid_date(epoch_dates),
    intrinsic_severity_dates = intrinsic_severity_dates,
    strain_epochs = strain_epochs,
    pillar2_age_bands = pillar2_age_bands,
    
    beta_date = sircovid::sircovid_date(beta_date),
    beta_names = beta_names,
    severity_data = severity_data,
    progression_data = progression_data,
    sens_and_spec = sens_and_spec,
    initial_seed_size = initial_seed_size,
    initial_seed_pattern = initial_seed_pattern,
    natural_waning_rate = natural_waning_rate,
    
    cross_immunity_wildtype = cross_immunity_wildtype,
    cross_immunity_alpha = cross_immunity_alpha,
    cross_immunity_delta = cross_immunity_delta,
    cross_immunity_omicron = cross_immunity_omicron,
    rel_severity_wildtype_alpha = rel_severity_wildtype_alpha,
    rel_severity_alpha_delta = rel_severity_alpha_delta,
    rel_severity_delta_omicron = rel_severity_delta_omicron,
    severity_cross_multiplier_alpha = severity_cross_multiplier_alpha,
    severity_cross_multiplier_delta = severity_cross_multiplier_delta,
    severity_cross_multiplier_omicron = severity_cross_multiplier_omicron,
    rel_gamma_wildtype = rel_gamma_wildtype,
    rel_gamma_wildtype_alpha = rel_gamma_wildtype_alpha,
    rel_gamma_alpha_delta = rel_gamma_alpha_delta,
    rel_gamma_delta_omicron = rel_gamma_delta_omicron,
    
    vaccine_eligibility_min_age = vaccine_eligibility_min_age,
    vaccine_progression_rate = vaccine_progression_rate,
    ## TODO: change this to vaccine_schedule_real when safe to do so
    vaccine_schedule = vaccine_schedule_real,
    vaccine_schedule_effect = vaccine_schedule_effect,
    vaccine_uptake = vaccine_uptake,
    vaccine_mean_days_between_doses = mean_days_between_doses,
    vaccine_index_dose2 = vaccine_index_dose2,
    vaccine_index_booster = vaccine_index_booster,
    vaccine_days_to_effect = vaccine_days_to_effect,
    vacc_skip_to = vacc_skip_to,
    vacc_skip_weight = vacc_skip_weight,
    vacc_skip_progression_rate = vacc_skip_progression_rate,
    n_doses = n_doses,
    strain_seed_size = strain_seed_size,
    strain_seed_pattern = strain_seed_pattern)
  
  message("  - Creating transformation function")
  tr <- make_transform(apply_assumptions(ret, assumptions))
  message("  - Testing transformation function")
  p <- tr(pars_info)
  message("  - Testing creating model with transformed parameters")
  for (i in seq_along(p)) {
    m <- sircovid::lancelot$new(p[[i]]$pars, 0, 1)
  }
  
  ret
}


rel_severity_strains <- function(vacc_rel_severity_strain1,
                                 vacc_rel_severity_strain2) {
  strain_severity_modifier <- rep(list(list(
    rel_susceptibility = 1,
    rel_p_sympt = 1,
    rel_p_hosp_if_sympt = 1,
    rel_infectivity = 1,
    rel_p_death = 1
  )), 4)
  
  ## modify_severity modifies severity to make it relative to alpha/wildtype according to;
  ## efficacy, efficacy of strain 2 and the strain severity modifier
  ## see `sircovid\R\vaccination.R` for more details
  sircovid::modify_severity(vacc_rel_severity_strain1,
                            vacc_rel_severity_strain2,
                            strain_severity_modifier)
}


supl_fig_vac_age <- function() {
  
  data_vaccination <- read_csv("data_vaccination.csv")
  # Figure of proportion of each vaccine by age
  fig_sup_vacc_age <- data_vaccination %>%
    dplyr::mutate(age_band_min = replace(age_band_min, age_band_min == 16, 15)) %>%
    dplyr::filter(!is.na(age_band_min)) %>%
    dplyr::group_by(vaccine, age_band_min) %>%
    dplyr::summarise(n = sum(first_dose, na.rm = TRUE)) %>%
    dplyr::group_by(age_band_min) %>%
    dplyr::mutate(freq = n / sum(n, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = c(age_band_min),
                       values_from = freq, names_from = vaccine) %>%
    mutate(
      age_band = factor(case_when(
        age_band_min != 80 ~ paste0(age_band_min, " to ", age_band_min + 5),
        TRUE ~ "80 plus")),
      Mod = Moderna,
      PF = Pfizer) %>%
    ungroup() %>%
    select(!c(age_band_min, Janssen, Novavax, Moderna, Pfizer, `NA`)) %>%
    pivot_longer(!age_band, names_to = "Vaccine")
  

  fig_sup_vacc_age %>%
    ggplot2::ggplot(., ggplot2::aes(age_band, value, fill = Vaccine)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "", y = "Proportion") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

}
