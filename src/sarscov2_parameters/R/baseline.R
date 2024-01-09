
## WARNING: not everything in here actually has an effect;
## vaccine_eligibility is a notable example.  Look in transform to see
## which things are actually used, and a further warning about this!

create_baseline <- function(region, date, model_type,
                            restart_date, epoch_dates, pars_info, assumptions) {
  message(sprintf("Creating baseline for '%s'", region))
  pars_info <- pars_info[pars_info$region == region | is.na(pars_info$region), ]
  pars_info <- setNames(pars_info$initial, pars_info$name)

  spimalot::spim_check_model_type(model_type)

  ## 1. Set up basic variables ----
  date <- as.Date(date)
  dt <- 0.25
  vaccine_days_to_effect_dose1 <- 21
  vaccine_days_to_effect_dose2 <- 7
  vaccine_days_to_effect_booster <- 7
  vaccine_days_to_effect <- c(vaccine_days_to_effect_dose1,
                              vaccine_days_to_effect_dose2,
                              vaccine_days_to_effect_booster)

  ## 2. Read in data ----

  # a. Local resources
  
  # Vaccines efficacy
  vaccine_efficacy_alpha <- read_csv("vaccine_data/vaccine_efficacy_alpha.csv")
  vaccine_efficacy_delta <- read_csv("vaccine_data/vaccine_efficacy_delta.csv")
  vaccine_efficacy_omicron <- read_csv("vaccine_data/vaccine_efficacy_omicron.csv")

  # Vaccine uptake
  uptake_by_age <- read.csv("vaccine_data/vaccine_uptake.csv", row.names = "group")

  # b. Dependencies
  data_vaccination <- read_csv("data_vaccination.csv")
  weighted_severity_priors <- read_csv("weighted_prior_ranges.csv")

  # c. External data (default parameters)
  # TODO: might make more sense to put these back into this task at some point
  progression_data <- read_csv(spimalot_file("extdata/support_progression.csv"))
  severity_data <- read_csv(spimalot_file("extdata/support_severity.csv"))

  # Time in exposed compartment will be halved to 1.71 days

  # we want to change gamma_E_step to set new discrete latency 1.1 day
  # shorter than before. Current mean GT is 6.3 and new one should be 5.2)
  ret_si <- data.frame(
    parameter = c("k_E", "gamma_E"),
    value = c(2, 0.865))
  progression_data <- rbind(progression_data, ret_si)
  
  ## Assume Alpha and Delta have same SI
  rel_gamma_alpha <- list(E = 1,
                          A = 1,
                          P = 1,
                          C_1 = 1,
                          C_2 = 1)
  rel_gamma_alpha_delta <- list(E = 1,
                                A = 1,
                                P = 1,
                                C_1 = 1,
                                C_2 = 1)
  
  ## We assume Omicron SI is 35% shorter than Delta by multiplying
  ## the relevant gammas by a factor of 1 / 0.65. Then we ensure that
  ## mean[T_I_C_1 + T_I_C_2] is unchanged
  rel_si <- 0.65
  mean_C_1 <- 2.14
  mean_C_2 <- 1.86
  rel_C_2 <- (mean_C_2 + (1 - rel_si) * mean_C_1) / mean_C_2 
  rel_gamma_delta_omicron <- list(E = c(1, 1 / rel_si),
                                  A = c(1, 1 / rel_si),
                                  P = c(1, 1 / rel_si),
                                  C_1 = c(1, 1 / rel_si),
                                  C_2 = c(1, 1 / rel_C_2))
    
  
  ## 3. Set-up basic model parameters and assumptions ----
  # Beta change points
  ## beta_date - A vector of date (strings) for the beta parameters.
  ## parent_end_date - end date of the parent fit (single strain)
  last_beta_days_ago <- 7
  ## Dates are as follows
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
      "2021-09-13") ## 27. End of fits
  ## Validate beta_date
  beta_date <- spimalot::spim_pars_check_beta_date(beta_date)
  beta_names <- sprintf("beta%d", seq_along(beta_date))

  # Set of parameters that will be fitted for each model type
  to_fit_common <- c(
    # observation
    "alpha_D", "alpha_H", "alpha_death_hosp",
    # direct
    "eps", "m_CHR", "m_CHW", "start_date", beta_names,
    # severity
    "mu_D", "mu_D_2", "p_G_D", "p_G_D_CHR", "p_H",
    "p_H_2", "p_H_3", "p_H_D", "p_ICU", "p_ICU_2",
    "p_ICU_D", "p_W_D",
    # progression
    "mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4",
    # multistrain, direct
    "ta_delta", "seed_date_delta",
    "rel_p_H_delta", "rel_p_ICU_delta", "rel_p_D_delta")

  pillar2_age_bands <- c("15_24", "25_49", "50_64", "65_79", "80_plus")
  
  to_fit_bb <-  c("p_NC", "p_NC_weekend", "rho_pillar2_tests")
  
  to_fit_nb <- c("alpha_pillar2_cases", "phi_pillar2_cases",
                 "phi_pillar2_cases_weekend",
                 paste0("phi_pillar2_cases_", pillar2_age_bands),
                 paste0("phi_pillar2_cases_weekend_", pillar2_age_bands))

  to_fit_all <- c(to_fit_common,
                  if (model_type == "BB") to_fit_bb else to_fit_nb)

  ## TODO: better messages here, particularly the last one. This would
  ## indicate that the parameters we think we're fitting is out of sync
  ## with those that we have information on, and that could come either
  ## way.
  stopifnot(!any(to_fit_bb %in% to_fit_common),
            !any(to_fit_nb %in% to_fit_common),
            !any(to_fit_bb %in% to_fit_nb),
            setequal(to_fit_all, names(pars_info)))

  # Waning of natural immunity rate (rate of exponential distribution)
  natural_waning_rate <- list(
    optimistic = 0, # no waning
    central = 1 / (6 * 365), # 6 years
    pessimistic = 1 / (3 * 365) # 3 years
  )

  # Cross-immunity between the two strains
  cross_immunity_wildtype_alpha <- c(1, 1)
  cross_immunity_alpha_delta <- c(0.85, 1)
  
  cross_immunity_delta_omicron <- list(
    central = c(0.75, 1))

  strain_seed_size <- 7 * 2e-6 * sum(sircovid:::sircovid_population(region))
  strain_seed_pattern <- rep(1, 7 * 4)

  ## 4. Set-up vaccination parameters and assumptions ----

  # Vaccine eligibility by age; min_18 only here, though we never
  # actually use it in the fits which are driven by data.
  ## TODO: remove from here; best placed in parameters simulation task/will remove form here soon.  
  vaccine_eligibility_min_age <- 5

  # Average duration of stay in each vaccinated compartment
  mean_days_between_doses <- 7 * 11 # second dose starts 12 weeks after first
  # mean_time_to_full_effect <- c(vaccine_days_to_effect_dose1,
  #                               vaccine_days_to_effect_dose2)
  
  stopifnot(length(unique(vaccine_efficacy_delta$week_wane)) == 1,
            length(unique(vaccine_efficacy_omicron$week_wane)) == 1,
            vaccine_efficacy_delta$week_wane[1] == vaccine_efficacy_omicron$week_wane[1])
  
  mean_time_to_waned <- (vaccine_efficacy_delta$week_wane[1]) * 7
  
  vaccine_efficacy_delta$week_wane <- NULL
  vaccine_efficacy_omicron$week_wane <- NULL
  
  # TODO: specify this parameter once booster class is added
  # mean_time_to_boosting <- 365 / 2
  # mean_duration <- c(mean_time_to_full_effect[1],
  #                    mean_days_between_doses -
  #                    mean_time_to_full_effect[1] +
  #                    mean_time_to_full_effect[2],
  #                    mean_time_to_waned)

  # Compartment progression rates
  # unvaccinated / first dose / second dose / waned / booster
  # TODO: add a fifth compartment for boosted class
  # vaccine_progression_rate <- c(0, 1 / mean_duration, 0)
  vaccine_progression_rate <- c(0, 0, 1/ mean_time_to_waned, 0, 0)
  
  # Proportion eligible for booster - all over 50s and CEV proportion from ONS
  vaccine_booster_proportion <- c(0, 0, 0, 0.004, # under 20s
                                  0.004, 0.004, # 20 - 29
                                  0.005, 0.005, # 30 - 39
                                  0.008, 0.008, # 40 - 49
                                  rep(1L, 9)) # 50 - 80+ general pop, CHW, CHR)

  ## 5. Set-up vaccine efficacy from report table ----

  # Data on vaccine doses by age, date and type of vaccine
  vacc_doses_by_age_date_vaccine <-
    data_vaccination %>%
    dplyr::mutate(age_band_min = replace(age_band_min, age_band_min == 16, 15)) %>%
    dplyr::filter(!is.na(age_band_min)) %>%
    dplyr::group_by(vaccine, age_band_min) %>%
    dplyr::summarise(n = sum(first_dose, na.rm = TRUE)) %>%
    dplyr::group_by(age_band_min) %>%
    dplyr::mutate(freq = n / sum(n, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = c(age_band_min), values_from = freq, names_from = vaccine)

  ## check all proportions sum to 1
  stopifnot(
    all(
      abs(rowSums(vacc_doses_by_age_date_vaccine[3:17, -1], na.rm = TRUE) - 1) < 1e-8
    )
  )

  prop_pfizer <- vacc_doses_by_age_date_vaccine$Pfizer + vacc_doses_by_age_date_vaccine$Moderna

  ## set %PF in CHW / CHR <-  80+
  prop_pfizer <- c(prop_pfizer, rep(prop_pfizer[17], 2))
  
  # 100% PF in 0-10
  prop_pfizer[sircovid:::sircovid_age_bins()$end < 10] <- 1

  ## Note that we do not change efficacy to allow for differences in future uptake
  average_vacc_efficacy_alpha <-
    calculate_average_vacc_efficacy(vaccine_efficacy_alpha, prop_pfizer)

  ## format that is readable by sircovid
  ##
  ## TODO: Simplify the csv file too and the loading here as we will
  ## only use a central estimate now (these apply only to alpha)
  rel_severity_alpha <- lapply(average_vacc_efficacy_alpha, function(e)
    get_vaccine_conditional_prob(e$death,
                                 e$severe_disease,
                                 e$disease,
                                 e$infection,
                                 e$transmission))$central

  
    ## This all can be refactored as it's not actually doing anything
  ## nearly as complicated here as we are doing as this was
  ## previously much more general.
  average_vacc_efficacy_delta <-
    calculate_average_vacc_efficacy(vaccine_efficacy_delta, prop_pfizer)
  rel_severity_delta <- lapply(
    average_vacc_efficacy_delta, function(e)
      get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                   e$infection, e$transmission))$central

  rel_severity_alpha_delta <-
    rel_severity_strains(rel_severity_alpha,
                         rel_severity_delta)
  
  rel_severity_wildtype_alpha <-
    rel_severity_strains(rel_severity_alpha, rel_severity_alpha)
  
  ## make sure infection with Alpha confers some protection against
  ## hospitalisation with Delta similar to 2 doses pfizer 
  severity_cross_multiplier_delta <- 
    (1 - vaccine_efficacy_delta$central[
      vaccine_efficacy_delta$vaccine == "PF" &
        vaccine_efficacy_delta$type == "severe_disease" &
        vaccine_efficacy_delta$dose == 2]) /
    (1 - vaccine_efficacy_delta$central[
      vaccine_efficacy_delta$vaccine == "PF" &
        vaccine_efficacy_delta$type == "disease" &
        vaccine_efficacy_delta$dose == 2])
  
                         
  average_vacc_efficacy_omicron <-
    calculate_average_vacc_efficacy(vaccine_efficacy_omicron, prop_pfizer)
  rel_severity_omicron <- lapply(
    average_vacc_efficacy_omicron, function(e)
      get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                   e$infection, e$transmission))
  
  rel_severity_delta_omicron <- 
    lapply(names(rel_severity_omicron),
           function(x) 
             rel_severity_strains(rel_severity_delta,
                                  rel_severity_omicron[[x]]))
  names(rel_severity_delta_omicron) <- names(rel_severity_omicron)
  
  ## make sure infection with Delta confers some protection against
  ## hospitalisation with Omicron similar to 2 doses pfizer
  severity_cross_multiplier_omicron <-
    lapply(names(rel_severity_omicron),
           function(nm) 
             (1 - vaccine_efficacy_omicron[[nm]][
               vaccine_efficacy_omicron$vaccine == "PF" &
                 vaccine_efficacy_omicron$type == "severe_disease" &
                 vaccine_efficacy_omicron$dose == 2]) /
             (1 - vaccine_efficacy_omicron[[nm]][
               vaccine_efficacy_omicron$vaccine == "PF" &
                 vaccine_efficacy_omicron$type == "disease" &
                 vaccine_efficacy_omicron$dose == 2]))
  names(severity_cross_multiplier_omicron) <- 
    names(rel_severity_omicron)
  
  sens_and_spec <-
    sircovid::lancelot_parameters_sens_and_spec(sero_sensitivity_1 = 1,
                                                sero_specificity_1 = 0.99,
                                                sero_sensitivity_2 = 1,
                                                sero_specificity_2 = 0.99,
                                                pillar2_sensitivity = 1,
                                                pillar2_specificity = 1,
                                                react_sensitivity = 1,
                                                react_specificity = 1)
  
  ## Note that vaccine_uptake[i, j] is proportional uptake of dose j for group i 
  vaccine_uptake <- 
    array(uptake_by_age$central, c(length(uptake_by_age$central), 2))
  
  data_vaccination <- data_vaccination %>%
    dplyr::rename(dose1 = first_dose,
                  dose2 = second_dose,
                  dose3 = third_dose,
                  dose4 = fourth_dose)
  n_doses <- 2
  dose_start_dates <- c("2020-12-08",
                        "2020-12-08")
  vaccination <- 
    spimalot::spim_vaccination_data(date, region, vaccine_uptake, 
                                    vaccine_days_to_effect, data_vaccination,
                                    n_doses, dose_start_dates,
                                    carehomes = TRUE)

  vaccine_schedule_real <- vaccination$schedule
  
  ## shift doses to account for time between dose and effect of dose
  vaccine_schedule_effect <- shift_doses(vaccine_schedule_real,
                                         vaccine_days_to_effect)

  n_doses <- vaccination$schedule$n_doses
  vaccine_index_booster <- 4L
  vaccine_index_dose2 <- 2L
  
  ## Initial seeding: seed 10 per million over 1 day (4 steps)
  initial_seed_size <- 10e-6 * sum(sircovid:::sircovid_population(region))
  initial_seed_pattern <- c(1, 1, 1, 1)
  
  ## TODO: We save epoch date as a sircovid date, but restart date as
  ## a string; we should pick one, but this does have some
  ## implications for onward functions.

  ## TODO: Some of these need to be present for the parameter construction
  ## to work but we will replace later.
  ret <- list(
    date = date,
    model_type = model_type,
    region = region,
    restart_date = restart_date,
    epoch_dates = sircovid::sircovid_date(epoch_dates),
    pillar2_age_bands = pillar2_age_bands,
    
    beta_date = sircovid::sircovid_date(beta_date),
    beta_names = beta_names,
    severity_data = severity_data,
    progression_data = progression_data,
    sens_and_spec = sens_and_spec,
    initial_seed_size = initial_seed_size,
    initial_seed_pattern = initial_seed_pattern,
    natural_waning_rate = natural_waning_rate,
    
    cross_immunity_wildtype_alpha = cross_immunity_wildtype_alpha,
    cross_immunity_alpha_delta = cross_immunity_alpha_delta,
    cross_immunity_delta_omicron = cross_immunity_delta_omicron,
    rel_severity_alpha_delta = rel_severity_alpha_delta,
    rel_severity_wildtype_alpha = rel_severity_wildtype_alpha,
    rel_severity_delta_omicron = rel_severity_delta_omicron,
    severity_cross_multiplier_delta = severity_cross_multiplier_delta,
    severity_cross_multiplier_omicron = severity_cross_multiplier_omicron,
    rel_gamma_alpha = rel_gamma_alpha,
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
    ##
    #Thom adds
    prop_pfizer = prop_pfizer,
    average_vacc_efficacy_alpha = average_vacc_efficacy_alpha,
    average_vacc_efficacy_delta = average_vacc_efficacy_delta,
    #
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
