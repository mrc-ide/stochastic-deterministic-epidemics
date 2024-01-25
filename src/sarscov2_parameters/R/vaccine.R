shift_doses <- function(vaccine_schedule, vaccine_days_to_effect) {
  ## This function shifts, for each dose j, all jth doses in vaccine_schedule by
  ## vaccine_days_to_effect[j]
  n_groups <- dim(vaccine_schedule$doses)[1]
  n_doses <- dim(vaccine_schedule$doses)[2]
  n_days <- dim(vaccine_schedule$doses)[3]
  
  shift1 <- function(i, j) {
    c(rep(0, vaccine_days_to_effect[j]),
      vaccine_schedule$doses[i,j,],
      rep(0, max(vaccine_days_to_effect) - vaccine_days_to_effect[j]))
  }
  
  schedule_effect <-
    vapply(seq_len(n_groups),
           function (i) {
             vapply(seq_len(n_doses), function(j) shift1(i, j),
                    numeric(n_days + max(vaccine_days_to_effect)))},
           array(0, c(n_days + max(vaccine_days_to_effect), n_doses)))
  schedule_effect <- aperm(schedule_effect, c(3, 2, 1))
  
  vaccine_schedule$doses <- schedule_effect
  vaccine_schedule
}

calculate_average_vacc_efficacy <- function(vaccine_efficacy, prop_pfizer) {
  
  # unvaccinated / 1st dose no eff / 1st dose full eff / 2nd dose full eff / waned / boosted
  vaccine_efficacy$week_wane <- NULL
  ret <- vaccine_efficacy %>%
    tidyr::pivot_longer(-c(type, vaccine, dose), names_to = "analysis") %>%
    tidyr::expand_grid(group = seq_len(19)) %>%
    # code dose numbers according to their vaccine strata
    dplyr::mutate(prop_pfizer = prop_pfizer[group],
                  stratum = forcats::fct_recode(as.character(dose),
                                                "stratum_3" = "1",
                                                "stratum_4" = "2",
                                                "stratum_5" = "waned",
                                                "stratum_6" = "booster",
                                                "stratum_7" = "booster_waned")) %>%
    dplyr::select(-dose) %>%
    tidyr::pivot_wider(names_from = stratum) %>%
    # add in strata 1 in which there is 0 vaccine protection
    dplyr::mutate(stratum_1 = 0, .before = stratum_3) %>%
    # add in strata 3 in which there is 0 vaccine protection for 3 weeks after dose 1
    dplyr::mutate(stratum_2 = 0, .before = stratum_3) %>%
    tidyr::pivot_longer(dplyr::starts_with("stratum"), names_to = "stratum") %>%
    tidyr::pivot_wider(names_from = vaccine) %>%
    # calculate combined vaccine efficacy based on % of PF
    dplyr::mutate(efficacy = AZ * (1 - prop_pfizer) + PF * prop_pfizer) %>%
    dplyr::select(-c(AZ, PF)) %>%
    tidyr::pivot_wider(names_from = type, values_from = efficacy) %>%
    dplyr::arrange( analysis, stratum, group)
  
  ## split by analysis type
  split(ret, ret$analysis)
}

get_vaccine_conditional_prob <- function(eff_death,
                                         eff_severe_disease,
                                         eff_disease, eff_infection,
                                         eff_onwards_transmission = NULL) {
  
  n_group <- 19
  rel_susceptibility <- matrix(1 - eff_infection, n_group)
  rel_p_sympt <- matrix(1 - eff_disease, n_group) / rel_susceptibility
  rel_p_hosp_if_sympt <-
    matrix(1 - eff_severe_disease, n_group) / (rel_susceptibility * rel_p_sympt)
  rel_p_death <-
    matrix(1 - eff_death, n_group) / (rel_susceptibility * rel_p_sympt * rel_p_hosp_if_sympt)
  res <- list(rel_susceptibility = rel_susceptibility,
              rel_p_sympt = rel_p_sympt,
              rel_p_hosp_if_sympt = rel_p_hosp_if_sympt,
              rel_p_death = rel_p_death)
  
  if(!is.null(eff_onwards_transmission)) {
    res$rel_infectivity <- matrix(1 - eff_onwards_transmission, n_group)
  } else {
    res$rel_infectivity <- res$rel_susceptibility
    res$rel_infectivity[] <- 1
  }
  
  ## This caused problems because of rounding.
  ## never interested in efficacy beyond 0.1% hence the rounding
  if (any(signif(unlist(res), 3) < 0) |
      any(signif(unlist(res), 3) > 1)) {
    stop("incompatible vaccine efficacy parameters")
  }
  
  res
}

add_vaccine_efficacy_booster_waned <- function(vaccine_efficacy) {
  ve_booster_waned <- vaccine_efficacy %>% 
    pivot_longer(!c("type", "vaccine", "dose"), names_to = "scenario") %>%
    mutate(dose = case_when(dose == "1" ~ "dose1",
                            dose == "2" ~ "dose2",
                            TRUE ~ dose)) %>%
    pivot_wider(names_from = dose) %>%
    filter(vaccine == "PF") %>%
    ## We want to ensure ((1 - VE[booster_waned]) / VE[booster_waned]) / 
    ##                     ((1 - VE[booster]) / VE[booster]) =
    ##                       ((1 - VE[waned]) / VE[waned]) / 
    ##                         ((1 - VE[dose2]) / VE[dose2])
    mutate(booster_waned = 1 / (1 + 
                                  ((1 - booster) / booster) * 
                                  ((1 - waned) / waned) / ((1 - dose2) / dose2))) %>%
    pivot_longer(!c(type, vaccine, scenario), names_to = "dose") %>%
    pivot_wider(names_from = scenario) %>%
    filter(dose == "booster_waned")
  
  ve_booster_waned <- as.data.frame(ve_booster_waned)
  
  ## In case this has produced VE against a later stage of infection being LOWER
  ## than VE against an earlier stage, we make sure VE against each stage is
  ## at least equal to the VE against the previous stage
  scenarios <- !(names(ve_booster_waned) %in% c("type", "vaccine", "dose"))
  ve_booster_waned[ve_booster_waned$type == "disease", scenarios] <-
    pmax(ve_booster_waned[ve_booster_waned$type == "infection", scenarios], 
         ve_booster_waned[ve_booster_waned$type == "disease", scenarios])
  ve_booster_waned[ve_booster_waned$type == "severe_disease", scenarios] <-
    pmax(ve_booster_waned[ve_booster_waned$type == "disease", scenarios], 
         ve_booster_waned[ve_booster_waned$type == "severe_disease", scenarios])
  ve_booster_waned[ve_booster_waned$type == "death", scenarios] <-
    pmax(ve_booster_waned[ve_booster_waned$type == "death", scenarios], 
         ve_booster_waned[ve_booster_waned$type == "severe_disease", scenarios])
  
  ## We assume all booster efficacy (including waned) is that of PF
  vaccine_efficacy <- rbind(vaccine_efficacy, ve_booster_waned)
  ve_booster_waned$vaccine <- "AZ"
  vaccine_efficacy <- rbind(vaccine_efficacy, ve_booster_waned)
  vaccine_efficacy
}
