
get_age_band_labels <- function(default, bands = NULL) {
  
  if (default) {
    
    age_bands <- paste0(c(seq(0, 80, 5)), "_", c(seq(4, 79, 5), 120))
    
  } else {
    
    age_bands <- bands
  }
  
  age_bands <- gsub("_", " to ", age_bands)
  age_bands <- gsub("to 120", "plus", age_bands)
  
  age_bands
}


get_strain_timelines <- function(fits_data) {
  
  # Extract VOC frequency we fitted to
  ret <- NULL
  
  for (i in names(fits_data)) {
    tmp <- fits_data[[i]]$fitted %>%
      select(date = date_string, strain_tot, strain_non_variant) %>%
      filter(date >= as.Date("2020-09-17")) %>%
      mutate(region = paste(i)) %>%
      pivot_longer(!c(date, region), names_to = "type", values_to = "value") %>%
      pivot_wider(names_from = type, values_from = value) %>%
      mutate(value = (strain_tot - strain_non_variant) / strain_tot) %>%
      select(date, region, value) %>%
      mutate(strain = case_when(
        date >= as.Date("2020-09-17") & date <= as.Date("2021-03-01") ~ "Alpha",
        date >= as.Date("2021-03-08") & date <= as.Date("2021-07-31") ~ "Delta",
        date >= as.Date("2021-11-20") & date <= as.Date("2022-01-15") ~ "Omicron"
      )) %>%
      filter(!is.na(strain)) %>%
      pivot_wider(names_from = strain, values_from = value)
    
    ret <- rbind(ret, tmp)
  }
  ret
}


get_model_demography <- function(dat, age_bands, r, what) {
 
  model_demo <- dat[["model_demography"]][[r]]
  date <- sircovid::sircovid_date_as_date(model_demo[[1]]$date)
  
  if (grepl("death", what)) {
    
    tmp <- model_demo[["deaths_hosp"]]$mean_t + model_demo[["deaths_comm"]]$mean_t
    
    df <- tmp %>%
      data.frame(date, .) %>%
      pivot_longer(!date) %>%
      mutate(name = gsub("D_hosp_", "", name)) %>%
      mutate(name = sub("_[^_]+$", "", name)) %>%
      filter(!name %in% c("CHW", "CHR")) %>%
      group_by(date, name) %>%
      summarise(value = sum(value, na.rm = TRUE)) %>%
      mutate(name = factor(name,
                           levels = c("0", "5", "10", "15", "20", "25", "30",
                                      "35", "40", "45", "50", "55", "60",
                                      "65", "70", "75", "80"),
                           labels = age_bands)) %>%
      arrange(date, name) %>%
      ungroup()
    
  } else {
    tmp <- model_demo[[what]]$mean_t[, c(1:17)]
    df <- data.frame(
      date = date,
      tmp)
    colnames(df) <- c("date", age_bands)
    
    df <- df %>%
      pivot_longer(!date)
    
    df$name <- factor(df$name,
                      levels = rev(unique(df$name)),
                      labels = rev(age_bands),
                      ordered = is.ordered(df$name))
  }
  
  df
}


get_model_trajectory <- function(dat, r, what, mean_quantiles = FALSE) {
  
  trajectories <- dat$samples[[r]]$trajectories
  state <- trajectories$state[, , -1]
  dates <- sircovid::sircovid_date_as_date(trajectories$date[-1])
  
  if (what == "admissions") {
    
    df <- state["admitted_inc", , ] + state["diagnoses_inc", , ]
    
  } else if (what %in% c("vaccine_status", "vacc_uptake")) {
    
    rows <- grep(what, rownames(state))
    df <- state[rows, , ]
    df <- apply(df, c(1, 3), mean)
    
  } else if (what == "hosp_occup") {
    
    df <- state["general", , ] + state["icu", , ]
    
  } else if (what == "diagnoses_vacc") {
    
    df <- list()
    df[["Unvaccinated"]] <- state["diagnoses_admitted_vacc_0", , ]
    df[["Dose 1"]] <- state["diagnoses_admitted_vacc_1", , ] +
      state["diagnoses_admitted_vacc_2", , ]
    df[["Dose 2"]] <- state["diagnoses_admitted_vacc_3", , ] +
      state["diagnoses_admitted_vacc_4", , ]
    df[["Booster"]] <- state["diagnoses_admitted_vacc_5", , ] +
      state["diagnoses_admitted_vacc_6", , ]
    
  } else if (what == "D_hosp_vacc") {
    
    df <- list()
    df[["Unvaccinated"]] <- state["D_hosp_vacc_0", , ]
    df[["Dose 1"]] <- state["D_hosp_vacc_1", , ] +
      state["D_hosp_vacc_2", , ]
    df[["Dose 2"]] <- state["D_hosp_vacc_3", , ] +
      state["D_hosp_vacc_4", , ]
    df[["Booster"]] <- state["D_hosp_vacc_5", , ] +
      state["D_hosp_vacc_6", , ]
    
  } else {
    
    df <- state[what, , ]
    
  }
  
  if (mean_quantiles) {
    
    get_mean_quant <- function (x, diff = FALSE) {
      if (diff) {
        data.frame(
          date = dates,
          mean = c(0, diff(colMeans(x))),
          lb = c(0, diff(matrixStats::colQuantiles(x, probs = 0.025))),
               
          ub = c(0, diff(matrixStats::colQuantiles(x, probs = 0.975))))
      } else {
        data.frame(
          date = dates,
          mean = colMeans(x),
          lb = matrixStats::colQuantiles(x, probs = 0.025),
          ub = matrixStats::colQuantiles(x, probs = 0.975)) 
      }
    }
    
    if (is.list(df)) {
      nm <- names(df)
      df <- lapply(df, get_mean_quant, TRUE)
      names(df) <- nm
    } else {
      df <- get_mean_quant(df)
    }
  }
  df
}


get_severity_excl_immunity <- function(dat, what, r, dates_vect) {

  # Vector of model steps/dates
  step <- dat$severity[[1]]$step
  dates <- step / 4
  
  out <- NULL
  
  if (r == "england") {
    reg <- sircovid::regions("england")
  } else {
    reg <- r
  }
  
  for (r in reg) {
    
    pars <- dat$samples[[r]]$pars
    transform <- dat$samples[[r]]$predict$transform
    
    step_vect <- step[match(dates_vect, dates)]
    
    sev_vector <- function(x) {
      mean <- mean(x)
      lb <- quantile(x, 0.025)
      ub <- quantile(x, 0.975)
      
      unname(c(mean, lb, ub))
    }
    
    get_variants <- function(j) {
      pars_model <- lapply(spimalot:::seq_rows(pars),
                           function(i) transform(pars[i, ])[[j]]$pars)
      
      sircovid::lancelot_ifr_excl_immunity(step_vect, pars_model)
    }
    
    # Get variants' intrinsic severity
    Wildtype_Alpha <- get_variants(3)
    Delta_Omicron <- get_variants(6)
    
    tmp <- NULL
    for (i in seq_along(dates_vect)) {
      for (w in what) {
        tmp_w <- data.frame(
          estimate = c("mean", "lb", "ub"),
          Wildtype = sev_vector(Wildtype_Alpha[[w]][i, 1, ]),
          Alpha = sev_vector(Wildtype_Alpha[[w]][i, 2, ]),
          Delta = sev_vector(Delta_Omicron[[w]][i, 1, ]),
          Omicron = sev_vector(Delta_Omicron[[w]][i, 2, ])
        ) %>%
          pivot_longer(!estimate) %>%
          mutate(period = names(dates_vect[i])) %>%
          mutate(region = r,
                 source = w) 
      
        tmp <- rbind(tmp, tmp_w)
      }
    }
      
    out <- rbind(out, tmp)
  }
  out$period <- factor(out$period, levels = unique(names(dates_vect)))
  
  out
  
}


get_r0_region <- function(dat, r) {
  
  get_region_variants <- function(r) {
    pars <- dat$samples[[r]]$pars
    
    # Take regional Rt_general at the start of the model for Wildtype R0
    # then parametrically calculate for the VOCs
    Wildtype <- dat$rt[[r]]$Rt_general[1, "weighted", ]
    Alpha <- Wildtype * pars[, "ta_alpha"]
    Delta <- Alpha * pars[, "ta_delta"]
    Omicron <- Delta * pars[, "ta_omicron"]
    
    list(Wildtype = Wildtype,
         Alpha = Alpha,
         Delta = Delta,
         Omicron = Omicron)
  }
  
  if (r == "england") {
    reg <- sircovid::regions("england")
    r0 <- lapply(reg, get_region_variants)
    r0 <- spimalot:::list_transpose(r0)
    
    get_reg_pop <- function(r) {
      p <- dat$samples[[r]]$predict$transform(dat$samples[[r]]$pars[1, ])
      sum(p[[1]]$pars$population)
    }
    
    wts <- vapply(reg, get_reg_pop, numeric(1))
    
    r0 <- lapply(r0, function(x) {
      apply(spimalot:::abind_quiet(x, along = 2), 1, weighted.mean, w = wts)})
  } else {
    r0 <- get_region_variants(r)
  }
  
  summary_vector <- function(x) {
    mean <- mean(x)
    lb <- quantile(x, 0.025)
    ub <- quantile(x, 0.975)
    
    unname(c(mean, lb, ub))
  }
  
  data.frame(
    Wildtype = summary_vector(r0$Wildtype),
    Alpha = summary_vector(r0$Alpha),
    Delta = summary_vector(r0$Delta),
    Omicron = summary_vector(r0$Omicron)
  ) %>%
  mutate(value = c("mean", "lb", "ub"),
         region = r,
         source = "r0")
    
}


get_model_severity <- function(severity, vam_data, strain_epochs, what,
                               by_strain = FALSE, by_age = FALSE,
                               by_vacc_class = FALSE, classes = NULL,
                               cut_off = 0.5) {
  age_bands <- seq(0, 80, 5)
  
  # Define severity trajectories to extract
  if (by_strain) {
    get <- c(what, paste0(what, c("_strain_1", "_strain_2")))
  } else if (by_age) {
    
    get <- paste0(what, "_age_", age_bands)
    
  } else if (by_vacc_class) {
    
    get <- paste0(what, "_disag_", age_bands)
    get <- outer(get, c("", paste0("_", seq(1, 6, 1))), paste0)
    colnames(get) <- c("unvacc", "1_no_eff", "1_full_eff", "2_dose",
                       "2_dose_wane", "boost", "boost_wane")
    
    get <- as.vector(get[, classes])
    
  } else {
    get <- what
  }
  
  waves <- c("Wildtype", "Alpha", "Delta", "Omicron")
  dates <- data.frame(
    date = sircovid::sircovid_date_as_date(severity[[1]]$date[-1]))
  
  ret <- NULL
  
  suppressMessages(
    for (i in names(severity)) {
      
      voc <- vam_data %>%
        filter(region == i) %>%
        select(date, Alpha, Delta, Omicron)

      voc_dates <- c(voc$date[min(which(voc$Alpha > cut_off))],
                     voc$date[min(which(voc$Delta > cut_off))],
                     voc$date[min(which(voc$Omicron > cut_off))]
                     )
      
      tmp <- data.frame(date = dates,
                        region = i) %>%
        mutate(dominant_voc = case_when(
          date < voc_dates[1] ~ waves[[1]],
          date >= voc_dates[1] & date < voc_dates[2] ~ waves[[2]],
          date >= voc_dates[2] & date < voc_dates[3] ~ waves[[3]],
          date >= voc_dates[3] ~ waves[[4]]))
      
      for (g in get) {
        sev <- severity[[i]][[g]][-1, ]
        
        tmp1 <- data.frame(t(apply(sev, 1, quantile, c(0.5, 0.025, 0.975),
                                   na.rm = TRUE))) %>%
          setNames(., c(paste0(g, "_mean"), paste0(g, "_lb"), paste0(g, "_ub")))
        
        tmp <- cbind(tmp, tmp1)
      }
     
      ret <- rbind(ret, tmp)
    }
  )
  ret$dominant_voc <- factor(ret$dominant_voc, levels = waves)
  
  # NA strain_1 values on strain rotation epochs
  strain_1_columns <- grep("strain_1", colnames(ret), value = TRUE)
  ret[ret$date %in% strain_epochs, strain_1_columns] <- NA_real_
  
  # VOCs are seeded in the 15_19 age susceptible compartment. NA stain_2 values
  # while new VOC gets properly seeded across the population to avoid underestimating
  # severity. Hardcoded as this is VOC specific.
  strain_2_columns <- grep("strain_2", colnames(ret), value = TRUE)
  alpha_date_na <- seq.Date(strain_epochs[1], strain_epochs[1] + floor(7 * 6.5), by = "day")
  delta_date_na <- seq.Date(strain_epochs[2], strain_epochs[2] + floor(7 * 4.5), by = "day")
  omicron_date_na <- seq.Date(strain_epochs[3], strain_epochs[3] + floor(7 * 3.5), by = "day")
  
  ret[ret$date %in% c(alpha_date_na, delta_date_na, omicron_date_na),
      strain_2_columns] <- NA_real_
  
  ret
}


get_national_intrinsic_values <- function(dat, p) {
  
  int_sev <- dat$intrinsic_severity %>%
    filter(period == p) %>%
    filter(region == "england")
  
  r0 <- get_r0_region(dat, "england") %>%
    mutate(estimate = value) %>%
    select(!value) %>%
    pivot_longer(!c(estimate, region, source)) %>%
    select(estimate, name, value, region, source)
    
  tmp <- int_sev %>%
    select(!period) %>%
    rbind(., r0)
    
  out <- NULL
  for (v in unique(tmp$name)) {
    
    ret <- NULL
    for (i in unique(tmp$source)) {
      mean <- round(sum((tmp %>% filter(source == i, estimate == "mean", name == v))$value), 4)
      lb <- round(sum((tmp %>% filter(source == i, estimate == "lb", name == v))$value), 4)
      ub <- round(sum((tmp %>% filter(source == i, estimate == "ub", name == v))$value), 4)
      
      x <- data.frame(mean = mean, lb = lb, ub = ub)
      
      ret[[i]] <- unlist(x)
    }
    out[[v]] <- ret
  }
  out
}


vnapply <- function (.x, .fun, ...) 
{
  vapply(.x, .fun, numeric(1), ...)
}


get_fitted_par <- function(dat, regions = NULL, pars_to_plot) {
  
  if (is.null(regions)) {regions <- names(dat$samples)}
  
  samples <- dat$samples[regions]
  date <- dat$info$date
  
  if ("start_date" %in% colnames(samples[[1]]$pars)) {
    # order by start_date
    mean_start_date <-
      vnapply(samples, function(x) mean(x$pars[, "start_date"]))
    
    ordered_regions <- c(names(sort(mean_start_date[regions],
                                    decreasing = TRUE)))
  } else {
    ordered_regions <-sort(regions, decreasing = TRUE)
  }
  
  samples <- samples[ordered_regions]
  par_names <- colnames(samples[[1]]$pars)
  
  ## Setting maximum value to plot
  par_max <- rep(NA, length(par_names))
  names(par_max) <- par_names
  pars_p_NC <- grep("p_NC", names(par_max), value = TRUE)
  par_max[pars_p_NC] <- 0.01
  par_max["rho_pillar2_tests"] <- 0.02
  
  n_regions <- length(samples)
  
  extract_sample <- function(par_name) {
    lapply(samples, function(x) as.numeric(x$pars[, par_name]))
  }
  
  
  get_par <- function(par_name) {
    
    par <- extract_sample(par_name)
    par_info <- subset(pars_info, pars_info$name == par_name)
    
    if (is.na(par_max[[par_name]])) {
      xmax <- max(par_info$max)
    } else {
      xmax <- par_max[[par_name]]
    }
    
    xmin <- min(par_info$min)
    
    par <- t(sapply(par, function(p) c(mean(p), quantile(p, 0.025), quantile(p, 0.975))))
    par <- as.data.frame(par) %>% tibble::rownames_to_column(.) %>%
      mutate(prior_min = xmin, prior_max = xmax) %>%
      magrittr::set_colnames(c("region", "mean", "lb", "ub", "prior_min", "prior_max"))
    
    if (grepl("date", par_name)) {
      par <- par %>% mutate_if(is.numeric, sircovid::sircovid_date_as_date)
    }
    
    par
  }
  
  hps <- dat$parameters$prior
  pars_info <- dat$parameters$info
  
  if ("start_date" %in% pars_to_plot) {
    pars_to_plot <- c("start_date", pars_to_plot[pars_to_plot != "start_date"])
  }
  
  x <- lapply(pars_to_plot, get_par)
  names(x) <- pars_to_plot
  
  x
}

forest_plot_labels <- function(dat) {
  par_names <- unique(dat$parameters$proposal$name)
  beta_date <- sircovid::sircovid_date_as_date(dat$samples[[1]]$info$beta_date)
  
  n_betas <- length(beta_date)
  
  labels <- lapply(seq_len(n_betas), 
                   function(x) {
                     date_x <- as.character(beta_date[x])
                     bquote(beta[.(x)] ~ (.(date_x)))
                   })
  names(labels) <- paste0("beta", seq_len(n_betas))
  
  plus <- "+"
  labels <- c(labels,
              alpha_admission = expression(alpha[A]),
              alpha_D = expression(alpha[D[comm]]),
              alpha_death_hosp = expression(alpha[D[hosp]]),
              alpha_H = expression(alpha[H]),
              mu_D = expression(mu["D,1"]),
              mu_D_2 = expression(mu["D,2"]),
              mu_D_3 = expression(mu["D,3"]),
              mu_D_4 = expression(mu["D,4"]),
              mu_D_5 = expression(mu["D,5"]),
              mu_gamma_H = expression(mu[paste(gamma[H], ",", 1)]),
              mu_gamma_H_2 = expression(mu[paste(gamma[H], ",", 2)]),
              mu_gamma_H_3 = expression(mu[paste(gamma[H], ",", 3)]),
              mu_gamma_H_4 = expression(mu[paste(gamma[H], ",", 4)]),
              p_G_D = expression(p[paste(G[D], ",", 1)]^max),
              p_G_D_2 = expression(p[paste(G[D], ",", 2)]^max),
              p_H = expression(p[H]^max),
              p_H_2 = expression(p["H,2"]^max),
              p_H_D = expression(p[H[D]]^max),
              p_ICU = expression(p["ICU,1"]^max),
              p_ICU_2 = expression(p["ICU,2"]^max),            
              p_ICU_D = expression(p[ICU[D]]^max),
              p_NC_15_24 = expression(p[NC]^"[15,25)"),
              p_NC_25_49 = expression(p[NC]^"[25,50)"),
              p_NC_50_64 = expression(p[NC]^"[50,65)"),
              p_NC_65_79 = expression(p[NC]^"[65,80)"),
              p_NC_80_plus = expression(p[NC]^"[80+)"),
              p_NC_weekend_15_24 = expression(p[NCW]^"[15,25)"),
              p_NC_weekend_25_49 = expression(p[NCW]^"[25,50)"),
              p_NC_weekend_50_64 = expression(p[NCW]^"[50,65)"),
              p_NC_weekend_65_79 = expression(p[NCW]^"[65,80)"),
              p_NC_weekend_80_plus = expression(p[NCW]^"[80+)"),
              p_W_D = expression(p[W[D]]^max),
              rel_p_D_alpha = expression(pi[D]^"Alpha/Wildtype"),
              rel_p_D_delta = expression(pi[D]^"Delta/Alpha"),
              rel_p_D_omicron = expression(pi[D]^"Omicron/Delta"),
              rel_p_H_alpha = expression(pi[H]^"Alpha/Wildtype"),
              rel_p_H_delta = expression(pi[H]^"Delta/Alpha"),
              rel_p_H_omicron = expression(pi[H]^"Omicron/Delta"),
              rel_p_ICU_alpha = expression(pi[ICU]^"Alpha/Wildtype"),
              rel_p_ICU_delta = expression(pi[ICU]^"Delta/Alpha"),
              rel_p_ICU_omicron = expression(pi[ICU]^"Omicron/Delta"),
              rho_pillar2_tests = expression(rho[P2[test]]),
              seed_date_alpha = expression(t["Alpha"]),
              seed_date_delta = expression(t["Delta"]),
              seed_date_omicron = expression(t["Omicron"]),
              start_date = expression(t["Wildtype"]),
              ta_alpha = expression(sigma["Alpha/Wildtype"]),
              ta_delta = expression(sigma["Delta/Alpha"]),
              ta_omicron = expression(sigma["Omicron/Delta"])
  )
  labels
}

get_convergence_diagnostic <- function(dat) {
  
  conv_dx <- function(sample) {
    n_full_pars <- nrow(sample$full_pars)
    n_chains <- max(sample$chain)
    
    sample$chain_full <- rep(seq_len(n_chains), each = n_full_pars / n_chains)
    
    chains <- lapply(unname(split(data.frame(sample$full_pars),
                                  sample$chain_full)), coda::as.mcmc)
    
    rhat <- tryCatch(coda::gelman.diag(chains), error = function(e) NULL)
    if (!is.null(rhat)) {
      rhat <- round(max(rhat$psrf[, "Point est."]), 2)
    } else {
      rhat <- NA_real_
    }
    
    ess <- function(p) {
      traces <- matrix(p, ncol = n_chains)
      sum(coda::effectiveSize(coda::as.mcmc(traces)))
    }
    
    pars <- sample$full_pars
    nms <- colnames(pars)
    pars_ess <- lapply(nms, function (nm) {
      ess(pars[, nm])
    })
    pars_ess <- round(min(unlist(pars_ess)))
    
    data.frame(rhat, pars_ess)
  }
  
  ret <- NULL
  
  for (r in sircovid::regions("england")) {
    tmp <- conv_dx(dat$samples[[r]])
    tmp <- cbind(data.frame(region = r), tmp)
    ret <- rbind(ret, tmp)
  }
  
  ret
}
