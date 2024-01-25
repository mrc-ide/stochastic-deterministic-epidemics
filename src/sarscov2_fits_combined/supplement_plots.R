suppl_regional_intrinsic <- function(dat, r, which, strain_epochs,
                                     vam_data, cut_off = 0.9) {
  
  # Plot intrinsic severity by region
  int_sev <- plot_intrinsic_severity_vs_transmissibility(
    dat, r, c("IHR", "HFR", "IFR"), which)
  
  # Make sure bottom panel has x axis text
  int_sev[["IFR"]] <- int_sev[["IFR"]] +
    labs(x = expression(R[0])) +
    theme(axis.text.x = element_text(size = 10))
  
  (int_sev[["IHR"]] / int_sev[["HFR"]] / int_sev[["IFR"]]) +
    
    # (col_1 | p2) +
    plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 10))
  
}


suppl_compare_demography <- function(dat, r) {
  
  age_bands <- get_age_band_labels(TRUE)
  vax_ro <- as.Date("2020-12-08")
  
  ## Process data and model admissions by age
  data_adm <- dat$data[[r]]$full %>%
    select(date_string, contains("all_admission_")) %>%
    mutate(date = lubridate::floor_date(date_string, "weeks", week_start = 3)) %>%
    group_by(date) %>%
    summarise_if(is.numeric, sum)
  colnames(data_adm) <- gsub("all_admission_", "", colnames(data_adm))
  
  model_demo_adm <- get_model_demography(dat, age_bands, r, "admissions") %>%
    pivot_wider() %>%
    mutate(`0_9` = rowSums(.[2:3])) %>%
    mutate(`10_19` = rowSums(.[4:5])) %>%
    mutate(`20_29` = rowSums(.[6:7])) %>%
    mutate(`30_39` = rowSums(.[8:9])) %>%
    mutate(`40_49` = rowSums(.[10:11])) %>%
    mutate(`50_59` = rowSums(.[12:13])) %>%
    mutate(`60_69` = rowSums(.[14:15])) %>%
    mutate(`70_79` = rowSums(.[16:17])) %>%
    mutate(`80_plus` = rowSums(.[18])) %>%
    mutate(week = lubridate::floor_date(date, "weeks", week_start = 3)) %>%
    group_by(week) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(date = as.Date(week))
  
  model_demo_adm <- model_demo_adm[, colnames(data_adm)] %>%
    pivot_longer(!date)
  model_demo_adm$name <- forcats::fct_rev(model_demo_adm$name)
  
  data_adm <- data_adm %>% pivot_longer(!date)
  data_adm$name <- forcats::fct_rev(data_adm$name)
  
  
  ## Process data and model deaths by age
  
  # Add hospital and community deaths together
  data_death <- dat$data[[r]]$full %>%
    select(date = date_string, contains("deaths_hosp_"))
  data_comm <- dat$data[[r]]$full %>%
    select(date = date_string, contains("deaths_comm_"))
  
  data_death[, 2:9] <- data_death[, 2:9] + data_comm[, 2:9]
  
  data_death <- data_death %>%
    mutate(date = lubridate::floor_date(date, "weeks", week_start = 3)) %>%
    group_by(date) %>%
    summarise_if(is.numeric, sum)
  colnames(data_death) <- gsub("deaths_hosp_", "", colnames(data_death))
  
  
  model_demo_death <- get_model_demography(dat, age_bands, r, "death") %>%
    pivot_wider()
  colnames(model_demo_death) <- gsub(" to ", "_", colnames(model_demo_death))
  colnames(model_demo_death) <- gsub(" ", "_", colnames(model_demo_death))
  model_demo_death <- model_demo_death %>%
    mutate(`0_49` = rowSums(.[2:11])) %>%
    mutate(week = lubridate::floor_date(date, "weeks", week_start = 3)) %>%
    group_by(week) %>%
    summarise_if(is.numeric, sum) %>%
    mutate(date = as.Date(week))
  model_demo_death <- model_demo_death[, colnames(data_death)] %>%
    pivot_longer(!date)
  model_demo_death$name <- forcats::fct_rev(model_demo_death$name)
  
  data_death <- data_death %>% pivot_longer(!date)
  data_death$name <- forcats::fct_rev(data_death$name)
  
  
  # Age distribution plots
  adm_age <- gg_trajectory_age(model_demo_adm, "Model (%)", vax_ro,
                               legend.tit = "Hospital admissions") +
    geom_vline(xintercept = seq(as.Date("2020-03-01", format = "%Y-%m-%d"),
                                by = "month", length.out = 25),
               col = "grey95", linewidth = 0.1) +
    geom_vline(xintercept = vax_ro, lty = 2) +
    theme(axis.ticks = element_line())
  
  
  data_adm_age <- gg_trajectory_age(data_adm, "Data (%)", vax_ro) +
    geom_vline(xintercept = seq(as.Date("2020-03-01", format = "%Y-%m-%d"),
                                by = "month", length.out = 25),
               col = "grey95", linewidth = 0.1) +
    geom_vline(xintercept = vax_ro, lty = 2) +
    theme(legend.position = "none", axis.ticks = element_line())
  
  death_age <- gg_trajectory_age(model_demo_death, "Model (%)",
                                 vax_ro, option = 'A', legend.tit = "All deaths (hosp + comm)") +
    geom_vline(xintercept = seq(as.Date("2020-03-01",
                                        format = "%Y-%m-%d"),
                                by = "month", length.out = 25), col = "grey95", size = 0.1) +
    geom_vline(xintercept = vax_ro, lty = 2) +
    theme(axis.ticks = element_line())
  
  data_death_age <- gg_trajectory_age(data_death, "Data (%)", vax_ro, option = 'A') +
    geom_vline(xintercept = seq(as.Date("2020-03-01",
                                        format = "%Y-%m-%d"),
                                by = "month", length.out = 25), col = "grey95", size = 0.1) +
    geom_vline(xintercept = vax_ro, lty = 2) +
    theme(legend.position = "none", axis.ticks = element_line())
  
  (adm_age + data_adm_age) / (death_age + data_death_age)
}


suppl_emergence_demography <- function(dat, r) {
  
  age_bands <- get_age_band_labels(TRUE)
  emergence_date <- dat$parameters$base[[1]]$intrinsic_severity_dates["Emergence3"]
  emergence_date <- sircovid:::sircovid_date_as_date(emergence_date)
  
  ## Extract model demography on date of Emergence3
  model_demo_adm <- get_model_demography(dat, age_bands, r, "admissions") %>%
    filter(date == emergence_date) %>%
    mutate(prop = value/sum(value))
  model_demo_adm$name <- factor(model_demo_adm$name,
                                levels = unique(model_demo_adm$name))
  
  model_demo_death <- get_model_demography(dat, age_bands, r, "death") %>%
    filter(date == emergence_date) %>%
    mutate(prop = value/sum(value))
  model_demo_death$name <- factor(model_demo_death$name,
                                levels = unique(model_demo_death$name))
  
  # Age distribution plots
  adm <- ggplot(model_demo_adm, aes(name, prop, fill = name)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_fill_viridis_d() +
    labs(y = "Hospital admissions", x = "") +
    theme_minimal() +
    theme(axis.line = element_line(),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  death <- ggplot(model_demo_death, aes(name, prop, fill = name)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_fill_viridis_d(option = "H") +
    labs(y = "Hospital deaths", x = "") +
    theme_minimal() +
    theme(axis.line = element_line(),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 0.7))
  
  (adm + death) +
    plot_annotation(
      title = paste0("Model age distribution on date of assessment of variants' baseline severity (",
                   emergence_date, ")"))
  
}


suppl_compare_hfr <- function(dat, hfr_week, vam_data, strain_epochs) {
  
  tmp <- hfr_week %>%
    select(date, region, data = HFR, lb, ub) %>%
    mutate(data = na_if(data, 0))
  tmp$date <- as.Date(tmp$date)
  
  mod_hfr <- get_model_severity(dat$severity, vam_data, strain_epochs,
                                "hfr", cut_off = 0.75)
  
  tmp <- left_join(tmp, mod_hfr, by = c("date", "region")) %>%
    mutate(region = gsub("_", " ", region)) %>%
    mutate(region = stringr::str_to_title(region))
  
  tmp$region <- factor(tmp$region,
                       levels = c("England", unique(tmp$region)[-2]))
  
  p <- ggplot(tmp, aes(date, data)) +
    geom_point(size = 0.7, alpha = 0.5, aes(col = region)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, col = region), size = 0.1) +
    geom_ribbon(aes(ymin = hfr_lb * 100, ymax = hfr_ub * 100, fill = region),
                alpha = 0.5, show.legend = FALSE) +
    labs(x = "", y = "Effective HFR (%)") +
    facet_wrap(~region, nrow = 2) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(expand = c(0, 0))
  
  apply_theme(p) +
    theme(axis.text.x = element_text(angle = 45),
          legend.position = "none")
  
}


gg_par_plot <- function(dat, regions = NULL, pars_to_plot) {
  
  pars <- get_fitted_par(dat, regions, pars_to_plot)
  
  # Get model R0 by region and range to plot
  r0_reg <- get_r0_region(dat, "england")
  
  start_wildtype <- r0_reg %>%
    select(region, value, variant = Wildtype) %>%
    pivot_wider(names_from = value, values_from = variant) %>%
    magrittr::set_colnames(c("region", "r0_mean", "r0_lb", "r0_ub")) %>%
    left_join(., pars$start_date)
  
  seed_alpha <- r0_reg %>%
    select(region, value, variant = Alpha) %>%
    pivot_wider(names_from = value, values_from = variant) %>%
    magrittr::set_colnames(c("region", "r0_mean", "r0_lb", "r0_ub")) %>%
    left_join(., pars$seed_date_alpha)
  
  seed_delta <- r0_reg %>%
    select(region, value, variant = Delta) %>%
    pivot_wider(names_from = value, values_from = variant) %>%
    magrittr::set_colnames(c("region", "r0_mean", "r0_lb", "r0_ub")) %>%
    left_join(., pars$seed_date_delta)
  
  seed_omicron <- r0_reg %>%
    select(region, value, variant = Delta) %>%
    pivot_wider(names_from = value, values_from = variant) %>%
    magrittr::set_colnames(c("region", "r0_mean", "r0_lb", "r0_ub")) %>%
    left_join(., pars$seed_date_omicron)
  
  gg_seed_vs_r0 <- function(df, tit) {
    priors <- c(df$prior_min[1], df$prior_max[2])
    ggplot(df, aes(mean, r0_mean, col = region)) +
      geom_point(size = 2) +
      geom_errorbar(aes(xmin = lb, xmax = ub), linewidth = 0.5) +
      geom_errorbar(aes(ymin = r0_lb, ymax = r0_ub), linewidth = 0.5) +
      geom_vline(xintercept = priors, lty = 2, col = "red3") +
      scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %y",
                   limits = c(priors[1] - 5, priors[2] + 5)) +
      labs(x = "Seed date", y = expression(R[0]), title = tit)
  }
  
  p1 <- apply_theme(gg_seed_vs_r0(start_wildtype, "Wildtype")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, size = 10, vjust = 0.7))
  p2 <- apply_theme(gg_seed_vs_r0(seed_alpha, "Alpha")) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, size = 10, vjust = 0.7))
  p3 <- apply_theme(gg_seed_vs_r0(seed_delta, "Delta")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, size = 10, vjust = 0.7))
  p4 <- apply_theme(gg_seed_vs_r0(seed_omicron, "Omicron")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, size = 10, vjust = 0.7))
  
  (p1 + p2) / (p3 + p4)
}


suppl_sev_winter_20_21 <- function(dat, england_data, hfr_data, strain_epochs,
                          what, xmin, xmax) {
  
  periods <- as.vector(unique(dat$intrinsic_severity$period)[-c(1:4, 12:15)])
  national <- NULL
  for (p in periods) {
    tmp2 <- get_national_intrinsic_values(dat, p)
    
    tmp2 <- plyr::ldply(tmp2, data.frame)
  }
  
  national <- lapply(periods, function(x) {
    get_national_intrinsic_values(dat, x)})
  names(national) <- periods
  
  national <- plyr::ldply(national, data.frame)
  national$estimate <- rep(c("mean", "lb", "ub"), nrow(national) / 3)
  
  national <- national %>%
    pivot_longer(!c('.id', estimate)) %>%
    setNames(c("period", "estimate", "name", "value")) %>%
    mutate(source = stringr::str_extract(name, "[^.]+$")) %>%
    mutate(name = stringr::str_extract(name, "[^.]+")) %>%
    filter(name %in% c("Wildtype", "Alpha"),
           source == "HFR") %>%
    pivot_wider(names_from = estimate, values_from = value)
  national$name <- factor(national$name, levels = c("Wildtype", "Alpha"))
  national$period <- factor(national$period, levels = periods)
  
  
  icu_traj <- get_model_trajectory(dat, "england", "icu", TRUE)
  icu_traj$data <- dat$data[["england"]]$full$icu
  icu_traj <- icu_traj %>%
    filter(date >= xmin, date <= xmax)

  # Some graphical parameters
  ymax <- max(icu_traj$ub)
  sec_ymax <- max(icu_traj$data)
  
  
  # Panel A: intrinsic HFR by variant (Wt & Alpha) at Oct 20, Jan 21 & Feb 21
  hfr_t <- national %>%
    ggplot(., aes(x = period, y = mean, col = name, fill = name)) +
    geom_point(shape = 18, size = 3, position = position_dodge(0.3)) +
    geom_crossbar(aes(ymin = lb, ymax = ub),
                  width = 0.5, position = position_dodge(0.3), alpha = 0.3) +
    scale_x_discrete(expand = expansion(mult = 0.001)) +
    labs(x = "", y = "Basic HFR *") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = all_variant_colours) +
    scale_color_manual(values = all_variant_colours) +
    scale_shape_manual(values = seq(1:7))
  
  hfr_t <- apply_theme(hfr_t) +
    theme(axis.text.x = element_blank(),
          legend.position = "top")
  
  
  # Panel B: MV occupancy epidemic curve
  mv_occup <- england_data %>%
    select(date, region, grep("mech_", colnames(england_data), value = TRUE)) %>%
    mutate(date = as.Date(date)) %>%
    pivot_longer(!c(date, region))
  
  mv_occup_national <- mv_occup %>%
    group_by(date, name) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(occupancy = (mech_vent_open - mech_vent_unoccupied) / mech_vent_open,
           region = "england")
  
  mv_plot <- ggplot(mv_occup_national, aes(x = date)) +
    geom_point(aes(y = occupancy, col = "Data"), alpha = 0.6, size = 1) +
    scale_color_manual(name = "", values = fit_cols[1]) +
    scale_x_date(date_breaks = "month", date_labels = "%b %y",
                 limits = c(xmin, xmax), minor_breaks = NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(y = "MV occupancy data **", x = "")
  
  mv_plot <- apply_theme(mv_plot) +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.7),
          legend.position = "none")
  
  
  eff_hfr <- get_model_severity(dat$severity, vam_data = vam_data,
                                strain_epochs = strain_epochs, "hfr",
                                by_strain = TRUE, cut_off = 0.5) %>%
    select(!c(hfr_mean, hfr_lb, hfr_ub)) %>%
    pivot_longer(!c(date, region, dominant_voc)) %>%
    mutate(variant = case_when(
      str_detect(name, "strain_1") ~ "Wildtype",
      TRUE ~ "Alpha")) %>%
    mutate(name = gsub("hfr_strain_1_", "", name)) %>%
    mutate(name = gsub("hfr_strain_2_", "", name)) %>%
    pivot_wider()
  
  eff_hfr_plot <- eff_hfr %>%
    filter(region == "england", date >= xmin, date <= xmax) %>%
    ggplot(., aes(x = date)) +
    geom_line(aes(y = mean, col = variant)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = variant), alpha = 0.3) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5), 
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_date(date_breaks = "month", date_labels = "%b %y",
                 limits = c(xmin, xmax), minor_breaks = NULL) +
    labs(y = "Effective HFR", x = "") +
    scale_color_manual(values = all_variant_colours[2:1]) +
    scale_fill_manual(values = all_variant_colours[2:1])
  
  eff_hfr_plot <- apply_theme(eff_hfr_plot) +
    theme(axis.text.x = element_blank(),
          legend.position = "none")
  
  
  mv_occup_reg <- mv_occup %>%
    group_by(date, name, region) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(occupancy = (mech_vent_open - mech_vent_unoccupied) / mech_vent_open) %>%
    select(date, region, occupancy) %>%
    rbind(., mv_occup_national[, c("date", "region", "occupancy")])
  
  hfr_occup_data <- hfr_data[, c("date", "region", "HFR")]
  hfr_occup_data$date <- as.Date(hfr_occup_data$date)
  hfr_occup_data <- hfr_occup_data %>%
    left_join(., mv_occup_reg) %>%
    left_join(., eff_hfr[, c("date", "region", "dominant_voc")]) %>%
    filter(dominant_voc %in% c("Wildtype", "Alpha")) %>%
    mutate_at(c("occupancy", "HFR"), ~na_if(., 0)) %>%
    mutate(region = gsub("_", " ", region)) %>%
    mutate(region = stringr::str_to_title(region))
    
  hfr_occup_data$Occupancy_logit <- car::logit(hfr_occup_data$occupancy)
  hfr_occup_data$HFR_logit <- car::logit(hfr_occup_data$HFR)
  hfr_occup_data$region <- relevel(factor(hfr_occup_data$region), "England")
  
  reg_plot <- ggplot(hfr_occup_data, aes(x = Occupancy_logit, y = HFR_logit, col = dominant_voc)) +
    geom_point(alpha = 0.1) +
    scale_color_manual(values = all_variant_colours) +
    geom_smooth(method = "lm") +
    labs(x = "Occupancy (logit) ***", y = "HFR (logit) ***") +
    facet_grid(rows = vars(dominant_voc), cols = vars(region)) +
    theme(legend.position = "none")
  
  fit <- lm(data = hfr_occup_data, HFR_logit ~ Occupancy_logit * dominant_voc * region)
  
  col1 <- (hfr_t / eff_hfr_plot / mv_plot)
  
  plot <- (col1 | reg_plot) +
    plot_layout(widths = c(0.4, 1)) +
    plot_annotation(caption =
"* Basic HFR assumes healthcare at the start of each monthly interval.
 ** Only data available to classify beds as with mechanical ventilation (MV) or not.
 *** HFR and occupancy data calculated from line-list (patient level) and Sit-Rep data from NHS England, respectively, and applying logit transformation of proportions.",
   tag_levels = 'A') & theme(plot.tag = element_text(size = 10))
  
  list(plot = plot,
       linear_model = fit)
}


suppl_age_heatmaps <- function (dat, r, vam_data) {
  
  # Process trajectory for vaccine status
  vax_ro <- as.Date("2020-12-08")
  date <- sircovid::sircovid_date_as_date(dat$samples[[r]]$trajectories$date)[-1L]
  mod_traj_status <- get_model_trajectory(dat, r, "vaccine_status")
  vacc_status <- apply(mod_traj_status, 2, cumsum)
  vacc_status <- rbind(rep(0, dim(vacc_status)[2]), vacc_status)
  vacc_status <- data.frame(date = date, t(vacc_status) / vacc_status[nrow(vacc_status), ])
  
  # # Dates at which coverage ≥50% for 1st, 2nd and 3rd doses
  slice_dates <- vacc_status$date[c(min(which(1 - vacc_status$vaccine_status_1 >= 0.5)),
                                    max(which(1 - vacc_status$vaccine_status_3 <= 0.5)),
                                    max(which(1 - vacc_status$vaccine_status_5 <= 0.5)))]
  
  # # Get aggregated trajectories and append FULL data
  # admissions <- get_model_trajectory(dat, r, "admissions", mean_quantiles = TRUE)
  # admissions$data <- dat$data[[r]]$full$all_admission
  # 
  # deaths <- get_model_trajectory(dat, r, "deaths_inc", mean_quantiles = TRUE)
  # deaths$data <- dat$data[[r]]$full$deaths
  # 
  # hosp_deaths <- get_model_trajectory(dat, r, "deaths_hosp_inc", mean_quantiles = TRUE)
  # hosp_deaths$data <- dat$data[[r]]$full$deaths_hosp
  # 
  # 
  # # Fitted trajectories plots
  # adm <- gg_trajectory(admissions, ymax = max(admissions$ub), 
  #                      ylab = "Hosp Admissions") +
  #   theme(legend.position = "top")
  # 
  # all_deaths <- gg_trajectory(deaths, ymax = max(deaths$ub),
  #                             ylab = "All deaths (community + hospital)") +
  #   theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
  #         legend.position = "none")
  # 
  # h_deaths <- gg_trajectory(hosp_deaths, ymax = max(hosp_deaths$ub),
  #                             ylab = "Hosp deaths") +
  #   theme(legend.position = "none")
  
  # Effective severity by age
  age_bands <- get_age_band_labels(TRUE)
  
  ihr_age <- get_eff_sev_age(dat, r, "ihr", vam_data, age_bands,
                             age_bands_select = FALSE, vacc_uptake = TRUE,
                             by_class = FALSE, slice = FALSE)

  hfr_age <- get_eff_sev_age(dat, r, "hfr", vam_data, age_bands,
                             age_bands_select = FALSE, vacc_uptake = TRUE,
                             by_class = FALSE, slice = FALSE)

  ifr_age <- get_eff_sev_age(dat, r, "ifr", vam_data, age_bands,
                             age_bands_select = FALSE, vacc_uptake = TRUE,
                             by_class = FALSE, slice = FALSE)
  
  (ihr_age / hfr_age / ifr_age) +
    plot_annotation(tag_levels = "A")
}


suppl_outcomes_vacc_status <- function(dat, data, what) {
  
  model <- NULL
  for (i in sircovid::regions("england")) {
    tmp <- get_model_trajectory(dat, i, what, TRUE)
    tmp <- data.table::rbindlist(tmp, fill = TRUE, idcol = "name") %>%
      mutate(region = i)
    model <- rbind(model, tmp)
  }
  
  title <- ifelse(what == "diagnoses_vacc",
                  "Daily hospital admissions by vaccination status",
                  "Daily hospital deaths by vaccination status")
  
  dat <- if (what == "diagnoses_vacc") {data$admissions} else {data$deaths}
  data <- data %>%
    select(name, date, region) %>%
    mutate(date = as.Date(date),
           data = dat)
  
  df <- left_join(model, data, by = c("date", "region", "name"))
  df$name <- factor(df$name, levels = c("Unvaccinated", "Dose 1",
                                        "Dose 2", "Booster"))
  
  p <- ggplot(df, aes(x = date)) +
    geom_point(aes(y = data, col = name), size = 0.7, alpha = 0.5) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = name),
                alpha = 0.5, show.legend = FALSE) +
    labs(title = title, y = "", x = "") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y",
                 limits = c(as.Date("2020-12-08"), as.Date("2022-02-24")),
                 minor_breaks = NULL) +
    facet_grid(rows = vars(name), cols = vars(region), scales = "free_y")
  
  apply_theme(p) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, size = 11, vjust = 0.7))
  
}


plot_traceplots <- function(samples) {
  
  n_full_pars <- nrow(samples$full_pars)
  n_chains <- max(samples$chain)
  cols <- rev(viridisLite::viridis(n_chains))
  samples$chain_full <- rep(seq_len(n_chains), each = n_full_pars / n_chains)
  
  stopifnot(
    identical(samples$chain_full,
              rep(seq_len(n_chains),
                  each = length(samples$chain_full) / n_chains)))

  i <- spimalot:::reorder_beta(colnames(samples$full_pars))
  pars <- samples$full_pars[, i]
  
  nms <- colnames(pars)
  probs <- samples$full_probabilities
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  new_grid <- function(n, title) {
    par(mfrow = rep(ceiling(sqrt(n + 1)), 2),
        mar = c(3, 3, 2, 1),
        mgp = c(2, 0.5, 0),
        oma = c(1, 1, 1 + as.integer(title), 1))
  }
  
  plot_traces1 <- function(p, name) {
    traces <- matrix(p, ncol = n_chains)
    ess <- coda::effectiveSize(coda::as.mcmc(traces))
    
    matplot(traces, type = "l", lty = 1,
            xlab = "Iteration", bty = "n",
            ylab = name, col = cols,
            main = "",
            font.main = 1)
  }
  
  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$full_pars[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
  
}

