
all_variant_colours <- c("#723e8a", "#1ea1eb", "#f76002", "#eb2f97")
fit_cols <- c("green4", "red3")



gg_trajectory <- function(traj, ymax, ylab, shape = 1, annotate = FALSE) {
  
  ld_dates <- as.Date(c("2020-03-25", "2020-05-11",
                        "2020-11-05", "2020-12-02",
                        "2021-01-05", "2021-03-08",
                        "2021-03-08", "2021-04-19", "2021-05-17", "2021-07-19",
                        "2021-12-12", "2022-01-27"))
  
  p <- ggplot(traj, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(alpha = 0.3, aes(ymin = lb, ymax = ub, fill = "Model"), show.legend = FALSE) +
    geom_point(aes(y = data, col = "Data"), size = 1, shape = shape) +
    scale_color_manual(values = fit_cols,
                       guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                                size = c(2, 1),
                                                                shape = c(shape, NA)))) +
    scale_fill_manual(values = fit_cols[2]) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax * 1.35)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    labs(y = ylab, x = "")
  
  p <- apply_theme(p) +
    theme(axis.text.x = element_blank())
  
  if (annotate) {
    annotate_npis(p, ymax * 1.35)
  } else {
    p
  }
}


gg_trajectory_age <- function(traj, ylab, vax_ro, mute_axis_text = FALSE,
                              option = 'D', legend.tit = NULL) {
  
  p <- ggplot(traj, aes(date, value, fill = name, col = name)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_vline(xintercept = vax_ro, lty = 2) +
    labs(x = "", y = ylab) +
    scale_fill_viridis_d(name = legend.tit, direction = -1, alpha = 0.5, option = option) +
    scale_color_viridis_d(name = legend.tit, direction = -1, alpha = 0.5, option = option) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    guides(fill = guide_legend(ncol = 2))
  
  p <- apply_theme(p) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10),
          legend.title = element_text(size = 12))
  
  if (mute_axis_text) {
    p <- p +
      theme(axis.text.x = element_blank())
  }
  
  p
}


plot_effective_immunity_by_strain <- function(dat, r, na_variants = NULL) {
  
  ## Get trajectories
  # Susceptible population by strain
  date <- data.frame(date = as.Date(dat$data[[r]]$full$date_string))
  
  protected_strain_1 <- get_model_trajectory(
    dat, r, "prop_protected_1", TRUE)

  protected_strain_2 <-get_model_trajectory(
    dat, r, "prop_protected_2", TRUE)

  # Proportion of population not in S compartments
  S_vaccinated <- 
    get_model_trajectory(dat, r, "prop_protected_S_vaccinated", TRUE)
  
  # R that have been vaccinated
  R_vaccinated <- 
    get_model_trajectory(dat, r, "prop_protected_R_vaccinated", TRUE)
  
  # R that are unvaccinated
  R_unvaccinated <-
    get_model_trajectory(dat, r, "prop_protected_R_unvaccinated", TRUE)
  
  ## Process trajectories for plotting
  strain_epochs <- as.Date(c("2019-12-31", "2020-09-18", "2021-03-09", "2021-11-02"))
  
  wildtype <- 
    protected_strain_1[protected_strain_1$date < strain_epochs[3], ] %>%
    right_join(., date) %>%
    arrange(date)
  
  
  alpha <- rbind(
    protected_strain_2[protected_strain_2$date >= strain_epochs[2] &
                         protected_strain_2$date < strain_epochs[3], ],
    protected_strain_1[protected_strain_1$date >= strain_epochs[3] &
                         protected_strain_1$date < strain_epochs[4], ]) %>%
    right_join(., date) %>%
    arrange(date)
  
  delta <- rbind(
    protected_strain_2[protected_strain_2$date >= strain_epochs[3] &
                         protected_strain_2$date < strain_epochs[4], ],
    protected_strain_1[protected_strain_1$date >= strain_epochs[4], ]) %>%
    right_join(., date) %>%
    arrange(date)
  
  omicron <-
    protected_strain_2[protected_strain_2$date >= strain_epochs[4], ] %>%
    right_join(., date) %>%
    arrange(date)

  
  if(!is.null(na_variants)) {

    wildtype[wildtype$date >= na_variants[1], c(2:4)] <- NA_integer_
    alpha[alpha$date >= na_variants[2], c(2:4)] <- NA_integer_
    delta[delta$date >= na_variants[3], c(2:4)] <- NA_integer_
  }
  
  # Calculate prop protected by vaccine as the difference between
  # S and variant-specific immunity
  w_date <- max(wildtype$date[!is.na(wildtype$mean)])
  a_date <- max(wildtype$date[!is.na(alpha$mean)])
  d_date <- max(wildtype$date[!is.na(delta$mean)])
  diff_dates <- c(min(S_vaccinated$date),
                  w_date, a_date, d_date,
                  max(S_vaccinated$date))
  
  S_vaccinated$vaccine <- c(
    wildtype$mean[wildtype$date >= diff_dates[1] &
                    wildtype$date < diff_dates[2]],
    alpha$mean[alpha$date >= diff_dates[2] &
                 alpha$date < diff_dates[3]],
    delta$mean[delta$date >= diff_dates[3] &
                 delta$date < diff_dates[4]],
    omicron$mean[omicron$date >= diff_dates[4] &
                   omicron$date <= diff_dates[5]]
    )
  
  # Plot
  strain_names <- c("Wildtype", "Alpha", "Delta", "Omicron")
  
  p <- ggplot(S_vaccinated, aes(date, seq(0, 1, length.out = nrow(S_vaccinated))))
  
  p <- apply_theme(p)
  
  p +
    geom_vline(xintercept = as.Date("2020-12-08"), linetype = 2) +
    
    ggtext::geom_richtext(x = as.Date("2020-12-08"), y = 0.625,  size = 4.4,
                          angle = 90, label = "Vaccination rollout") +
    
    # Prop WITH nat. immunity
    geom_ribbon(alpha = 0.3,
                aes(ymin = 0, ymax = R_unvaccinated$mean,
                    fill = "% Inf. immunity")) +
    
    # Prop with both vacc and natural immunity
    geom_ribbon(alpha = 0.3,
                aes(ymin = R_unvaccinated$mean,
                    ymax = R_unvaccinated$mean + R_vaccinated$mean,
                    fill = "% Hybrid (vacc. + inf.)")) +
    
    # Prop with vacc. immunity
    geom_ribbon(alpha = 0.3,
                aes(ymin = R_unvaccinated$mean + R_vaccinated$mean,
                    ymax = mean + R_unvaccinated$mean + R_vaccinated$mean,
                    fill = "% Vacc. immunity")) +
    
    # Wildtype
    geom_line(aes(y = wildtype$mean, col = strain_names[1]),
              linewidth = 1, show.legend = FALSE) +
    
    # Alpha
    geom_line(aes(y = alpha$mean, col = strain_names[2]),
              linewidth = 1, show.legend = FALSE) +
    
    # Delta
    geom_line(aes(y = delta$mean, col = strain_names[3]),
              linewidth = 1, show.legend = FALSE) +
    
    # Omicron
    geom_line(aes(y = omicron$mean, col = strain_names[4]),
              linewidth = 1, show.legend = FALSE) +
    
    scale_color_manual(values = all_variant_colours,
                       breaks = strain_names) +
    
    scale_fill_manual(values = c("green4", "cyan3", "orange4"),
                      breaks = c("% Vacc. immunity",
                                 "% Hybrid (vacc. + inf.)",
                                 "% Inf. immunity"),
                      guide = guide_legend(override.aes = list(alpha = rep(0.2, 3)))) +
    
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0), limits = c(0, 1)) +
    labs(x = "", y = "Immunity against infection (%)") +
    theme(axis.text.x = element_blank())
  
}


get_eff_sev_age <- function(dat, r, w, vam_data, age_bands, age_bands_select,
                            vacc_uptake, by_class = FALSE, plot = TRUE,
                            slice = FALSE, slice_dates = NULL, slice_only = FALSE) {
  
  national <- get_national_intrinsic_values(dat, "Emergence3")
  
  tmp1 <- get_model_severity(dat$severity, vam_data, strain_epochs, w,
                            by_age = TRUE) %>%
    pivot_longer(!c(date, region, dominant_voc), names_to = "age") %>%
    mutate(estimate = case_when(stringr::str_detect(age, "mean") ~ "mean",
                                stringr::str_detect(age, "lb") ~ "lb",
                                stringr::str_detect(age, "ub") ~ "ub")) %>%
    rowwise() %>%
    mutate(age = gsub(paste0("_", estimate), "", age)) %>%
    mutate(age = gsub(paste0(w, "_age_"), "", age)) %>%
    mutate(vacc_class = "Overall")
  
  if (by_class) {
    
    tmp_vacc_class <- get_model_severity(dat$severity, vam_data, strain_epochs, w, # replace here
                                         by_vacc_class = TRUE,
                                         classes = c("unvacc", "1_no_eff",
                                                     "1_full_eff", "2_dose",
                                                     "2_dose_wane", "boost",
                                                     "boost_wane")) %>%
      pivot_longer(!c(date, region, dominant_voc), names_to = "age") %>%
      mutate(estimate = case_when(stringr::str_detect(age, "mean") ~ "mean",
                                  stringr::str_detect(age, "lb") ~ "lb",
                                  stringr::str_detect(age, "ub") ~ "ub")) %>%
      mutate(age = gsub("_mean|_lb|_ub", "", age))
    
    
    # Temporary filter out unvaccinated to rename vacc classes
    filter_out <- paste0(w, "_disag_", seq(0, 80, 5))
    
    unvacc <- tmp_vacc_class %>%
      filter(age %in% filter_out) %>%
      mutate(vacc_class = "Unvaccinated")
    
    tmp_vacc_class <- tmp_vacc_class %>%
      filter(!age %in% filter_out) %>%
      mutate(vacc_class = case_when(stringr::str_ends(age, "_1") ~ "1 dose no eff",
                                    stringr::str_ends(age, "_2") ~ "1 dose eff",
                                    stringr::str_ends(age, "_3") ~ "2 doses",
                                    stringr::str_ends(age, "_4") ~ "2 doses waned",
                                    stringr::str_ends(age, "_5") ~ "Booster",
                                    stringr::str_ends(age, "_6") ~ "Booster waned",
                                    TRUE ~ "Unvaccinated"))
    
    tmp_vacc_class <- rbind(tmp_vacc_class, unvacc) %>%
      mutate(age = gsub(paste0(w, "_disag_"), "", age),
             age = gsub("_1|_2|_3|_4|_5|_6", "", age))
    
    tmp1 <- rbind(tmp1, tmp_vacc_class)
    
    tmp1$vacc_class <- factor(tmp1$vacc_class,
                             levels = c(
                               "Overall", "Unvaccinated",
                               "1 dose no eff", "1 dose eff", "2 doses",
                               "2 doses waned", "Booster", "Booster waned"))
    
    tmp1$value[tmp1$date < as.Date("2020-12-08") &
                tmp1$vacc_class != "Overall"] <- NA_real_
    
  }
  
  tmp <- tmp1 %>%
    ungroup() %>%
    filter(region == r) %>%
    pivot_wider(names_from = estimate, values_from = value)
  
  tmp$age <- factor(tmp$age,
                    levels = unique(tmp$age),
                    labels = age_bands,
                    ordered = is.ordered(tmp$age))
  
  scaleFUN <- function(x) paste0(sprintf("%.1f", x * 100), "%")
  
  # Create colour scales and breaks
  if (w == "ihr") {
    scale <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.35, 0.4, 1)
    log_breaks <- sort(c(scale[-which(scale == 0.4)], 0.001, 0.0025, 0.005))
    cols <- c("cyan3", "chartreuse4", "yellow2", "orange1", "red4")
    
  } else if (w == "hfr") {
    scale <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1)
    log_breaks <- sort(c(0.015, 0.025, scale))
    cols <- c("cyan1", "blue3", "darkviolet", "grey20")
    
  } else if (w == "ifr") {
    scale <- c(0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.3, 1)
    log_breaks <- scale
    cols <- RColorBrewer::brewer.pal(length(scale), "Spectral")
  }
  
  if (slice) {
    
    stopifnot(!is.null(slice_dates))
    xintrcpt <- get_model_severity(dat$severity, vam_data, strain_epochs, w) %>%
      filter(region == "england", date == as.Date("2020-12-08")) %>%
      select(grep("mean", colnames(.), value = TRUE))
    
    xintrcpt <- as.numeric(unname(xintrcpt))
    
    
    slice_df <- tmp[tmp$date %in% slice_dates, ] %>%
      mutate(epoch = factor(case_when(date == slice_dates[1] ~ "Rollout",
                                      date == slice_dates[2] ~ "1st dose 50%",
                                      date == slice_dates[3] ~ "2nd dose 50%",
                                      date == slice_dates[4] ~ "Booster 50%"),
                            levels = c("Rollout", "1st dose 50%",
                                       "2nd dose 50%", "Booster 50%")))
    if (age_bands_select) {
      slice_df <- slice_df %>%
        filter(age %in% age_bands[9:17])
    }
    
    p_slice <- ggplot(slice_df) +
      geom_pointrange(shape = 18, linetype = 1, size = 0.8,
                      aes(x = mean, xmin = lb, xmax = ub, y = age, col = epoch)) +
      geom_vline(xintercept = xintrcpt, linetype = 2, linewidth = 0.7,
                 col = "red4") +
      scale_x_continuous(labels = scaleFUN,
                         breaks = log_breaks,
                         trans = "log") +
      labs(x = paste0(toupper(w), " (% log scale) *"),
           y = "Age") +
      scale_color_manual(name = "Milestone",
                         values = c("#BBBBBB", "#332288", "#DDCC77", "#117733"))
    
    p_slice <- apply_theme(p_slice) +
      theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.7),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.position = "top",
            legend.title = element_blank())
    
    # if (w != "ifr") {
    #   p_slice <- p_slice + xlab("")
    # }
    
    if (w != "ihr") {
      p_slice <- p_slice + theme(legend.position = "none")
    }
    
    if (w == "ifr") {
      p_slice <- p_slice +
        labs(caption =
               "
             Red dashed lines: date of vaccination rollout in effective severity plots (A-C), and national effective severity (weighted) at vaccination rollout in severity by age plots (D-F).
             
                *Effective severity by age at the time vaccination dose coverage was 50% of the entire national population." ) +
        theme(plot.caption = element_text(size = 10))
    }
  }
  
  
  if (age_bands_select) {
    age_bands <- age_bands[6:17]
  }
  
  p <- tmp %>%
    select(date, age, mean, vacc_class) %>%
    filter(age %in% c(age_bands)) %>%
    ggplot(aes(date, age)) +
    geom_tile(aes(fill = mean)) +
    desplot::geom_tileborder(aes(group = 1, grp = age), col = "white") +
    labs(x = "", y = "Age group", title = paste("Effective", toupper(w))) +
    scale_fill_gradientn(labels = scaleFUN,
                         na.value = "transparent",
                         colours = cols,
                         values = scale,
                         name = "",
                         breaks = scale) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.ticks = element_line(),
          legend.position = "right",
          legend.key.height = unit(2, "cm"),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 18))
  
  if (by_class) {
    p <- p +
      facet_wrap(~vacc_class)
  }
  
  if (plot) {
    if (slice) {
      
      if (slice_only) {
        p_slice
      } else {
        p <- p +
          theme(axis.text.x = element_blank())
        
        p + p_slice + plot_layout(widths = c(2.5, 1))
      }
      
    } else {
      p
    }
  } else {
    tmp1
  }
}


## TODO: check if this is used
plot_traj_age_vacc_class <- function(dat, vacc_class, voc_dates, r, what) {
  
  vacc_class_names <- c(
    "Unvaccinated", "First dose - no effect",
    "First dose - full effect", "Second dose",
    "Second dose - waned", "Booster", "Booster - waned")
  
  vacc_class <- vacc_class[[what]]
  vacc_class$value <- vacc_class[, what]
  
  sample <- dat$samples[[r]]
  
  state <- sample$trajectories$state
  state <- state[, , -1L]
  
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]
  
  vacc_class_number <- c("", paste0("_", seq_len(length(vacc_class_names) - 1)))
  
  if (what == "admissions") {
    traj_names <- grep("diagnoses_admitted_vacc", rownames(state), value = TRUE)
  } else if (what == "deaths") {
    traj_names <- grep("D_hosp_vacc", rownames(state), value = TRUE)
  }
  
  state <- state[traj_names, , ]
  state <- aperm(apply(state, c(1, 2), diff), c(2, 3, 1))
  x <- x[-1L]
  
  get_quantile <- function(dates_vect, p, bound) {
    q <- t(apply(state, c(1, 3), quantile, p))
    y <- data.frame(date = x, q)
      
    colnames(y) <- c("date", vacc_class_names)
    
    y$`First dose` <- y$`First dose - no effect` +
      y$`First dose - full effect`
    y$`Second dose` <- y$`Second dose` + y$`Second dose - waned`
    y$`Booster` <- y$`Booster` + y$`Booster - waned`
    
    y %>% 
      select(date, Unvaccinated, `First dose`, `Second dose`, Booster) %>%
      melt(., id = "date") %>%
      group_by(date, variable) %>%
      filter(date >= as.Date("2020-12-08")) %>%
      summarise(n = sum(value)) %>%
      mutate(percentage = n / sum(n)) %>%
      mutate(bound = bound) %>%
      ungroup() %>%
      select(date, name = variable, bound, percentage)
  }
  
  y <- rbind(get_quantile(y, 0.0275, "lower"),
             get_quantile(y, 0.975, "upper"),
             get_quantile(y, 0.5, "mean"))
  
  y <- y %>%
    pivot_wider(names_from = bound, values_from = percentage) %>%
    rowwise() %>%
    mutate(ymax = max(upper, lower),
           ymin = min(upper, lower)) %>%
    mutate(across(where(is.numeric), ~na_if(., 0)))
  
  vacc_classes_keep <- c("Unvaccinated", "First dose", "Second dose", "Booster")
  
  y$source <- factor(y$name, levels = vacc_classes_keep)
  
  tmp_data <- vacc_class %>%
    filter(region == r) %>%
    select(date, age_group, variable = vacc_status, value) %>%
    filter(!is.na(variable)) %>%
    mutate(variable = case_when(
      variable == 0 ~ "Unvaccinated",
      variable == 1 ~ "First dose",
      variable == 2 ~ "Second dose",
      variable == 3 ~ "Booster"
    )) %>%
    group_by(date, variable) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>%
    mutate(Total = rowSums(.[grep("date", names(.), invert = TRUE)], na.rm = TRUE)) %>%
    group_by(year = year(date), week = week(date)) %>% 
    summarise_if(is.numeric, sum) %>%
    select(!c(year, week))
  
  tmp_data$date <- seq.Date(min(vacc_class$date), max(vacc_class$date + 7), by = "week")
  
  vacc_names_plot <- vacc_classes_keep
  data <- NULL
  for (i in 1:nrow(tmp_data)) {
    
    d <- tmp_data[i, "date"]
    
    n_vect <- as.numeric(as.vector(tmp_data[i, vacc_names_plot]))
    cis <- DescTools::MultinomCI(n_vect, method = "wald")
    cis <- data.frame(date = d, source = vacc_names_plot, cis)
    
    colnames(cis) <- c("date", "source", "mean", "lb", "ub")
    data <- rbind(data, cis)
    
    data <- data %>%
      naniar::replace_with_na(replace = list(mean = c(0, 1),
                                             lb = c(0, 1),
                                             ub = c(0, 1)))
  }
  
  data$source <- factor(data$source, levels = vacc_names_plot)
  data[data$source == "Booster" & data$date < as.Date("2021-09-14"),
       c("mean", "lb", "ub")] <- NA_real_
  
  breaks <- unique(data$source)
  cols <- khroma::colour("muted")(length(breaks))
  names(cols) <- breaks
  
  p <- ggplot() +
    geom_point(aes(x = data$date, y = data$mean, col = data$source), size = 0.8) +
    geom_line(aes(x = as.Date(y$date), y = y$mean, col = y$source), linewidth = 0.8) +
    geom_errorbar(aes(x = data$date, ymin = data$lb, ymax = data$ub,
                      col = data$source), alpha = 0.3, width = 1,
                  position = position_dodge(0.1)) +
    labs(title = paste("Hospital", what), y = "% daily total", x = "") +
    scale_color_manual(values = cols) +
    scale_x_date(limits = c(min(x), max(x)), date_breaks = "2 months",
                 date_labels = "%b %y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.line = element_line(),
          axis.text.x = element_blank(),
          legend.text = element_text(size = 14),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12))
  
  annotate_strains(p, voc_dates = voc_dates, ymax = 1, TRUE)
}


plot_eff_severity_by_strain <- function(dat, r, what = c("ihr", "hfr", "ifr"),
                                        strain_epochs, ymax_values = NULL,
                                        na_variants = NULL, national_box = FALSE,
                                        national = NULL) {
  

  # Arrange national intrinsic severity for plotting
  national <- switch_levels(national)
  
  # Get VOC frequency data we fitted to (VAM)
  vam_data <- get_strain_timelines(dat$data)
  
  tmp <- lapply(what, function(w)
    get_model_severity(dat$severity, vam_data, strain_epochs, w, by_strain = TRUE))
  names(tmp) <- toupper(what)
  
  # Some labels for plotting later
  strain_labels <- data.frame(
    date = strain_epochs,
    label = c("Alpha emerges", "Delta emerges", "Omicron emerges"))
  
  out <- NULL
  for (w in names(tmp)) {
    
    ret <- tmp[[w]] %>%
      filter(region == r) %>%
      pivot_longer(!c(date, region, dominant_voc)) %>%
      mutate(strain = case_when(
        grepl("strain_1", name) ~ "Strain 1",
        grepl("strain_2", name) ~ "Strain 2",
        TRUE ~ "Weighted")) %>%
      mutate(name = case_when(
        grepl("mean", name) ~ "mean",
        grepl("lb", name) ~ "lb",
        grepl("ub", name) ~ "ub")) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      mutate(strain = case_when(
        strain == "Strain 1" ~ case_when(
          date < strain_epochs[1] ~ NA_character_,
          date >= strain_epochs[1] & date < strain_epochs[2] ~ "Wildtype",
          date >= strain_epochs[2] & date < strain_epochs[3] ~ "Alpha",
          date >= strain_epochs[3] ~ "Delta"),
        strain == "Strain 2" ~ case_when(
          date < strain_epochs[1] ~ NA_character_,
          date >= strain_epochs[1] & date < strain_epochs[2] ~ "Alpha",
          date >= strain_epochs[2] & date < strain_epochs[3] ~ "Delta",
          date >= strain_epochs[3] ~ "Omicron"),
        TRUE ~ "Weighted"
      )) %>%
      filter(!is.na(strain))
    
    ret$strain[ret$date < strain_epochs[1]] <- "Wildtype"
    ret$strain <- factor(ret$strain,
                              levels = c("Weighted", "Wildtype", "Alpha", "Delta", "Omicron"))
    
    # If dates to NA variants provided, remove from trajectories to plot
    if (!is.null(na_variants)) {
      
      na_wt_dates <- ret$strain == "Wildtype" & ret$date >= na_variants[1] & ret$date <= na_variants[2]
      ret[na_wt_dates, c("mean", "lb", "ub")] <- NA_integer_
      
      na_al_dates <- ret$strain == "Alpha" & ret$date >= na_variants[2] & ret$date <= na_variants[3]
      ret[na_al_dates, c("mean", "lb", "ub")] <- NA_integer_
       
      na_dl_dates <- ret$strain == "Delta" & ret$date >= na_variants[3]
      ret[na_dl_dates, c("mean", "lb", "ub")] <- NA_integer_
      
    }
    
    ymax <- ifelse(!is.null(ymax_values), ymax_values[[w]], max(ret$ub, na.rm = TRUE) * 1.1)
    
    ntl <- NULL
    if (national_box) {
      ntl <- national[[w]]
      ymax <- max(ymax, max(unlist(ntl)))
      df <- data.frame(matrix(unlist(ntl), nrow = length(ntl), byrow = TRUE))
      colnames(df) <- names(ntl[[1]])
      df <- df %>%
        mutate(variant = factor(names(ntl), levels = names(ntl)))
      
      ntl <- ggplot(df, aes(x = NA, y = mean, col = variant, fill = variant)) +
        geom_point(shape = 18, size = 2, alpha = 0.7, position = position_dodge(width = 0.5)) +
        geom_crossbar(aes(ymin = lb, ymax = ub), position = "dodge", width = 0.5, alpha = 0.15) +
        scale_color_manual(values = all_variant_colours) +
        scale_fill_manual(values = all_variant_colours) +
        scale_y_continuous(labels = scales::percent_format(accuracy = ifelse(w == "HFR", 1, 0.1)),
                           expand = c(0, 0),
                           limits = c(0, ymax)) +
        scale_x_discrete(breaks = NULL) +
        labs(x = "", y = paste("Basic ", w))
      
      ntl <- apply_theme(ntl) + theme(legend.position = "none",
                                      axis.text.x = element_blank())
    }
    
    p <- ret %>%
      ggplot(., aes(date)) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = strain), alpha = 0.2) +
      geom_line(aes(y = mean, col = strain)) +
      geom_vline(xintercept = as.Date("2020-12-08"), linetype = 2, linewidth = 0.7,
                 col = "red4") +
      scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
      scale_y_continuous(labels = scales::percent_format(accuracy = ifelse(w == "HFR", 1, 0.1)),
                         expand = c(0, 0),
                         limits = c(0, ymax)) +
      scale_color_manual(values = c("#000000", all_variant_colours)) +
      scale_fill_manual(values = c("#000000", all_variant_colours)) +
      labs(x = "", y = paste("Effective", w))
    
    p <- apply_theme(p) + theme(legend.position = "none")
    
    if (w %in% c("IHR", "HFR")) {
       p <- p +
         theme(axis.text.x = element_text(colour = "white", size = 10, vjust = 0.7, angle = 45))
    } else if (w == "IFR") {
      p <- p +
        theme(axis.text.x = element_text(size = 10, vjust = 0.7, angle = 45),
              plot.caption = element_text(size = 10))
    }
    
    if (w == "IHR") {
      p <- p + theme(legend.position = "top")
    }
    
    if (!is.null(ntl)) {
      
      p <- p + 
        scale_y_continuous(labels = scales::percent_format(accuracy = ifelse(w == "HFR", 1, 0.1)),
                           expand = c(0, 0),
                           limits = c(0, ymax),
                           position = "right")
      p <- ntl + p + plot_layout(widths = c(0.25, 1))
    }
    
    out[[w]] <- p
      
  }
  
  out[[1]] / out[[2]] / out[[3]]
  
}


paper_plot_2 <- function(dat, r, strain_epochs, vam_data,
                         cut_off = 0.9, age_bands_select = FALSE) {
  
  # Get cut-off date to NA replaced variants (i.e. new variant frequency >=90%)
  vam_data <- vam_data %>%
    filter(region == r)
  wildtype_date <- min(vam_data$date[vam_data$Alpha >= cut_off], na.rm = TRUE)
  alpha_date <- min(vam_data$date[vam_data$Delta >= cut_off], na.rm = TRUE)
  delta_date <- min(vam_data$date[vam_data$Omicron >= cut_off], na.rm = TRUE)
  na_variants <- c(wildtype_date, alpha_date, delta_date)
  
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
  
  # Get national intrinsic severity
  national <- get_national_intrinsic_values(dat, "Emergence3")

  # Effective severity by age
  age_bands <- get_age_band_labels(TRUE)
  ihr_age <- get_eff_sev_age(dat, r, "ihr", vam_data, age_bands,
                             age_bands_select = age_bands_select, vacc_uptake = TRUE,
                             by_class = FALSE, slice = TRUE,
                             slice_dates = c(vax_ro, slice_dates),
                             slice_only = TRUE)

  hfr_age <- get_eff_sev_age(dat, r, "hfr", vam_data, age_bands,
                             age_bands_select = age_bands_select, vacc_uptake = TRUE,
                             by_class = FALSE, slice = TRUE,
                             slice_dates = c(vax_ro, slice_dates),
                             slice_only = TRUE)

  ifr_age <- get_eff_sev_age(dat, r, "ifr", vam_data, age_bands,
                             age_bands_select = age_bands_select, vacc_uptake = TRUE,
                             by_class = FALSE, slice = TRUE,
                             slice_dates = c(vax_ro, slice_dates),
                             slice_only = TRUE)
  
  sev_age <- (ihr_age / hfr_age / ifr_age)
  
  # Plot effective severity
  sev <- plot_eff_severity_by_strain(dat, r, strain_epochs = strain_epochs,
                                     na_variants = na_variants, national_box = TRUE,
                                     national = national, ymax_values = c(IHR = 0.05,
                                                                          HFR = 0.5,
                                                                          IFR = 0.03))
  
  (sev | sev_age) +
    plot_layout(widths = c(1, 0.55)) +
    plot_annotation(tag_levels = list(c('A', '', 'B', '', 'C', '', 'D', 'E', 'F')))
  
}


paper_plot_1 <- function(dat, r, vam_data) {
  
  date <- sircovid::sircovid_date_as_date(dat$samples[[r]]$trajectories$date)[-1L]
  
  # R0 by variant
  variants <- c("Wildtype", "Alpha", "Delta", "Omicron")
  r0 <- get_national_intrinsic_values(dat, "Emergence3")
  
  r0 <- unlist(switch_levels(r0)$r0)
  r0 <- matrix(r0, nrow = length(r0), byrow = TRUE) %>%
    as.data.frame(.) %>%
    magrittr::set_colnames(., "value") %>%
    mutate(name = rep(c("mean", "lb", "ub"), times = 4),
           variant = factor(rep(variants, each = 3),
                            levels = variants)) %>%
    pivot_wider()
  
  ymax_r <- ceiling(max(r0$ub))
  
  r0 <- ggplot(r0, aes(x = NA, y = mean, col = variant, fill = variant, group = variant)) +
    geom_point(size = 3, shape = 18, position = position_dodge(width = 0.5)) +
    geom_crossbar(aes(ymin = lb, ymax = ub),
                  width = 0.5, alpha = 0.15, position = "dodge") +
    labs(x = "", y = "Intrinsic R0") +
    scale_color_manual(values = all_variant_colours) +
    scale_fill_manual(values = all_variant_colours) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, ymax_r),
                       breaks = seq(0, ymax_r, 1))
  
  r0 <- apply_theme(r0) + theme(axis.text.y = element_text(size = 12),
                                axis.text.x = element_blank(),
                                legend.position = "top",
                                legend.text = element_text(size = 12))
  
  # Get cut-off date to NA replaced variants (i.e. new variant frequency >=90%)
  vam_data <- vam_data %>% filter(region == r)
  wildtype_date <- min(vam_data$date[vam_data$Alpha >= 0.9], na.rm = TRUE)
  alpha_date <- min(vam_data$date[vam_data$Delta >= 0.9], na.rm = TRUE)
  delta_date <- min(vam_data$date[vam_data$Omicron >= 0.9], na.rm = TRUE)
  na_variants <- c(wildtype_date, alpha_date, delta_date)
  
  # Get model Rt for plot 1
  mod_rt <- list(dat$rt[[r]]$eff_Rt_general[-1L, "weighted", ],
                 dat$rt[[r]]$Rt_general[-1L, "weighted", ])
  eff_rt <- data.frame(date,
                       mean = rowMeans(mod_rt[[1]]),
                       lb = matrixStats::rowQuantiles(mod_rt[[1]], probs = 0.025),
                       ub = matrixStats::rowQuantiles(mod_rt[[1]], probs = 0.975))
  rt <- data.frame(date,
                   mean = rowMeans(mod_rt[[2]]),
                   lb = matrixStats::rowQuantiles(mod_rt[[2]], probs = 0.025),
                   ub = matrixStats::rowQuantiles(mod_rt[[2]], probs = 0.975))
  
  p1 <- ggplot(eff_rt, aes(x = date)) +
    geom_line(aes(y = mean, col = "Effective Rt")) +
    geom_ribbon(show.legend = FALSE, aes(ymin = lb, ymax = ub, fill = "Effective Rt"), alpha = 0.3) +
    geom_line(aes(y = rt$mean, col = "Rt")) +
    geom_ribbon(show.legend = FALSE, aes(ymin = rt$lb, ymax = rt$ub, fill = "Rt"), alpha = 0.3) +
    geom_hline(yintercept = 1, lty = 2) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, ymax_r),
                       breaks = seq(0, ymax_r, 1)) +
    labs(x = "", y = "Reproduction number") +
    scale_color_manual(values = c("black", "red"), breaks = c("Rt", "Effective Rt")) +
    scale_fill_manual(values = c("black", "red"), breaks = c("Rt", "Effective Rt"))
  
  p1 <- annotate_npis(p1, ymax = ymax_r, TRUE,
                      ldwn_y = c(3, 3.5, 5, 3.5, 6, rep(7.5, 2)), vax_y = 3.3)
  p1 <- apply_theme(p1) +
    theme(axis.text.x = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 12))
  
  p2 <- plot_effective_immunity_by_strain(dat, r, na_variants) +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.7),
          legend.position = "bottom",
          legend.text = element_text(size = 12))
  
  
  # Fitted trajectories plots
  pillar2 <- ggplot_pillar_2(dat, r)
  
  react <- ggplot_infection_prevalence(dat, r, "react") +
    theme(axis.text.x = element_blank(),
          legend.position = "none")

  vacc_status <- plot_vacc_status(dat, r)
  
  r0 <- r0 / vacc_status
         
  
  inf_strain <- plot_inf_strain(dat, r, variants)
  
  timeseries <- pillar2 / react / inf_strain
  
  # Arrange all with patchwork
  (timeseries | r0 | (p1 / p2)) +
    plot_layout(widths = c(1, 0.4, 1)) +
    plot_annotation(tag_levels = 'A')
}


ggplot_infection_prevalence <- function(dat, r, which, shape = 1) {
  
  sample <- dat$samples[[r]]$trajectories$state
  pars <- dat$parameters$base[[1]]$sens_and_spec
  data <- dat$data[[r]]
  date <- dat$info$date
  min_tot <- 50
  
  pop <- sircovid:::sircovid_population(r)
  
  if (which == "ons") {
    
    N_tot <- pop[1] * 3 / 5 + sum(pop[2:17])
    pos <- sample["ons_pos", , -1]
    sens <- pars$ons_sensitivity
    spec <- pars$ons_specificity
    ylab <- "ONS proportion positive (%)"
    cols <- c("#E48C2A", fit_cols[2])
    
    # data
    npos <- data$full[, "ons_pos"]
    ntot <- data$full[, "ons_tot"]
    
  } else if (which == "react") {
    
    N_tot <- sum(pop[2:17])
    pos <- sample["react_pos", , -1]
    sens <- pars$react_sensitivity
    spec <- pars$react_specificity
    ylab <- "Infection prevalence (%)"
    cols <- fit_cols
    
    # data
    npos <- data$full[, "react_pos"]
    ntot <- data$full[, "react_tot"]
    
  }
  
  neg <- N_tot - pos
  
  res <- (pos * sens +
            neg * (1 - spec)) / (pos + neg) * 100
  
  res <- data.frame(
    date = as.Date(as.Date(data$fitted$date_string)),
    mean = colMeans(res),
    lb = matrixStats::colQuantiles(res, probs = 0.025),
    ub = matrixStats::colQuantiles(res, probs = 0.975))
  
  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0
  npos[ntot < min_tot] <- 0
  ntot[ntot < min_tot] <- 0
  
  cis <- Hmisc::binconf(x = npos, n = ntot) * 100
  dy <- cis[, "PointEst"]
  lower <- cis[, "Lower"]
  upper <- cis[, "Upper"]
  dy[ntot == 0] <- NA
  
  xmin <- min(which(!is.na(dy)))
  
  res[c(1:xmin), c("mean", "lb", "ub")] <- NA
  res$data_mean <- dy
  res$data_lb <- lower
  res$data_ub <- upper
  
  ymax <- max(res$ub, na.rm = TRUE) * 1.1
  
  p <- ggplot(res, aes(x = date)) +
    geom_line(aes(y = mean, col = "Model")) +
    geom_ribbon(alpha = 0.3, aes(ymin = lb, ymax = ub, fill = "Model"), show.legend = FALSE) +
    geom_point(aes(y = data_mean, col = "Data"), shape = shape) +
    geom_errorbar(aes(ymin = data_lb, ymax = data_ub, col = "Data")) +
    scale_color_manual(values = cols,
                       guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                                shape = c(shape, NA)))) +
    scale_fill_manual(values = cols[2]) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    labs(y = ylab, x = "")
  
  p <- apply_theme(p) +
    theme(axis.text.x = element_blank())
  
  
  annotate_npis(p, ymax, vax_y = ymax * 0.6, ldwn_y = rep(ymax * 0.95, 5), label = FALSE)
  
}


ggplot_pillar_2 <- function(dat, r) {

  sample <- dat$samples[[r]]
  pars <- dat$parameters$base[[1]]$sens_and_spec
  data <- dat$data[[r]]
  date <- dat$info$date
  
  xmin <- as.Date("2020-05-15")
  ylab <- "Pillar 2 over 15 positivity (%)"
  
  # data
  npos <- data$full[, "pillar2_over25_pos"] + data$full[, "pillar2_15_24_pos"]
  ntot <- data$full[, "pillar2_over25_tot"] + data$full[, "pillar2_15_24_tot"]
  npos[is.na(npos)] <- 0
  ntot[is.na(ntot)] <- 0
  
  # trajectory
  traj_pos <- sample$trajectories$state["sympt_cases_15_24_inc", , -1] +
    sample$trajectories$state["sympt_cases_over25_inc", , -1]
  traj_neg <- sample$trajectories$state["pillar2_negs_15_24", , -1] +
    sample$trajectories$state["pillar2_negs_over25", , -1]
  
  traj <-
    (traj_pos * pars$pillar2_sensitivity +
       traj_neg * (1 - pars$pillar2_specificity)) / (traj_pos + traj_neg) * 100
  
  
  traj_date <- 
    sircovid::sircovid_date_as_date(sample$trajectories$date[-1])
  
  res <- data.frame(
    date = traj_date,
    mean = colMeans(traj),
    lb = matrixStats::colQuantiles(traj, probs = 0.025),
    ub = matrixStats::colQuantiles(traj, probs = 0.975)
  )
  
  res[res$date < xmin, c("mean", "lb", "ub")] <- NA
  
  cis <- Hmisc::binconf(x = npos, n = ntot) * 100
  dy <- cis[, "PointEst"]
  dy[ntot == 0] <- NA
  
  res$data <- dy
  
  
  ymax <- max(res$ub, na.rm = TRUE) * 1.05
  res1 <- res %>%
    `colnames<-`(., c("date", "Model", "lb", "ub", "Data")) %>%
    pivot_longer(-c(date, lb, ub))
  
  p <- ggplot(res1, aes(x = date, y = value, col = name)) +
    geom_line() +
    geom_ribbon(alpha = 0.3, aes(ymin = lb, ymax = ub, fill = name)) +
    geom_point(aes(col = name), shape = NA) +
    scale_color_manual(values = fit_cols,
                       guide = guide_legend(override.aes = list(linetype = c(1, 1),
                                                                size = c(2, 1),
                                                                shape = c(1, NA)))) +
    scale_fill_manual(values = fit_cols) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    labs(y = ylab, x = "") 
  
  p <- apply_theme(p) +
    theme(axis.text.x = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 12))
  
  annotate_npis(p, ymax)
  
}


plot_intrinsic_severity_vs_transmissibility <- function(dat, r, what, which) {

  intrinsic_sev <- dat$intrinsic_severity %>%
    filter(region %in% sircovid::regions(r))
  
  # Get model R0 by region and range to plot
  r0_reg <- dplyr::bind_rows(
    lapply(sircovid::regions(r), function (r) get_r0_region(dat, r)))

  r0_plot <- r0_reg %>%
    filter(value == "mean") %>%
    select(!value)
  
  r0_range <- matrix(
    c(c(min(r0_plot$Wildtype), max(r0_plot$Wildtype)),
      c(min(r0_plot$Alpha), max(r0_plot$Alpha)),
      c(min(r0_plot$Delta), max(r0_plot$Delta)),
      c(min(r0_plot$Omicron), max(r0_plot$Omicron))),
    byrow = TRUE, ncol = 2)
  colnames(r0_range) <- c("min", "max")
  
  r0_reg <- r0_reg %>%
    mutate(estimate = value) %>%
    select(!value)

  out <- NULL
  # Get model intrinsic severity by region and range to plot
  for (w in what) {
    
    sev <- intrinsic_sev %>%
      pivot_wider(names_from = name, values_from = value) %>%
      filter(source == w, period == which)
    
    ret <- sev %>%
      mutate(source = "sev") %>%
      select(colnames(r0_reg)) %>%
      rbind(., r0_reg) %>%
      mutate(char = case_when(
        region == "east_of_england" ~ "EE",
        region == "london" ~ "LON",
        region == "midlands" ~ "MID",
        region == "north_east_and_yorkshire" ~ "NEY",
        region == "north_west" ~ "NW",
        region == "south_east" ~ "SE",
        region == "south_west" ~ "SW")) %>%
      pivot_longer(!c(region, source, char, estimate)) %>%
      mutate(estimate = paste0(estimate, "_", source)) %>%
      select(!source) %>%
      pivot_wider(names_from = estimate, values_from = value) %>%
      mutate(strain = factor(name, levels = c("Wildtype", "Alpha", "Delta", "Omicron")))
    
    sev_plot <- sev %>%
      filter(estimate == "mean") %>%
      select(Wildtype, Alpha, Delta, Omicron)
    
    sev_range <- matrix(c(c(min(sev_plot$Wildtype), max(sev_plot$Wildtype)),
                          c(min(sev_plot$Alpha), max(sev_plot$Alpha)),
                          c(min(sev_plot$Delta), max(sev_plot$Delta)),
                          c(min(sev_plot$Omicron), max(sev_plot$Omicron))),
                        byrow = TRUE, ncol = 2)
    colnames(sev_range) <- c("min", "max")
    
    
    ymax <- max(sev_range) * 1.1
    
    p <- ggplot(ret, aes(col = strain, mean_r0, mean_sev)) +
      geom_point() +
      geom_errorbar(aes(ymin = lb_sev, ymax = ub_sev)) +
      geom_errorbar(aes(xmin = lb_r0, xmax = ub_r0)) +
      annotate(geom = "rect", xmin = r0_range[1, "min"], xmax = r0_range[1, "max"],
               ymin = sev_range[1, "min"], ymax = sev_range[1, "max"],
               alpha = 0.1, fill = all_variant_colours[1], col = all_variant_colours[1], lty = 3) +
      
      annotate(geom = "rect", xmin = r0_range[2, "min"], xmax = r0_range[2, "max"],
               ymin = sev_range[2, "min"], ymax = sev_range[2, "max"],
               alpha = 0.1, fill = all_variant_colours[2], col = all_variant_colours[2], lty = 3) +
      
      annotate(geom = "rect", xmin = r0_range[3, "min"], xmax = r0_range[3, "max"],
               ymin = sev_range[3, "min"], ymax = sev_range[3, "max"],
               alpha = 0.1, fill = all_variant_colours[3], col = all_variant_colours[3], lty = 3) +

      annotate(geom = "rect", xmin = r0_range[4, "min"], xmax = r0_range[4, "max"],
               ymin = sev_range[4, "min"], ymax = sev_range[4, "max"],
               alpha = 0.1, fill = all_variant_colours[4], col = all_variant_colours[4], lty = 3) +
      
      ggrepel::geom_text_repel(aes(label = char), size = 4, show.legend = FALSE) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ymax),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_shape_manual(values = c(1:7)) +
      labs(x = "", y = paste("Intrinsic", w)) +
      scale_x_continuous(limits = c(0, ceiling(max(r0_range))), expand = c(0, 0),
                         breaks = seq(0, ceiling(max(r0_range)), 1),
                         minor_breaks = seq(0, ceiling(max(r0_range)), 0.5)) +
      scale_color_manual(values = all_variant_colours) +
      theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.text = element_text(size = 10))
    
    out[[w]] <- apply_theme(p) +
      theme(axis.text.x = element_blank(),
            legend.position = "none")
    
    if (length(what) == 3 && w != "HFR") {
      out[[w]] <- out[[w]] +
        theme(legend.position = "none")
    }
    
    out[["ymax_values"]][[w]] <- ymax
    
  }

  out

}


plot_intrinsic_severity <- function(dat, r, what, when) {
 
  intrinsic_sev <- dat$intrinsic_severity
  
  out <- NULL
  for (w in what) {
    
    sev <- intrinsic_sev %>%
      pivot_wider(names_from = name, values_from = value) %>%
      filter(period == when, source == w)
    
    ret <- sev %>%
      mutate(region = case_when(
        region == "east_of_england" ~ "EE",
        region == "london" ~ "LON",
        region == "midlands" ~ "MID",
        region == "north_east_and_yorkshire" ~ "NEY",
        region == "north_west" ~ "NW",
        region == "south_east" ~ "SE",
        region == "south_west" ~ "SW"
      )) %>%
      select(!c(source, period)) %>%
      pivot_longer(!c(estimate, region)) %>%
      mutate(name = factor(name, levels = c("Wildtype", "Alpha", "Delta", "Omicron"))) %>%
      pivot_wider(names_from = estimate)
    
    
    ymax <- max(ret$ub) * 1.1
    
    p <- ggplot(ret, aes(x = name, col = name, group = region)) +
      geom_point(aes(y = mean), position = position_dodge(width = 1),
                 size = 0.8, show.legend = FALSE) +
      geom_errorbar(aes(ymin = lb, ymax = ub), position = position_dodge(width = 1),
                    size = 0.3, show.legend = FALSE) +
      geom_text(aes(y = mean * 1.05, label = region), position = position_dodge(width = 1),
                show.legend = FALSE) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, ymax),
                         labels = scales::percent_format(accuracy = 1)) +
      labs(x = "", y = paste("Intrinsic", w)) +
      scale_color_manual(values = all_variant_colours)
    
    p <- apply_theme(p) +
      theme(axis.text.x = element_blank())
    
    out[[w]] <- p
    
    out[["ymax_values"]][[w]] <- ymax
  }
  
  out
}


annotate_npis <- function(p, ymax, label = FALSE,
                        ldwn_y = NULL, vax_y = NULL) {
  
  ld_dates <- as.Date(c("2020-03-25", "2020-05-11", # lockdown 1
                        "2020-05-11", "2020-06-01", "2020-06-15", "2020-07-04", # gradual easing of 1st ldwn
                        "2020-10-14", "2020-11-05", # three-tier system
                        "2020-11-05", "2020-12-02", # lockdown 2
                        "2020-12-02", "2020-12-21", # 3-tiers
                        "2020-12-21", "2021-01-05", # 4-tiers
                        "2021-01-05", "2021-03-08", # lockdown 3
                        "2021-03-08", "2021-04-19", "2021-05-17", "2021-07-19", # roadmap
                        "2021-12-12", "2022-01-27")) # plan B
  vax_ro <- as.Date("2020-12-08")
  
  p <- p + 
    # Vaccine rollout
    geom_vline(xintercept = vax_ro, lty = 2) +
    
    # NPIs
    annotate("rect", xmin = ld_dates[1], xmax = ld_dates[2],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.55) +
    annotate("rect", xmin = ld_dates[3], xmax = ld_dates[4],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.35) +
    annotate("rect", xmin = ld_dates[4], xmax = ld_dates[5],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.2) +
    annotate("rect", xmin = ld_dates[5], xmax = ld_dates[6],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.1) +
    annotate("rect", xmin = ld_dates[7], xmax = ld_dates[8],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.2) +
    annotate("rect", xmin = ld_dates[9], xmax = ld_dates[10],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.55) +
    annotate("rect", xmin = ld_dates[11], xmax = ld_dates[12],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.2) +
    annotate("rect", xmin = ld_dates[13], xmax = ld_dates[14],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.35) +
    annotate("rect", xmin = ld_dates[15], xmax = ld_dates[16],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.55) +
    annotate("rect", xmin = ld_dates[17], xmax = ld_dates[18],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.35) +
    annotate("rect", xmin = ld_dates[18], xmax = ld_dates[19],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.2) +
    annotate("rect", xmin = ld_dates[19], xmax = ld_dates[20],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.1) +
    annotate("rect", xmin = ld_dates[21], xmax = ld_dates[22],
             ymin = 0, ymax = ymax,
             fill = "grey50", alpha = 0.2)
  
  if (label) {
    
    p <- p +
      annotate("label", size = 4, y = c(7, 7, 7),
               x = ld_dates[c(4, 16, 21)] + c(10, 70, 23),
               label = c(
" Easing of
 Lockdown 1", "Roadmap", "Plan B")) +
      ggtext::geom_richtext(x = ld_dates[1] + 24, y = 3.5, size = 4,
                            angle = 90, label = "Lockdown 1") +
      ggtext::geom_richtext(x = ld_dates[9] + 12, y = 3.5, size = 4,
                            angle = 90, label = "Lockdown 2") +
      ggtext::geom_richtext(x = ld_dates[14] + 31, y = 3.5, size = 4,
                            angle = 90, label = "Lockdown 3") +
      ggtext::geom_richtext(x = ld_dates[7] + 10, y = 6, size = 4,
                            angle = 90, label = "3 Tiers") +
      ggtext::geom_richtext(x = ld_dates[10] + 21, y = 6, size = 4,
                            angle = 90, label = "4 Tiers")
  }
  
  p
}

apply_theme <- function(p) {
  
  p +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.title.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.text.y = element_text(size = 12))
}


annotate_strains <- function(p, ymax, vam_data = NULL, voc_dates = NULL,
                             r = NULL, coff = NULL,
                             label = FALSE) {
  
  if (is.null(voc_dates)) {
    
    vam_data <- vam_data[vam_data$region == r, ]
    
    # time horizon of the severity analysis
    min_max <- sircovid::sircovid_date_as_date(c(76, 786))
    
    # establish dates at which vocs where x% dominant
    voc_dates <- c(
      min(min_max[1]),
      min(vam_data$date[which(vam_data$Alpha >= coff[1])]),
      min(vam_data$date[which(vam_data$Alpha >= coff[2])]),
      min(vam_data$date[which(vam_data$Alpha >= coff[3])]),
      min(vam_data$date[which(vam_data$Delta >= coff[1])]),
      min(vam_data$date[which(vam_data$Delta >= coff[2])]),
      min(vam_data$date[which(vam_data$Delta >= coff[3])]),
      min(vam_data$date[which(vam_data$Omicron >= coff[1])]),
      min(vam_data$date[which(vam_data$Omicron >= coff[2])]),
      min(vam_data$date[which(vam_data$Omicron >= coff[3])]),
      max(min_max[2])
    )
  }
  
  p <- p +
    # Wildtype
    annotate("rect", xmin = voc_dates[1], xmax = voc_dates[2] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[1], alpha = 0.1) +
    
    # Alpha period
    annotate("rect", xmin = voc_dates[2], xmax = voc_dates[3] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[2], alpha = 0.05) +
    annotate("rect", xmin = voc_dates[3], xmax = voc_dates[4] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[2], alpha = 0.15) +
    annotate("rect", xmin = voc_dates[4], xmax = voc_dates[5] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[2], alpha = 0.25) +
    
    # Delta rectangles
    annotate("rect", xmin = voc_dates[5], xmax = voc_dates[6] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[3], alpha = 0.05) +
    annotate("rect", xmin = voc_dates[6], xmax = voc_dates[7] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[3], alpha = 0.15) +
    annotate("rect", xmin = voc_dates[7], xmax = voc_dates[8] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[3], alpha = 0.25) +
    
    # Omicron
    annotate("rect", xmin = voc_dates[8], xmax = voc_dates[9] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[4], alpha = 0.05) +
    annotate("rect", xmin = voc_dates[9], xmax = voc_dates[10] - 1,
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[4], alpha = 0.15) +
    annotate("rect", xmin = voc_dates[10], xmax = voc_dates[11],
             ymin = 0, ymax = ymax,
             fill = all_variant_colours[4], alpha = 0.25)
  
  if (label) {
    
    p <- p +
      annotate("label", x = voc_dates[c(1, 2, 5, 8)] + c(145, 95, 120, 40), y = ymax * 0.9, size = 4,
               label = c("Wildtype wave", "Alpha wave", "Delta wave", "Omicron wave"))
  }
  p
  
}


## TODO: check if this is used
plot_healthcare_box <- function(dat, r, vam_data, strain_epochs) {
  
  hosp_occup <- get_model_trajectory(dat, r, "hosp_occup", mean_quantiles = TRUE)
  hosp_occup$data <- dat$data[[r]]$full$general + dat$data[[r]]$full$icu
  
  occup <- gg_trajectory(hosp_occup, ymax = max(hosp_occup$ub),
                         ylab = "Hospital occupancy") +
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 0.7),
          axis.title.y = element_text(size = 10),
          axis.text.y =  element_text(size = 10),
          legend.text = element_text(size = 10))
  
  tx_labels <- data.frame(
    date = as.Date(c("2020-05-26", "2020-06-16", "2021-02-04", "2021-09-17", "2021-11-04", "2021-12-02", "2022-02-10")),
    label = c("Remdesivir", "Dexamethasone", "Anti-IL6", "Casirivimab-Imdevimab *", "Molnupiravir *", "Sotrovimab *", "Nirmatrelvir/ritonavir *"))
  
  date_vect <- data.frame(date = sircovid::sircovid_date_as_date(dat$samples[[1]]$trajectories$date[-1]))
  
  healthcare <- left_join(date_vect, tx_labels) %>%
    mutate(value = case_when(date %in% tx_labels$date[c(1, 3, 5, 7)] ~ 4e4,
                             date %in% tx_labels$date[c(2, 6)] ~ 3e4,
                             date %in% tx_labels$date[4] ~ 2e4,
                             TRUE ~ 4e4))
  
  strain_labels <- data.frame(
    date = strain_epochs,
    label = c("Alpha emerges", "Delta emerges", "Omicron emerges"))
  
  hc_tl <- occup +
    geom_vline(xintercept = strain_labels$date, lty = 2) +
    geom_vline(xintercept = tx_labels$date, lty = 3) +
    geom_label(aes(y = healthcare$value, label = healthcare$label))
  
  hc_tl <- annotate_strains(hc_tl, ymax = 4.5e4,
                            vam_data, r = "england", coff = c(0.1, 0.7, 0.9))
  
  hc_tl <- hc_tl + annotate(geom = "label", y = rep(4.25e4, 3), x = strain_epochs, label = paste(c("Alpha", "Delta", "Omicron"), "emerges"),
                            fill = all_variant_colours[2:4], alpha = 0.6)
  hc_tl
}


plot_inf_strain <- function(dat, r, variants) {
  
  sample <- dat$samples[[r]]
  data <- dat$data[[r]]
  date <- dat$info$date
  
  state <- sample$trajectories$state
  state <- state[, , -1L]
  
  x <- sircovid::sircovid_date_as_date(sample$trajectories$date)
  x <- x[-1L]
  
  inf_names <- c(rbind(variants, paste0(variants, " (reinf)")))
  
  inf_inc <- state[paste0("infections_inc_strain_", seq_len(4)), , ]
  
  strain_dates <- as.Date(c("2019-12-31", "2020-09-18", "2021-03-09", "2021-11-02"))
  strain_dates <- c(strain_dates, as.Date(date) + 1)
  
  res <- array(0, dim = c(length(variants) * 2, dim(inf_inc)[2:3]))
  row.names(res) <- inf_names
  
  res[variants[1], ,  x < strain_dates[2]] <-
    inf_inc["infections_inc_strain_1", , x < strain_dates[2]]
  
  for (i in seq_along(variants)[-1L]) {
    phase_dates <- (x >= strain_dates[i] & x < strain_dates[i + 1])
    res[variants[i - 1], , phase_dates] <-
      state["infections_inc_strain_1", , phase_dates]
    res[variants[i], , phase_dates] <-
      state["infections_inc_strain_2", , phase_dates]
    res[paste0(variants[i], " (reinf)"), , phase_dates] <-
      state["infections_inc_strain_3", , phase_dates]
    res[paste0(variants[i - 1], " (reinf)"), , phase_dates] <-
      state["infections_inc_strain_4", , phase_dates]
  }
  
  res <- apply(res, c(1, 3), mean)
  res <- as.data.frame(t(res))
  
  df <- data.frame(Total = rowSums(res)) %>%
    cbind(res, .) %>%
    mutate_at(vars(-Total), funs(. / Total)) %>%
    cbind("date" = x, .) %>%
    select(!Total) %>%
    pivot_longer(!date) %>%
    mutate(name = factor(name, levels = c(inf_names))) %>%
    mutate(reinf = case_when(stringr::str_detect(name, "reinf") ~ TRUE,
                             TRUE ~ FALSE))
  
  p <- ggplot(df, aes(date, value, fill = name, alpha = reinf)) +
    geom_area() +
    scale_alpha_discrete(range = c(0.35, 0.75), guide = FALSE) +
    geom_vline(xintercept = as.Date("2020-12-08"), lty = 2) +
    scale_fill_manual(values = rep(all_variant_colours, each = 2),
                      guide = guide_legend(
                        override.aes = list(alpha = rep(c(0.35, 0.75), 4)))) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "", y = "% daily infections")
  
  apply_theme(p) +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.7),
          legend.position = "bottom")
  
}


plot_vacc_status <- function (dat, r) {
  
  # Process trajectory for vaccine status
  vax_ro <- as.Date("2020-12-08")
  date <- sircovid::sircovid_date_as_date(dat$samples[[r]]$trajectories$date)[-1L]
  mod_traj_status <- get_model_trajectory(dat, r, "vaccine_status")
  vacc_status <- apply(mod_traj_status, 2, cumsum)
  vacc_status <- rbind(rep(0, dim(vacc_status)[2]), vacc_status)
  vacc_status <- data.frame(date = date, t(vacc_status) / vacc_status[nrow(vacc_status), ]) %>%
    filter(date >= as.Date("2020-11-01"))
  
  # # Dates at which coverage ≥50% for 1st, 2nd and 3rd doses
  slice_dates <- vacc_status$date[c(min(which(1 - vacc_status$vaccine_status_1 >= 0.5)),
                                    max(which(1 - vacc_status$vaccine_status_3 <= 0.5)),
                                    max(which(1 - vacc_status$vaccine_status_5 <= 0.5)))]
  
  # Some graphical parameters
  vacc_names <- factor(c("Unvaccinated", "First dose", "Second dose", "Booster"))
  cols <- c("#BBBBBB", "goldenrod2", "goldenrod3", "goldenrod4")
  names(cols) <- vacc_names
  alpha <- 0.85
  
  p <- ggplot(vacc_status, aes(x = date)) +
    geom_ribbon(aes(ymax = 1 - vaccine_status_1, ymin = 1 - vaccine_status_3, fill = vacc_names[2]), alpha = alpha) +
    geom_ribbon(aes(ymax = 1 - vaccine_status_3, ymin = 1 - vaccine_status_5, fill = vacc_names[3]), alpha = alpha) +
    geom_ribbon(aes(ymax = 1 - vaccine_status_5, ymin = 1 - vaccine_status_7, fill = vacc_names[4]), alpha = alpha) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0),
                       limits = c(0, 1)) +
    scale_fill_manual(values = cols[-1], breaks = vacc_names[-1]) +
    
    # Vacc rollout
    geom_vline(xintercept = vax_ro, lty = 2) +
    labs(x = "", y = "Vaccination status (%)")
  
  p <- apply_theme(p) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7, size = 10),
          legend.position = "bottom")
  
  p
}
