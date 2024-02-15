load_mcmc_parameters <- function(assumptions, deterministic) {

  message(sprintf("Assumptions are with '%s' parameters", assumptions))
  message(sprintf("Model will be run in the '%s' mode",
                  if (deterministic) "deterministic" else "stochastic"))
  
  path_pars <- file.path("pars", assumptions,
                         if (deterministic) "deterministic" else "stochastic")

  info <- read_csv(file.path(path_pars, "info.csv"))
  ## "discrete" has been deprecated in mcstate and replaced by "integer"
  names(info)[names(info) == "discrete"] <- "integer"
  proposal <- read_csv(file.path(path_pars, "proposal.csv"))
  prior <- create_priors(info)
  proposal <- update_proposal(info, proposal)

  ret <- list(info = info,
              proposal = proposal,
              prior = prior)
  
  ret
}


add_parameter <- function(pars, name, initial, min, max, proposal, integer) {
  ## info , prior, proposal
  regions <- unique(pars$info$region)
  stopifnot(
    length(initial) %in% c(1, length(regions)),
    length(proposal) %in% c(1, length(regions)))
  info_new <- data.frame(
    region = regions,
    name = name,
    initial = initial,
    min = min,
    max = max,
    integer = integer,
    include = TRUE)

  info <- rbind(pars$info, info_new)
  info <- info[order(info$region, info$name), ]

  proposal_values <- proposal # TODO: rewrite to avoid churn later
  proposal <- pars$proposal
  proposal[[name]] <- 0
  k <- proposal$name == proposal$name[1]
  proposal_new <- proposal[k, ]
  proposal_new$name <- name
  proposal_new[, -c(1, 2)] <- 0
  proposal_new[[name]] <- proposal_values

  proposal <- rbind(proposal, proposal_new)
  proposal <- proposal[order(proposal$region, proposal$name), ]
  col_order <- c(1, 2, 2 + order(names(proposal)[-c(1, 2)]))
  proposal <- proposal[, col_order]

  prior_new <- data.frame(
    region = regions,
    type = "null",
    name = name,
    gamma_scale = NA,
    gamma_shape = NA,
    beta_shape1 = NA,
    beta_shape2 = NA)
  prior <- rbind(pars$prior, prior_new)
  prior <- prior[order(prior$region, prior$name), ]

  list(info = info, prior = prior, proposal = proposal)
}


spimalot_file <- function(path) {
  system.file(path, package = "spimalot", mustWork = TRUE)
}


get_daily_doses <- function(region, scenario,
                            daily_doses_per_region = NULL, date,
                            future_doses_by_region = NULL) {
  date_start <- as.character(date)
  date_start_plus_two <- as.character(date + 2)
  date_start_plus_three <- as.character(date + 3)
  if (scenario == "pessimistic") {
    ret <- c(daily_doses_per_region[[region]])
    names(ret) <- c(date_start_plus_two)
  } else if (scenario == "central") {
    if(as.Date(date_start_plus_two) < as.Date("2021-01-25")) {
      ret <- c(daily_doses_per_region[[region]],
               2 * daily_doses_per_region[[region]])
      names(ret) <- c(date_start_plus_two,
                      "2021-01-25")
    } else {
      ret <- c(2 * daily_doses_per_region[[region]])
      names(ret) <- c(date_start_plus_two)
    }
  } else if (scenario == "optimistic") {

    if(as.Date(date_start_plus_two) < as.Date("2021-01-25")) {
      ret <- c(daily_doses_per_region[[region]],
               2 * daily_doses_per_region[[region]],
               3 * daily_doses_per_region[[region]])
      names(ret) <- c(date_start_plus_two,
                      "2021-01-25",
                      "2021-02-01")
    } else if(as.Date(date_start_plus_two) < as.Date("2021-02-01")) {
      ret <- c(2 * daily_doses_per_region[[region]],
               3 * daily_doses_per_region[[region]])
      names(ret) <- c(date_start_plus_two,
                      "2021-02-01")
    } else {
      ret <- c(3 * daily_doses_per_region[[region]])
      names(ret) <- c(date_start_plus_two)
    }
  }  else if (scenario == "counterfactual") {
    ret <- c(0)
    names(ret) <- c(date_start_plus_two)
  } else if (scenario == "manual") {

    future <- round(future_doses_by_region[[region]] / 7)
    future_date <- future_doses_by_region$week_start
    if (any(diff(future_date) <= 0)) {
      stop("future_date$week_start must be strictly increasing")
    }
    min_future_date <- date + 2
    ret <- check_future_dates(future, future_date, min_future_date)
  } else {
    stop("Unknown scenario")
  }

  ## All dates must be increasing; just checks the logic above
  if (length(ret) > 1) {
    stopifnot(all(diff(as.Date(names(ret))) > 0))
  }


  ret
}


check_future_dates <- function(future, future_date, min_future_date) {
  if (any(future_date < min_future_date)) {
    i <- max(which(future_date <= min_future_date))
    future <- future[i:length(future)]
    future_date <- future_date[i:length(future_date)]
    future_date[[1]] <- min_future_date
  }
  names(future) <- as.character(future_date)
  future
}


get_efficacy_against_disease_conditional_on_infection <-
  function(efficacy_against_disease_total, efficacy_against_acquisition) {
    if(any(efficacy_against_disease_total < efficacy_against_acquisition))
      stop("efficacy_against_disease_total cannot be lower than efficacy_against_acquisition")
    (efficacy_against_disease_total - efficacy_against_acquisition) / (1 - efficacy_against_acquisition)
}


split_doses_by_region <- function(future_doses,
                                  mean_time_between_doses,
                                  regional_prop_pop) {
  week_start <- as.Date(future_doses$week_start, "%d-%B-%y")
  future_doses$week_start <- week_start

  message("Using following weekly vaccine roll-out schedule(s): ")
  msg <- apply(future_doses, 1, paste, collapse= ": ")
  message(paste(msg, collapse = "\n"))

  doses_per_week_cols <- grep("doses_per_week", names(future_doses), value = TRUE)

  if (length(doses_per_week_cols) > 1) {
    doses_per_week <- lapply(doses_per_week_cols, function(e) future_doses[[e]])
    names(doses_per_week) <- gsub("doses_per_week_", "", doses_per_week_cols)
  } else {
    doses_per_week <- list(main = future_doses$doses_per_week)
  }

  ### split future_doses by region
  future_doses_by_region <- lapply(doses_per_week, function(e) {
    x <- as.data.frame(matrix(
      vapply(regional_prop_pop, function(r) e * r, numeric(length(e))),
      nrow = length(e)))
    names(x) <- names(regional_prop_pop)
    x$week_start <- future_doses$week_start
    x
  })

  future_doses_by_region
}


update_proposal <- function(info, proposal) {
  msg <- setdiff(info$name, proposal$name)
  if (length(msg) > 0) {
    message(sprintf("Adding %d parameters to proposal matrix", length(msg)))
    extra <- expand.grid(region = unique(proposal$region), name = msg)
    extra[names(proposal)[-(1:2)]] <- 0
    proposal <- rbind(proposal, extra)
    proposal[msg] <- 0
    proposal <- proposal[order(proposal$region, proposal$name), ]
    col_order <- c(1, 2, 2 + order(names(proposal)[-c(1, 2)]))
    proposal <- proposal[, col_order, drop = FALSE]
  }
  proposal
}
