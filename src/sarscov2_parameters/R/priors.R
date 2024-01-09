create_priors <- function(pars_info) {
  ##### derive priors for hospital progression
  regional_ps <- read_csv("weighted_prior_ranges.csv")
  regions <- c(sircovid::regions("england"))

  regional_ps[regional_ps$param == "p_ICU", "mean"] <-
    0.5 * regional_ps[regional_ps$param == "p_ICU", "mean"]

  additional_ps <- t(data.frame(c(param = "p_H", mean = 0.75),
                                c(param = "p_H_2", mean = 0.75),
                                c(param = "p_H_3", mean = 0.75)))
  rownames(additional_ps) <- c()

  regional_ps <- regional_ps %>%
    dplyr::filter(region == "england") %>%
    dplyr::select(param, mean) %>% rbind(additional_ps) %>%
    mutate(mean = as.numeric(mean))

  # Named vector of prior ranges for hospitalisation parameters
  ps_to_lower <- data.frame(
    param = c("p_ICU", "p_H_D", "p_ICU_D",
              "p_W_D", "p_G_D", "p_H", "p_H_2", "p_H_3"),
    to_lower = c(0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3)
  )

  regional_ps <- regional_ps %>%
    dplyr::inner_join(ps_to_lower)

  p_hps <- signif(mapply(FUN = fit_beta,
                         mean = regional_ps$mean,
                         lower = regional_ps$mean - regional_ps$to_lower,
                         ci = 0.95))
  p_hps <- as.data.frame(t(p_hps))
  colnames(p_hps) <- c("shape1","shape2")
  regional_ps <- cbind(regional_ps,p_hps)

  regional_ps <- 
    do.call("rbind", replicate(length(regions), regional_ps, simplify = FALSE))
  regional_ps <- regional_ps %>% 
    mutate(region = rep(regions, each = nrow(p_hps))) %>%
    dplyr::select(par = "param", region, mean, shape1, shape2)

  ### set beta_priors
  beta_hps <- data.frame(
    scale = rep(NA, 4),
    shape = rep(NA, 4)
  )
  row.names(beta_hps) <- c("beta1", "beta2", "beta3", "beta19")

  ## beta value that would give R0 = 1
  R0_fac <- 0.0367

  ## beta1 aim for 95% CI of [2.5, 3.5]
  beta_hps["beta1", ] <- fit_gamma(mean_D = 2.979,
                                   lower_D = 2.5,
                                   upper_D = 3.5,
                                   ci = 0.95)
  beta_hps["beta1", "scale"] <- beta_hps["beta1", "scale"] * R0_fac
  ## beta2 aim for 95% CI of [0.4, 3.5]
  beta_hps["beta2", ] <- fit_gamma(mean_D = 1.562,
                                   lower_D = 0.4,
                                   upper_D = 3.5,
                                   ci = 0.95)
  beta_hps["beta2", "scale"] <- beta_hps["beta2", "scale"] * R0_fac
  ## beta3 aim for 95% CI of [0.4, 3]
  beta_hps["beta3", ] <- fit_gamma(mean_D = 1.395,
                                   lower_D = 0.4,
                                   upper_D = 3,
                                   ci = 0.95)
  beta_hps["beta3", "scale"] <- beta_hps["beta3", "scale"] * R0_fac
  ## beta19 aim for 95% CI of [0.4, 5.4]
  beta_hps["beta19", ] <- fit_gamma(mean_D = 2.168,
                                    lower_D = 0.4,
                                    upper_D = 5.4,
                                    ci = 0.95)
  beta_hps["beta19", "scale"] <- beta_hps["beta19", "scale"] * R0_fac

  ## Between beta3 and beta18 we use the same prior distribution
  ## After beta19 we use the same prior distribution
  beta_hps <-
    beta_hps[c("beta1", "beta2", rep("beta3", 16), rep("beta19", 9)), ]
  rownames(beta_hps) <- paste0("beta", seq_len(nrow(beta_hps)))
  beta_names <- rownames(beta_hps)

  pars <- c(beta_names, unique(regional_ps$par))

  hps <- matrix(NA, nrow = length(pars), ncol = 7,
                dimnames = list(pars, c("par", "region", "scale", "shape",
                                        "shape1", "shape2", "correlation")))
  hps <- as.data.frame(hps)
  hps$par <- pars
  hps$region <- "england"
  hps[beta_names, colnames(beta_hps)] <- beta_hps
  hps[unique(regional_ps$par), c("shape1", "shape2")] <- as.matrix(p_hps)
  
  ## uniform beta priors
  pillar2_age_bands <- c("under15", "15_24", "25_49", "50_64", "65_79", "80_plus")
  par <- c("eps", "rho_pillar2_tests", "alpha_H", "alpha_D", "mu_D", "mu_D_2",
           "p_ICU_2", "p_G_D", "p_G_D_CHR", "p_NC",
           "p_NC_weekend", "alpha_pillar2_cases", "phi_pillar2_cases",
           "phi_pillar2_cases_weekend", "alpha_death_hosp",
           paste0("phi_pillar2_cases_", pillar2_age_bands),
           paste0("phi_pillar2_cases_weekend_", pillar2_age_bands))

  extra_beta <- data.frame(
    par = par,
    region = "england",
    scale = NA_real_,
    shape = NA_real_,
    shape1 = 1,
    shape2 = 1,
    correlation = NA_real_,
    stringsAsFactors = FALSE)
  rownames(extra_beta) <- extra_beta$par

  ## pseudouniform gamma priors
  extra_gamma_flat <- data.frame(
    par = c("mu_gamma_H", "mu_gamma_H_2", "mu_gamma_H_3", "mu_gamma_H_4"),
    region = "england",
    scale = 1000,
    shape = 1.001,
    shape1 = NA_real_,
    shape2 = NA_real_,
    correlation = NA_real_,
    stringsAsFactors = FALSE)
  rownames(extra_gamma_flat) <- extra_gamma_flat$par

  hps <- rbind(hps, extra_beta, extra_gamma_flat)

  ## For "start_date" we have a uniform prior on sircovid_dates
  ## 1,...,75. This is fine for restarting purposes as this parameter
  ## will inevitably not be fitted after a restart.

  ### add regional ps
  suppressMessages(
    hps <- hps %>%
      dplyr::full_join(regional_ps))

  ## gamma priors on m_CHW, m_CHR
  regions <- unique(hps$region)
  chr_pop <- sapply(regions,
                    function(x){sircovid::lancelot_parameters(1, x)$N_tot[19]})
  ch_mean <- 1e-1
  ch_shape <- 5
  m_CHW_hps <- data.frame(
    par = "m_CHW",
    region = regions,
    scale = ch_mean / (ch_shape * chr_pop), #scale by CHR population
    shape = ch_shape,
    shape1 = NA_real_,
    shape2 = NA_real_,
    correlation = NA_real_,
    mean = NA_real_,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  m_CHR_hps <- m_CHW_hps
  m_CHR_hps$par <- "m_CHR"
  hps <- rbind(hps, m_CHW_hps, m_CHR_hps)

  ret <- priors_wide_to_long(hps)

  extra_multistrain <-
    expand.grid(region = unique(pars_info$region),
                type = "null",
                name = c("seed_date_delta", "ta_delta",
                         "rel_p_H_delta", "rel_p_ICU_delta",
                         "rel_p_D_delta"),
                gamma_scale = NA_real_,
                gamma_shape = NA_real_,
                beta_shape1 = NA_real_,
                beta_shape2 = NA_real_,
                stringsAsFactors = FALSE)
  ret <- rbind(ret, extra_multistrain)

  nms_expected <- unique(pars_info$name)
  nms_found <- unique(ret$name)
  msg <- setdiff(nms_expected, nms_found)
  if (length(msg) > 0) {
    stop(sprintf("Missing parameters, update priors (missing %s)",
                 paste(msg, collapse = ", ")))
  }
  extra <- setdiff(nms_found, nms_expected)
  if (length(extra)) {
    message(sprintf("Dropping %d unused priors: %s",
                    length(extra), paste(extra, collapse = ", ")))
    ret <- ret[ret$name %in% nms_expected, ]
  }
  rownames(ret) <- NULL

  invisible(ret)
}


## specify gamma in terms of mean and variance rather than shape and scale

# convert mean and variance of gamma to shape and scale
mv_to_ss <- function(mean, var) {
  scale <- var / mean
  shape <- mean / scale
  list(shape = shape, scale = scale)
}

# gamma dist functions
qgammamv <- function(p, mean, var) {
  X <- mv_to_ss(mean, var)
  qgamma(p = p, shape = X$shape, scale = X$scale)
}

dgammav <- function(x, mean, var) {
  X <- mv_to_ss(mean, var)
  dgamma(x = x, shape = X$shape, scale = X$scale)
}


## fitting function by least-squares based on mean and CIs
fit_gamma <- function(mean_D, lower_D, upper_D, ci = 0.99) {

  alpha <- (1 - ci)/2
  p <- c(alpha , 1 - alpha)

  f <- function(v) {
    x <- qgammamv(p = p, mean = mean_D, var = v)
    sum((x[1] - c(lower_D))^2)
  }

  var_D <- optimise(f = f, interval = c(0,10), maximum = FALSE)$minimum
  X <- mv_to_ss(mean_D, var_D)
  message(paste(c("fitted qs =", round(qgamma(p, shape = X$shape, scale = X$scale),3)), collapse = " "))
  message(paste(c("target qs =", c(lower_D, upper_D)), collapse = " "))
  message(paste("fitted var =", round(var_D, 3)))
  c(scale = X$scale, shape = X$shape)

}


fit_beta <- function(mean, lower, ci = 0.99) {
  a <- (1 - ci)/2
  p <- c(a , 1 - a)

  f <- function(alpha) {
    beta <- alpha*(1 - mean) / mean
    x <- qbeta(p = p,shape1 = alpha, shape2 = beta)
    sum((x[1] - c(lower))^2)
  }

  alpha <- optimise(f = f, interval = c(0,1e3), maximum = FALSE)$minimum
  beta <- alpha*(1 - mean) / mean

  message(paste(c("fitted qs =", round(qbeta(p, shape1 = alpha, shape2 =  beta),3)), collapse = " "))

  c(alpha = alpha, beta = beta[[1]])

}


priors_wide_to_long <- function(d) {
  stopifnot(all(xor(is.na(d$shape1), is.na(d$shape))))
  d$type <- ifelse(is.na(d$shape1), "gamma", "beta")

  tr <- c(par = "name",
          scale = "gamma_scale",
          shape = "gamma_shape",
          shape1 = "beta_shape1",
          shape2 = "beta_shape2")
  d <- sircovid:::rename(d, names(tr), tr)
  d <- d[c("region", "type", tr)]

  extra <- data.frame(
    region = "england",
    type = "null",
    name = "start_date",
    gamma_scale = NA_real_,
    gamma_shape = NA_real_,
    beta_shape1 = NA_real_,
    beta_shape2 = NA_real_,
    stringsAsFactors = FALSE)

  d <- rbind(d, extra)

  ## Not all parameters are region-specific, let's fix that too:
  f <- function(p) {
    s <- d[d$name == p, ]
    msg <- setdiff(d$region, s$region)
    if (length(msg) > 0) {
      i <- match("england", s$region)
      extra <- s[rep(i, length(msg)), ]
      extra$region <- msg
      s <- rbind(s, extra)
    }
    s
  }

  res <- do.call(rbind, lapply(unique(d$name), f))

  ## We must have all parameters for all regions, and no doubles
  stopifnot(all(table(res$region, res$name) == 1))

  res <- res[order(res$region, res$name), ]
  rownames(res) <- NULL
  res
}
