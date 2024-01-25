create_priors <- function(pars_info) {
  ##### derive priors for hospital progression
  regional_ps <- read_csv("weighted_prior_ranges.csv")
  regions <- c(sircovid::regions("england"), "england")

  regional_ps[regional_ps$param == "p_ICU", "mean"] <-
    0.5 * regional_ps[regional_ps$param == "p_ICU", "mean"]

  regional_ps <- regional_ps %>%
    dplyr::filter(region == "england") %>%
    dplyr::select(param, mean) %>%
    mutate(mean = as.numeric(mean))

  # Named vector of prior ranges for hospitalisation parameters
  ps_to_lower <- data.frame(
    param = c("p_ICU", "p_H_D", "p_ICU_D",
              "p_W_D", "p_G_D", "p_H", "p_H_2"),
    to_lower = c(0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3)
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
    scale = rep(NA, 3),
    shape = rep(NA, 3)
  )
  row.names(beta_hps) <- c("beta1", "beta2", "beta3")

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
  
  ## Between beta3 and beta18 we use the same prior distribution
  ## After beta19 we use the same prior distribution
  beta_hps <-
    beta_hps[c("beta1", "beta2", rep("beta3", 34)), ]
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
  
  pillar2_age_bands <- c("15_24", "25_49", "50_64", "65_79", "80_plus")
  
  ### add regional ps
  suppressMessages(
    hps <- hps %>%
      dplyr::full_join(regional_ps))

  ret <- priors_wide_to_long(hps)
  
  par <- c("rho_pillar2_tests", "alpha_H", "alpha_admission", "alpha_D",
           "alpha_death_hosp", "mu_D", "mu_D_2", "mu_D_3", "mu_D_4", "mu_D_5",
           "p_ICU_2", "p_H", "p_H_2", "p_G_D", "p_G_D_2", "mu_gamma_H", "mu_gamma_H_2",
           "mu_gamma_H_3", "mu_gamma_H_4", "seed_date_alpha", "ta_alpha",
           "seed_date_delta", "ta_delta", "seed_date_omicron", "ta_omicron",
           "rel_p_H_alpha", "rel_p_H_delta", "rel_p_ICU_alpha",
           "rel_p_ICU_delta", "rel_p_D_alpha", "rel_p_D_delta",
           "rel_p_H_omicron", "rel_p_ICU_omicron", "rel_p_D_omicron",
           paste0("p_NC_", pillar2_age_bands),
           paste0("p_NC_weekend_", pillar2_age_bands))
  
  if (assumptions == "mu_d_summer") {
    par <- c(par, "mu_D_6")
  }
  
  extra_uniform <-
    expand.grid(region = regions,
                type = "null",
                name = par,
                gamma_scale = NA_real_,
                gamma_shape = NA_real_,
                beta_shape1 = NA_real_,
                beta_shape2 = NA_real_,
                stringsAsFactors = FALSE)
  ret <- rbind(ret, extra_uniform)

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
