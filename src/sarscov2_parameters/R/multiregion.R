## Need some utilities here; these functions likely headed to spimalot
## at some point:
assert_has_names <- spimalot:::assert_has_names
assert_unique <- spimalot:::assert_unique
assert_setequal <- spimalot:::assert_setequal

independent_to_multiregion <- function(pars) {
  info <- independent_to_multiregion_info(pars$info)
  prior <- independent_to_multiregion_prior(info, pars$prior)
  proposal <- independent_to_multiregion_proposal(info, pars$proposal)
  list(info = info, prior = prior, proposal = proposal)
}


independent_to_multiregion_info <- function(info) {
  regions <- unique(info$region)
  assert_has_names(info, c("region", "include",
                           "name", "initial", "max", "integer"))

  ## I have no idea why we have these, but knock them out now.
  info <- info[info$include, ]
  rownames(info) <- NULL

  nms <- split(info$name, info$region)
  if (length(unique(nms)) != 1) {
    stop("Names differ across regions")
  }
  assert_unique(nms[[1]], "info$name")

  re_varied <- "^(beta[0-9]+$|seed_date_|start_date$|p_ICU)"
  is_varied <- grepl(re_varied, info$name)

  info_varied <- info[is_varied, ]

  ## More work to do for fixed things:
  f <- function(x) {
    x$region <- NA
    x$initial <- mean(x$initial)
    for (v in c("min", "max", "integer", "varied")) {
      if (!all(x[[v]] == x[[v]][[1]])) {
        stop(sprintf("Inconsistent '%s' across regions for '%s'",
                     v, x$name[[1]]))
      }
    }
    x[1, ]
  }
  info_fixed <- info[!is_varied, ]
  info_fixed <- do.call(rbind, lapply(split(info_fixed, info_fixed$name), f))

  ret <- rbind(info_fixed, info_varied)
  rownames(ret) <- NULL

  ret
}


independent_to_multiregion_prior <- function(info, prior) {
  prior_cols <- c("region", "type", "name", "gamma_scale", "gamma_shape",
                  "beta_shape1", "beta_shape2")
  assert_has_names(prior, prior_cols)

  prior <- prior[prior$region %in% info$region, ]

  nms_varied <- unique(info$name[!is.na(info$region)])
  nms_fixed <- info$name[is.na(info$region)]
  assert_setequal(prior$name, c(nms_varied, nms_fixed), unique = FALSE)

  nms <- split(prior$name, prior$region)
  if (length(unique(nms)) != 1) {
    stop("Names differ across regions")
  }
  
  prior_varied <- prior[prior$name %in% nms_varied, ]

  ## Then for the fixed, more complicated:
  f <- function(x) {
    is_same <- function(x) {
      all(is.na(x)) || (!anyNA(x) && all(x == x[[1]]))
    }
    nm <- x$name[[1]]
    ok <- (nm %in% c("m_CHR", "m_CHW")) ||
      all(vapply(x[names(x) != "region"], is_same, TRUE))
    if (!ok) {
      stop(sprintf("prior disagreement for '%s'", nm))
    }
    x[1, ]
  }

  prior_fixed <- prior[prior$name %in% nms_fixed , ]
  prior_fixed <- dplyr::bind_rows(
    lapply(split(prior_fixed, prior_fixed$name), f))
  prior_fixed$region[] <- NA

  prior <- rbind(prior_fixed, prior_varied)
  rownames(prior) <- NULL

  prior
}


independent_to_multiregion_proposal <- function(info, proposal) {
  nms_varied <- unique(info$name[!is.na(info$region)])
  nms_fixed <- info$name[is.na(info$region)]
  nms_all <- c(nms_fixed, nms_varied)
  assert_setequal(proposal$name, nms_all, unique = FALSE)
  assert_setequal(setdiff(names(proposal), c("name", "region")), nms_all)

  regions <- unique(info$region[!is.na(info$region)])

  ## We can't easily create vcv matrices for subsets, so we'll just
  ## create a matrix of variances:

  f <- function(p) {
    diag(as.matrix(p[match(nms_all, p$name), nms_all]))
  }
  var <- vapply(split(proposal, proposal$region), f, numeric(length(nms_all)))
  rownames(var) <- nms_all

  f <- function(r) {
    m_varied <- diag(var[nms_varied, r])
    colnames(m_varied) <- nms_varied
    m_fixed <- matrix(0, length(nms_varied), length(nms_fixed))
    colnames(m_fixed) <- nms_fixed
    common <- data.frame(region = r, name = nms_varied,
                         stringsAsFactors = FALSE)
    cbind(common, m_varied, m_fixed)
  }
  proposal_varied <- do.call(rbind, lapply(colnames(var), f))

  m_varied <- matrix(0, length(nms_fixed), length(nms_varied))
  colnames(m_varied) <- nms_varied
  m_fixed <- diag(apply(var[nms_fixed, ], 1, function(x) mean(x[x > 0])))
  colnames(m_fixed) <- nms_fixed
  common <- data.frame(region = NA_character_, name = nms_fixed)
  proposal_fixed <- cbind(common, m_varied, m_fixed)

  ret <- rbind(proposal_fixed, proposal_varied)
  rownames(ret) <- NULL
  ret
}
