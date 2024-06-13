plot_traceplots <- function(samples) {
  if (is.null(samples$chain)) {
    n_chains <- 1L
  } else {
    ## We assume this below
    n_chains <- length(unique(samples$chain))
    stopifnot(
      identical(samples$chain,
                rep(seq_len(n_chains),
                    each = length(samples$chain) / n_chains)))
  }
  cols <- rev(viridisLite::viridis(n_chains))
  
  pars <- samples$pars_full

  nms <- colnames(pars)
  probs <- samples$probabilities_full
  
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
            font.main = 1)
  }
  
  
  new_grid(length(nms), FALSE)
  for (nm in nms) {
    plot_traces1(samples$pars_full[, nm], nm)
  }
  plot_traces1(probs[, "log_likelihood"], "log_likelihood")
}
