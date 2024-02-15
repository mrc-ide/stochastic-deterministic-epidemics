plot_mu_D <- function(dat, regions) {
  
  oo <- par(mfrow = c(2, ceiling(length(regions) / 2)), oma = c(2, 1, 2, 1),
            mar = c(3, 3, 3, 1))
  on.exit(par(oo))
  
  for (r in regions) {
    plot_mu_D_region(r, dat)
  }
}

plot_mu_D_region <- function(region, dat) {
  
  samples <- dat$samples[[region]]
  
  cols <- spimalot:::spim_colours()
  
  time <- samples$trajectories$time
  
  get_mu_D <- function(i) {
    p <- samples$predict$transform(samples$pars[i, ])
    p <- p[[length(p)]]$pars
    j_max <- which.max(p$p_H_D_step[1, 1:17])
    mu_D_step <- p$p_H_D_step[, j_max] / p$p_H_D_step[1, j_max]
    
    sircovid::sircovid_parameters_expand_step(time, mu_D_step)
  }
  
  res <- 
    t(vapply(seq_len(nrow(samples$pars)), get_mu_D, numeric(length(time))))
  
  x <- sircovid::sircovid_date_as_date(samples$trajectories$date)
  colnames(res) <- as.character(x)
  
  ps <- seq(0.025, 0.975, 0.005)
  qs <- apply(res, MARGIN = 2, FUN = quantile, ps, na.rm = TRUE)
  
  oo <- par(mgp = c(1.7, 0.5, 0), bty = "n")
  on.exit(oo)
  
  
  xlim <- c(min(x), max(x))
  ylim <- c(0, 1)
  
  plot(xlim[1], 0, type = "n",
       xlim = xlim,
       ylim = ylim,
       main = toupper(spimalot::spim_region_name(region)),
       font.main = 1,
       xlab = "", ylab = expression(mu[D]),
       axes = FALSE,
       xaxs = "i",
       yaxs = "i")
  
  mu_cols <- c(spimalot:::mix_cols(cols$now, "white", 0.7),
               spimalot:::mix_cols(cols$now, "white", 0.495))
  
  #Extract every first date of month from x
  firsts <- x[!duplicated(substring(x, 1, 7))]
  #Extract every first date of year from x
  year_firsts <- x[!duplicated(substring(x, 1, 4))]
  abline(v = firsts[-1], col = "ivory2") #Plot gray line on 1st of every month
  abline(v = year_firsts[-1], col = "gray")
  axis(side = 1, at = pretty(xlim), labels = format(pretty(xlim), "%b-%y"))
  axis(side = 2, at = pretty(ylim), labels = pretty(ylim))
  
  
  spimalot:::ci_bands(qs[c("2.5%", "25.0%", "75.0%", "97.5%"), ], x,
                      cols = mu_cols, horiz = FALSE, leg = FALSE)
  lines(x, qs["50.0%", ], col = cols$now, lty = 1, lwd = 1.5,
        lend = 1)
  
}
