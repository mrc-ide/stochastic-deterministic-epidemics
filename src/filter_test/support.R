colour_covariance <- function(x) {
  colors <- c("#A50F15", "#DE2D26", "#FB6A4A", "#FCAE91", "#FEE5D9", "white",
              "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")
  i <- !is.na(x)
  cols <- colorRamp(colors)((x[i] + 1) / 2)
  cols <- rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 255)
  ret <- rep(NA_character_, length(x))
  ret[i] <- cols
  ret
}

## run particle filter on parameter set
filtered_output <- function(pars_index,
							n_particles_seq,
							n_filter_iter,
							filter_data,
							mcmc_pars_deterministic,
							pmcmc_run_deterministic,
							mcmc_pars_stochastic,
							pmcmc_run_stochastic,
							compare,
							sir_data
) {
  
  ## run for the length of particles given
  for (i in 1:length(n_particles_seq)) {
		
	## Build stochastic filter
	filter_stochastic <- mcstate::particle_filter$new(data = sir_data,
													model = gen_sir,
													n_particles = n_particles_seq[i],
													n_threads = 4,
													compare = case_compare)
	
	## store deterministic samples filtered
	deterministic_filtered <- list()
	
	## run filter on deterministic parameter set
	deterministic_filtered <- lapply(1:n_filter_iter, function(i) { filter_stochastic$run(
	mcmc_pars_deterministic$model(
	pmcmc_run_deterministic$pars[pars_index, ]
	)
	) + mcmc_pars_deterministic$prior(mcmc_pars_deterministic$initial()
	)
	}
	)
	
	
	deterministic_filtered <- unlist(deterministic_filtered)
	
	## store stochastic samples filtered
	stochastic_filtered <- list()
	
	## run filter on stochastic parameter set
	stochastic_filtered <- lapply(1:n_filter_iter, function(i) { filter_stochastic$run(
	  mcmc_pars_stochastic$model(
		pmcmc_run_stochastic$pars[pars_index, ]
	  )
	) + mcmc_pars_stochastic$prior(mcmc_pars_stochastic$initial()
	)
	}
	)
	
	stochastic_filtered <- unlist(stochastic_filtered)
	
	## save filter outputs
	filter_data[[i]] <- list(deterministic_filtered, stochastic_filtered)
	names(filter_data[[i]]) <- c("deterministic_filtered", "stochastic_filtered")
  }
  
  names(filter_data) <- paste0("n_", n_particles_seq)
  filter_data
}
