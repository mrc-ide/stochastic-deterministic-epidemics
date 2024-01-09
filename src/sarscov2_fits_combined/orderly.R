orderly2::orderly_parameters(date = NULL,
                             pars_set = NULL,
                             deterministic = FALSE,
                             adaptive_proposal = FALSE,
                             short_run = FALSE,
                             assumptions = "central")


for (r in sircovid::regions("england")) {
  orderly2::orderly_dependency("sarscov2_fits",
                               quote(latest(parameter:region == environment:r && parameter:date == this:date && parameter:pars_set == this:pars_set && parameter:assumptions == this:assumptions && parameter:short_run == this:short_run && parameter:deterministic == this:deterministic && parameter:adaptive_proposal == this:adaptive_proposal)),
                               c("regional_results/${r}/fit.rds" = "outputs/fit.rds",
                                 "regional_figs/pmcmc_traceplots_${r}.pdf" =  "outputs/pmcmc_traceplots.pdf"))
}

orderly2::orderly_artefact(
  "Files for external reviews",
  c("outputs/parameters/proposal.csv",
    "outputs/parameters/prior.csv",
    "outputs/parameters/info.csv",
    "outputs/aggregated_data.rds",
    "regional_results/Rt_england.rds"))
orderly2::orderly_artefact(
  "regional fitting plots and projections for comparison",
  c("figs/beta.png",
  "figs/data_fits_regional.png",
  "figs/forest_plot.png",
  "figs/forest_plot_betas.png",
  "figs/incidence.png",
  "figs/incidence_per_1000.png",
  "figs/pillar2_all_ages.png",
  "figs/pillar2_over25.png",
  "figs/react.png",
  "figs/Rt_eff_general.png",
  "figs/Rt_general.png",
  "figs/serology_euroimmun.png",
  "figs/serology_roche_n.png",
  "figs/variant_Alpha_Delta.png",
  "figs_by_age/pillar2_0_14.png",
  "figs_by_age/pillar2_15_24.png",
  "figs_by_age/pillar2_25_49.png",
  "figs_by_age/pillar2_50_64.png",
  "figs_by_age/pillar2_65_79.png",
  "figs_by_age/pillar2_80_plus.png",
  "figs_by_age/deaths_hosp_0_49.png",
  "figs_by_age/deaths_hosp_50_54.png",
  "figs_by_age/deaths_hosp_55_59.png",
  "figs_by_age/deaths_hosp_60_64.png",
  "figs_by_age/deaths_hosp_65_69.png",
  "figs_by_age/deaths_hosp_70_74.png",
  "figs_by_age/deaths_hosp_75_79.png",
  "figs_by_age/deaths_hosp_80_plus.png",
  "spim_view/regions.png",
  "spim_view/prevalence.png",
  "spim_view/pillar2_all_ages.png",
  "spim_view/pillar2_over25.png",
  "spim_view/pillar2_0_14.png",
  "spim_view/pillar2_15_24.png",
  "spim_view/pillar2_25_49.png",
  "spim_view/pillar2_50_64.png",
  "spim_view/pillar2_65_79.png",
  "spim_view/pillar2_80_plus.png",
  "spim_view/deaths_hosp_0_49.png",
  "spim_view/deaths_hosp_50_54.png",
  "spim_view/deaths_hosp_55_59.png",
  "spim_view/deaths_hosp_60_64.png",
  "spim_view/deaths_hosp_65_69.png",
  "spim_view/deaths_hosp_70_74.png",
  "spim_view/deaths_hosp_75_79.png",
  "spim_view/deaths_hosp_80_plus.png"))

library(sircovid)
library(spimalot)

orderly2::orderly_shared_resource(util.R = "util.R")
orderly2::orderly_resource("support.R")
source("support.R")
source("util.R")

version_check("sircovid", "0.14.13")
version_check("spimalot", "0.8.24")

sircovid_model <- "lancelot"
model_type <- "BB"

dat <- spimalot::spim_combined_load("regional_results", "england",
                                    get_onward = FALSE)


dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)
dir.create("figs_by_age", FALSE, TRUE)
dir.create("spim_view", FALSE, TRUE)

saveRDS(dat$data, "outputs/aggregated_data.rds")

saveRDS(dat$rt$england, "regional_results/Rt_england.rds")

spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")

write_png("figs/forest_plot.png",
          width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "non_betas"))

write_png("figs/forest_plot_betas.png", width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "betas"))

write_png("figs/data_fits_regional.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu",
              "general", "hosp", "all_admission"), age_band = "all",
            with_forecast = FALSE, add_betas = FALSE))


write_png("figs/serology_euroimmun.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(dat, sircovid::regions("england"), 1, 40))

write_png("figs/serology_roche_n.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(dat, sircovid::regions("england"), 2, 40))


pillar2_age_bands <-
  c("over25", "under15", "15_24", "25_49", "50_64", "65_79", "80_plus")
if (model_type == "BB") {
  write_png("figs/pillar2_all_ages.png", width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_pillar2_positivity(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date("2020-05-15"), ymax = 50))
  for (i in pillar2_age_bands) {
    if (i == "over25") {
      fig_name <- "figs/pillar2_over25.png"
    } else if (i == "under15") {
      fig_name <- "figs_by_age/pillar2_0_14.png"
    } else {
      fig_name <- paste0("figs_by_age/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_positivity(
                dat, sircovid::regions("england"), i,
                date_min = as.Date("2020-05-15"), ymax = 50))
  }
} else if (model_type == "NB") {
  write_png("figs/pillar2_all_ages.png", width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_pillar2_cases(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date("2020-05-15")))
  for (i in pillar2_age_bands) {
    if (i == "over25") {
      fig_name <- "figs/pillar2_over25.png"
    } else if (i == "under15") {
      fig_name <- "figs_by_age/pillar2_0_14.png"
    } else {
      fig_name <- paste0("figs_by_age/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_cases(
                dat, sircovid::regions("england"), i,
                date_min = as.Date("2020-05-15")))
  }
}

write_png("figs/react.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 10))

write_png("figs/variant_Alpha_Delta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"), "Delta",
            date_min = as.Date("2021-03-01"),
            date_max = as.Date("2021-08-15")))

write_png("figs/incidence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england")))

write_png("figs/incidence_per_1000.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england"), per_1000 = TRUE))

write_png("figs/Rt_eff_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"),
            "eff_Rt_general"))

write_png("figs/Rt_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "Rt_general"))

write_png("figs/beta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "beta"))


## add (zoomed in) plots of SPI-M-relevant trajectories
write_png("spim_view/regions.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths", "deaths_hosp", "icu", "general",
              "hosp", "all_admission"),
            date_min = as.Date(dat$info$date) - 75, age_band = "all",
            with_forecast = FALSE, add_betas = TRUE))

if (model_type == "BB") {
  write_png("spim_view/pillar2_all_ages.png", width = 2400, height = 1200,
            res = 200,
            spimalot::spim_plot_pillar2_positivity(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date(dat$info$date) - 75,
              ymax = 50, add_betas = TRUE))
  for (i in pillar2_age_bands) {
    if (i == "under15") {
      fig_name <- "spim_view/pillar2_0_14.png"
    } else {
      fig_name <- paste0("spim_view/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_positivity(
                dat, sircovid::regions("england"), i,
                date_min = as.Date(dat$info$date) - 75,
                ymax = 50, add_betas = TRUE))
  }
} else if (model_type == "NB") {
  write_png("spim_view/pillar2_all_ages.png", width = 2400, height = 1200,
            res = 200,
            spimalot::spim_plot_pillar2_cases(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date(dat$info$date) - 75,
              add_betas = TRUE))
  for (i in pillar2_age_bands) {
    if (i == "under15") {
      fig_name <- "spim_view/pillar2_0_14.png"
    } else {
      fig_name <- paste0("spim_view/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_cases(
                dat, sircovid::regions("england"), i,
                date_min = as.Date(dat$info$date) - 75,
                add_betas = TRUE))
  }
}

write_png("spim_view/prevalence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"),
            date_min = as.Date(dat$info$date) - 75,
            ymax = 10, add_betas = TRUE))

## Plot outputs by age

deaths_hosp_age_bands <- c("0_49", "50_54", "55_59", "60_64", "65_69", "70_74",
                           "75_79", "80_plus")
for (i in deaths_hosp_age_bands) {
  fig_name <- paste0("figs_by_age/deaths_hosp_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, sircovid::regions("england"), "deaths_hosp", age_band = i,
              with_forecast = FALSE, add_betas = FALSE))
  
  fig_name <- paste0("spim_view/deaths_hosp_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, sircovid::regions("england"), "deaths_hosp", age_band = i,
              date_min = as.Date(dat$info$date) - 75,
              with_forecast = FALSE, add_betas = FALSE))
}
