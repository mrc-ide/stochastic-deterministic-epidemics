orderly2::orderly_parameters(short_run = FALSE, deterministic = FALSE, assumptions = "central")

orderly2::orderly_shared_resource(util.R = "util.R")

orderly2::orderly_resource("paper_numbers.Rmd")

for (r in sircovid::regions("england")) {
  orderly2::orderly_dependency(
    "sarscov2_fits",
    quote(latest(parameter:region == environment:r && parameter:assumptions == this:assumptions && parameter:short_run == this:short_run && parameter:deterministic == this:deterministic)),
    c("regional_results/${r}/fit.rds" = "outputs/fit.rds",
      "regional_figs/pmcmc_traceplots_${r}.pdf" = "outputs/pmcmc_traceplots.pdf",
      "regional_figs/multipage/pmcmc_traceplots_${r}_separate.pdf" = "outputs/pmcmc_traceplots_separate.pdf"))  
}

orderly2::orderly_dependency(
  "sarscov2_data",
  "latest",
  c(severity_data.csv = "outputs/severity_data.csv",
    hfr_week.csv = "outputs/hfr_week.csv",
    england_region_data.csv = "outputs/england_region_data.csv",
    outcomes_vacc_status.csv = "outputs/outcomes_vacc_status.csv"))

orderly2::orderly_artefact("Paper numbers", "paper_numbers.html")
orderly2::orderly_artefact(
  "Files for external reviews",
  c("outputs/parameters/proposal.csv", "outputs/parameters/prior.csv", "outputs/parameters/info.csv", "outputs/model_demography.rds", "regional_results/Rt_england.rds", "regional_results/Rt_all.rds", "outputs/aggregated_data.rds", "outputs/intrinsic_severity.csv", "outputs/national_severity.rds"))
orderly2::orderly_artefact("Files for forward simulation", "outputs/combined.rds")
orderly2::orderly_artefact(
  "regional fitting plots and projections for comparison",
  c("figs/beta.png", "figs/cumulative_attack_rate.png", "figs/data_fits.png", "figs/forest_plot_betas.png", "figs/forest_plot_misc.png", "figs/forest_plot_variants.png", "figs/forest_plot_tv_severity.png", "figs/incidence.png", "figs/incidence_per_1000.png", "figs/infections_per_strain.png", "figs/pillar2_all_ages.png", "figs/pillar2_over25.png", "figs/prevalence_react.png", "figs/prevalence_ons.png", "figs/Rt_eff_general.png", "figs/Rt_general.png", "figs/serology_euroimmun.png", "figs/serology_roche_n.png", "figs/status_effective_susceptible.png", "figs/status_infection.png", "figs/status_vaccine.png", "figs/variant_Wildtype_Alpha.png", "figs/variant_Alpha_Delta.png", "figs/variant_Delta_Omicron.png", "figs/mu_D.png", "figs/traceplot_east_of_england.png", "figs/traceplot_london.png", "figs/traceplot_midlands.png", "figs/traceplot_north_east_and_yorkshire.png", "figs/traceplot_north_west.png", "figs/traceplot_south_east.png", "figs/traceplot_south_west.png", "figs_by_age/pillar2_0_14.png", "figs_by_age/pillar2_15_24.png", "figs_by_age/pillar2_25_49.png", "figs_by_age/pillar2_50_64.png", "figs_by_age/pillar2_65_79.png", "figs_by_age/pillar2_80_plus.png", "figs_by_age/deaths_hosp_0_49.png", "figs_by_age/deaths_hosp_50_54.png", "figs_by_age/deaths_hosp_55_59.png", "figs_by_age/deaths_hosp_60_64.png", "figs_by_age/deaths_hosp_65_69.png", "figs_by_age/deaths_hosp_70_74.png", "figs_by_age/deaths_hosp_75_79.png", "figs_by_age/deaths_hosp_80_plus.png", "figs_by_age/deaths_comm_0_49.png", "figs_by_age/deaths_comm_50_54.png", "figs_by_age/deaths_comm_55_59.png", "figs_by_age/deaths_comm_60_64.png", "figs_by_age/deaths_comm_65_69.png", "figs_by_age/deaths_comm_70_74.png", "figs_by_age/deaths_comm_75_79.png", "figs_by_age/deaths_comm_80_plus.png", "figs_by_age/admissions_0_9.png", "figs_by_age/admissions_10_19.png", "figs_by_age/admissions_20_29.png", "figs_by_age/admissions_30_39.png", "figs_by_age/admissions_40_49.png", "figs_by_age/admissions_50_59.png", "figs_by_age/admissions_60_69.png", "figs_by_age/admissions_70_79.png", "figs_by_age/admissions_80_plus.png", "figs_by_age/react_5_24.png", "figs_by_age/react_25_34.png", "figs_by_age/react_35_44.png", "figs_by_age/react_45_54.png", "figs_by_age/react_55_64.png", "figs_by_age/react_65_plus.png", "zoomed_view/regions.png", "zoomed_view/prevalence_ons.png", "zoomed_view/prevalence_react.png", "zoomed_view/pillar2_over25.png", "paper_plots/paper_figure_1.png", "paper_plots/paper_figure_2.png", "paper_plots/suppl_age_heatmaps.png", "paper_plots/suppl_compare_demography.png", "paper_plots/suppl_emergence_demography.png", "paper_plots/suppl_compare_hfr.png", "paper_plots/suppl_sev_winter_20_21.png", "paper_plots/suppl_regional_intrinsic.png", "paper_plots/suppl_deaths_hosp_age.png", "paper_plots/suppl_deaths_comm_age.png", "paper_plots/suppl_hosp_adm_age.png", "paper_plots/suppl_inf_prev_age.png", "paper_plots/suppl_pillar2_age.png", "paper_plots/suppl_admissions_vacc.png", "paper_plots/suppl_deaths_vacc.png"))

library(sircovid)
library(spimalot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(jtools)
library(patchwork)
library(reshape2)
library(DescTools)
library(ggrepel)
library(stringr)
library(lubridate)
library(ggtext)
library(car)

orderly2::orderly_resource("support_paper_plots.R")
orderly2::orderly_resource("support_data.R")
orderly2::orderly_resource("supplement_plots.R")
orderly2::orderly_resource("plot_mu_D.R")
source("support_paper_plots.R")
source("support_data.R")
source("supplement_plots.R")
source("plot_mu_D.R")

source("util.R")

version_check("sircovid", "0.14.11")
version_check("spimalot", "0.8.23")

sircovid_model <- "lancelot"

dat <- spimalot::spim_combined_load("regional_results",
                                    regions = "england",
                                    get_severity = TRUE,
                                    get_onward = FALSE)

dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)
dir.create("figs_by_age", FALSE, TRUE)
dir.create("zoomed_view", FALSE, TRUE)
dir.create("thesis_plots", FALSE, TRUE)
dir.create("paper_plots", FALSE, TRUE)

saveRDS(dat$data, "outputs/aggregated_data.rds")

saveRDS(dat$rt$england, "regional_results/Rt_england.rds")
saveRDS(dat$rt, "regional_results/Rt_all.rds")
saveRDS(dat$onward, "outputs/combined.rds")
spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")

model_demography <-
  list(model_demography = dat$model_demography,
       severity_data = dat$parameters$base[[1]]$severity_data)
saveRDS(model_demography, "outputs/model_demography.rds")
write_csv(dat$intrinsic_severity, "outputs/intrinsic_severity.csv")

par_names <- colnames(dat$samples[[1]]$pars)

subset_variants <- intersect(par_names,
                             c("ta_alpha", "rel_p_H_alpha",
                               "rel_p_ICU_alpha", "rel_p_D_alpha",
                               "ta_delta", "rel_p_H_delta",
                               "rel_p_ICU_delta",  "rel_p_D_delta",
                               "ta_omicron", "rel_p_ICU_omicron",
                               "rel_p_H_omicron", "rel_p_D_omicron"))
subset_tv_severity <- intersect(par_names, 
                                c("mu_D", "mu_D_5", "p_H",
                                  "mu_D_2", "p_H_D", "p_H_2",
                                  "mu_D_3", "p_ICU_D", "p_ICU",
                                  "mu_D_4", "p_W_D", "p_ICU_2"))
subset_misc <-   
  setdiff(grep("^beta", par_names, value = TRUE, invert = TRUE),
          c(subset_variants, subset_tv_severity))
seed_dates <- c("start_date", "seed_date_alpha", "seed_date_delta",
                "seed_date_omicron")
subset_misc <- intersect(par_names,
                         c(seed_dates, setdiff(subset_misc, seed_dates)))
par_labels <- forest_plot_labels(dat)

## Add traceplots
for (r in sircovid::regions("england")) {
  fig_name <- paste0("figs/traceplot_", r, ".png")
  write_png(fig_name, width = 3000, height = 1800, res = 200,
            plot_traceplots(dat$samples[[r]]))
}

write_png("figs/forest_plot_variants.png", width = 1600, height = 1600, res = 200,
          spim_plot_forest(
            dat, plot_type = "subset", nrow = length(subset_variants) / 4,
            subset = subset_variants,
            par_labels = par_labels))

write_png("figs/forest_plot_tv_severity.png", width = 1600, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "subset",
                           subset = subset_tv_severity,
                           par_labels = par_labels))

write_png("figs/forest_plot_misc.png",
          width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "subset",
                           subset = subset_misc,
                           par_labels = par_labels))

write_png("figs/forest_plot_betas.png", width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "betas",
                           par_labels = par_labels))


write_png("figs/mu_D.png", width = 2400, height = 1200, res = 200,
          plot_mu_D(
            dat, sircovid::regions("england")))

write_png("figs/data_fits.png", width = 2400 / 5 * 8, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, c(sircovid::regions("england"), "england"),
            c("deaths_hosp", "deaths_comm", "deaths", "icu",
              "general", "hosp", "all_admission"), age_band = "all",
            with_forecast = FALSE, add_betas = FALSE))

write_png("figs/serology_euroimmun.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(
            dat, sircovid::regions("england"), 1, 60))

write_png("figs/serology_roche_n.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(dat, sircovid::regions("england"), 2, 60))

write_png("figs/status_effective_susceptible.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_effective_susceptible(
            dat, sircovid::regions("england"),
            strain_names = c("Wildtype", "Alpha", "Delta", "Omicron"),
            as.Date(c("2019-12-31", "2020-09-18", "2021-03-09", "2021-11-02"))))

write_png("figs/status_vaccine.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_vaccine_status(
            dat, sircovid::regions("england"),
            c("Unvaccinated", "First dose - no effect", "First dose - full effect",
              "Second dose", "Waned 2nd dose", "Booster", "Waned booster")))

write_png("figs/status_infection.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_infection_status(
            dat, sircovid::regions("england")))

write_png("figs/infections_per_strain.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_infections_per_strain(
            dat, sircovid::regions("england"),
            strain_names = c("Wildtype", "Alpha", "Delta", "Omicron"),
            as.Date(c("2019-12-31", "2020-09-18", "2021-03-09", "2021-11-02"))))

write_png("figs/cumulative_attack_rate.png",
          width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_cumulative_attack_rate(
            dat, sircovid::regions("england")))

pillar2_age_bands <-
  c("over25", "under15", "15_24", "25_49", "50_64", "65_79", "80_plus")

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

write_png("figs/prevalence_react.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 10))

write_png("figs/prevalence_ons.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_ons(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 10))

write_png("figs/variant_Wildtype_Alpha.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"), "Alpha",
            date_min = as.Date("2020-09-17"),
            date_max = as.Date("2021-03-01")))

write_png("figs/variant_Alpha_Delta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"), "Delta",
            date_min = as.Date("2021-03-01"),
            date_max = as.Date("2021-08-15")))

write_png("figs/variant_Delta_Omicron.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"), "Omicron",
            date_min = as.Date("2021-11-03"),
            date_max = as.Date("2022-02-01")))

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

## add (zoomed in) plots of main trajectories
write_png("zoomed_view/regions.png", width = 2400 / 5 * 8, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, c(sircovid::regions("england"), "england"),
            c("deaths", "deaths_hosp", "icu", "general",
              "hosp", "all_admission"),
            date_min = as.Date(dat$info$date) - 75, age_band = "all",
            with_forecast = FALSE, add_betas = TRUE))

write_png("zoomed_view/pillar2_over25.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_pillar2_positivity(
            dat, sircovid::regions("england"), "over25",
            date_min = as.Date(dat$info$date) - 75,
            ymax = 50, add_betas = TRUE))

write_png("zoomed_view/prevalence_react.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"),
            date_min = as.Date(dat$info$date) - 75,
            ymax = 10, add_betas = TRUE))

write_png("zoomed_view/prevalence_ons.png", width = 2400, height = 1200,
          res = 200,
          spimalot::spim_plot_ons(
            dat, sircovid::regions("england"),
            date_min = as.Date(dat$info$date) - 75,
            ymax = 10, add_betas = TRUE))

## Plot outputs by age

deaths_age_bands <- c("0_49", "50_54", "55_59", "60_64", "65_69", "70_74",
                      "75_79", "80_plus")
for (i in deaths_age_bands) {
  fig_name <- paste0("figs_by_age/deaths_hosp_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, regions = c(sircovid::regions("england"), "england"),
              "deaths_hosp", age_band = i, with_forecast = FALSE, add_betas = FALSE))
  
  fig_name <- paste0("figs_by_age/deaths_comm_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, regions = c(sircovid::regions("england"), "england"),
              "deaths_comm", age_band = i, with_forecast = FALSE, add_betas = FALSE))
}

admissions_age_bands <- c("0_9", "10_19", "20_29", "30_39", "40_49",
                          "50_59", "60_69", "70_79", "80_plus")
for (i in admissions_age_bands) {
  fig_name <- paste0("figs_by_age/admissions_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, regions = c(sircovid::regions("england"), "england"),
              "all_admission", age_band = i, with_forecast = FALSE, add_betas = FALSE))
}

react_age_bands <- c("5_24", "25_34", "35_44", "45_54", "55_64", "65_plus")
for (i in react_age_bands) {
  fig_name <- paste0("figs_by_age/react_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_react(
              dat, sircovid::regions("england"),
              date_min = as.Date("2020-05-15"), ymax = 10, age_band = i))
}

## Plot main figures for severity paper
# Bespoke data to plot
outcomes_vacc_status <- read_csv("outcomes_vacc_status.csv")
hfr_data <- read_csv("severity_data.csv")
england_data <- read.csv("england_region_data.csv")
vam_data <- get_strain_timelines(dat$data)
hfr_week <- read_csv("hfr_week.csv")
strain_epochs <- sircovid::sircovid_date_as_date(
  dat$parameters$base[[1]]$epoch_dates)[dat$parameters$base[[1]]$strain_epochs[-1L]]
r0_all_regions <- dplyr::bind_rows(
  lapply(sircovid::regions("england"), function(r) get_r0_region(dat, r)))


# Save national intrinsic values for SA post-processing
periods <- as.vector(unique(dat$intrinsic_severity$period))
national_severity <- purrr::map(periods, function(x) get_national_intrinsic_values(dat, x))
names(national_severity) <- periods
saveRDS(national_severity, "outputs/national_severity.rds")

# Paper plots
png("paper_plots/paper_figure_1.png", units = "in", width = 16.5, height = 10, res = 300)
paper_plot_1(dat, "england", vam_data)
dev.off()

png("paper_plots/paper_figure_2.png", units = "in", width = 16.5, height = 10, res = 300)
paper_plot_2(dat, "england", strain_epochs, vam_data, age_bands_select = TRUE)
dev.off()


## Supplement plots
png("paper_plots/suppl_age_heatmaps.png", units = "in", width = 18, height = 15, res = 300)
suppl_age_heatmaps(dat, "england", vam_data)
dev.off()

png("paper_plots/suppl_regional_intrinsic.png", units = "in", width = 11, height = 11, res = 300)
suppl_regional_intrinsic(dat, "england", "Emergence3", strain_epochs, vam_data)
dev.off()

png("paper_plots/suppl_compare_demography.png", units = "in", width = 10, height = 10, res = 300)
suppl_compare_demography(dat, "england")
dev.off()

png("paper_plots/suppl_emergence_demography.png", units = "in", width = 10, height = 10, res = 300)
suppl_emergence_demography(dat, "england")
dev.off()

png("paper_plots/suppl_compare_hfr.png", units = "in", width = 10, height = 10, res = 300)
suppl_compare_hfr(dat, hfr_week, vam_data, strain_epochs)
dev.off()

sev_winter <- suppl_sev_winter_20_21(dat, england_data, hfr_data, strain_epochs, "HFR",
                       as.Date("2020-09-01"), xmax = as.Date("2021-03-01"))
png("paper_plots/suppl_sev_winter_20_21.png", units = "in", width = 18, height = 10, res = 300)
sev_winter$plot
dev.off()



## Add supplementary fits by age plots with national aggregation
write_png("paper_plots/suppl_deaths_hosp_age.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_trajectories_by_age(
            dat, regions = "england", "deaths_hosp", age_band = deaths_age_bands,
            with_forecast = FALSE, add_betas = FALSE,
            title = "Hospital deaths by age - England"))

write_png("paper_plots/suppl_deaths_comm_age.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_trajectories_by_age(
            dat, regions = "england", "deaths_comm", age_band = deaths_age_bands,
            with_forecast = FALSE, add_betas = FALSE,
            title = "Community deaths by age - England"))


write_png("paper_plots/suppl_hosp_adm_age.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_trajectories_by_age(
            dat, regions = "england", "all_admission", age_band = admissions_age_bands,
            with_forecast = FALSE, add_betas = FALSE,
            title = "Hospital admissions by age - England"))

write_png("paper_plots/suppl_inf_prev_age.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, "england", date_min = as.Date("2020-05-15"),
            ymax = 10, age_band = react_age_bands,
            title = "Infection prevalence (REACT-1 study) - England"))

write_png("paper_plots/suppl_pillar2_age.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_pillar2_positivity(
            dat, "england", pillar2_age_bands[-c(1, 2)],
            date_min = as.Date("2020-05-15"), ymax = 60,
            title = "Pillar 2 positivity (new infections) - England"))

png("paper_plots/suppl_admissions_vacc.png", units = "in", width = 14, height = 10, res = 300)
suppl_outcomes_vacc_status(dat, outcomes_vacc_status, "diagnoses_vacc")
dev.off()

png("paper_plots/suppl_deaths_vacc.png", units = "in", width = 14, height = 10, res = 300)
suppl_outcomes_vacc_status(dat, outcomes_vacc_status, "D_hosp_vacc")
dev.off()

## Render rmd
rmarkdown::render("paper_numbers.Rmd")
