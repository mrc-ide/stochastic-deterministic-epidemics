orderly2::orderly_resource(c("outputs/data_vaccination.csv", "outputs/england_region_data.csv", "outputs/hfr_week.csv", "outputs/outcomes_vacc_status.csv", "outputs/serology_for_inference.csv", "outputs/severity_data.csv", "outputs/weighted_prior_ranges.csv"))

orderly2::orderly_artefact(
  "Data for use in fits",
  c("outputs/data_vaccination.csv", "outputs/england_region_data.csv", "outputs/hfr_week.csv", "outputs/outcomes_vacc_status.csv", "outputs/serology_for_inference.csv", "outputs/severity_data.csv", "outputs/weighted_prior_ranges.csv"))
