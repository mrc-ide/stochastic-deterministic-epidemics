orderly2::orderly_resource(
  c("data/data_vaccination.csv",
    "data/uk_rtm.csv",
    "data/serology.csv",
    "data/weighted_prior_ranges.csv"))

orderly2::orderly_artefact(
  "Data for use in fits",
  c("data/data_vaccination.csv",
    "data/uk_rtm.csv",
    "data/serology.csv",
    "data/weighted_prior_ranges.csv"))
