deterministic <- TRUE
short_run <- TRUE

# assumptions can also be any of:
# crim_infect_high, crim_infect_low
# crim_hospi_high, crim_hospi_low
# crim_death_high, crim_death_low
# booster_ve_high, booster_ve_low
# alpha_ve_high, alpha_ve_low
# delta_ve_high, delta_ve_low
# mu_d_winter, mu_d_summer
# fixed_si_high, fixed_si_low
assumptions <- "central"

## 1. sarscov2_data
orderly2::orderly_run("sarscov2_data")

## 2. sarscov2_parameters 
orderly2::orderly_run(
  "sarscov2_parameters",
  parameters = list(deterministic = deterministic, assumptions = assumptions))

## 3. sarscov2_fits
for (r in sircovid::regions("england")) {
  orderly2::orderly_run(
    "sarscov2_fits",
    parameters = list(region = r,
                      short_run = short_run,
                      deterministic = deterministic,
                      assumptions = assumptions))
}


## 4. sarscov2_fits_combined
orderly2::orderly_run(
  "sarscov2_fits_combined",
  parameters = list(short_run = short_run,
                    deterministic = deterministic,
                    assumptions = assumptions))
