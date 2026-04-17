itn_days <- c(20)
itn_cov <- c(0.4)

smc_days <- 30
smc_cov <- 0.8

params <- get_parameters() |>
  set_bednets(continuous_distribution = FALSE,
              days = itn_days,
              coverages = itn_cov) |>
  set_smc(coverages = smc_cov,
          days = smc_days,
          min_age = 0,
          max_age = 5*365,
          distribution_type = "correlated") |>
  set_equilibrium(init_EIR = 20)

sim <- run_simulation(params, full_output = TRUE)



