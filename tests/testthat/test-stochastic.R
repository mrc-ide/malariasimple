test_that("Detected malaria prevalence provides integer outputs", {
  params <- get_parameters(stochastic = TRUE,
                           human_pop = 1000,
                           prevalence_rendering_min_ages = 730,
                           prevalence_rendering_max_ages = 3650) |>
    set_equilibrium(init_EIR = 10)
  out <- run_simulation(params)
  detect <- out[,"n_detect_730_3650"]
  expect_true(all(detect == as.integer(detect)))
})
