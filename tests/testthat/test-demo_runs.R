test_that("Deterministic demo run produces incidence and prevalence as expected",{
  n_days <- 365
  init_EIR <- 10

  Thies_ff <- list(
  g0 = 0.28,
  g = c(-0.3, -0.03, 0.17),
  h = c(-0.35, 0.33, -0.08),
  theta_c = 0.28)

  params <- get_parameters(
    n_days = n_days,
    stochastic = FALSE,
    human_pop = 1,
    tsd = 3,
    het_brackets = 5,
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365,
    clin_inc_rendering_min_ages = 0,
    clin_inc_rendering_max_ages = Inf,
    age = c(0, 1, 2, 3, 4, 5, 10, 20, 40, 60, 80) *
      365
  ) |>
    set_seasonality(Thies_ff) |>
    set_equilibrium(init_EIR = init_EIR)

  out <- run_simulation(params) |> as.data.frame()

  expect_equal(round(out$n_detect_730_3650[300], 3), 0.134)
  expect_equal(round(out$n_clin_inc_0_Inf[10], 5), 0.00144)
})


