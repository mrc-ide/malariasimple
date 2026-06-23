test_that("get_custom_output() function is working for standard deterministic usecase", {
  n_days <- 70
  params <- get_parameters(
    n_days = n_days,
    age_vector = c(0, 0.25, 0.5, 1, 3, 5, 10, 50) *
                             365) |>
    set_equilibrium(init_EIR = 10)
  sim <- run_simulation(params, full_output = TRUE)
  custom <- get_custom_output(
    params,
    sim,
    output_variables = c("n", "detect"),
    ages = c(0, 0.5) * 365
  )
  expect_true(all(dim(custom) == c(n_days,3)))
})

test_that("get_custom_output() function is working for standard stochastic usecase", {
  n_days <- 70
  n_particles <- 2
  params <- get_parameters(
    stochastic = TRUE,
    n_days = n_days,
    age_vector = c(0, 0.25, 0.5, 1, 3, 5, 10, 50) *
                             365) |>
    set_equilibrium(init_EIR = 10)
  sim <- run_simulation(params, n_particles = n_particles, full_output = TRUE)
  custom <- get_custom_output(
    params,
    sim,
    output_variables = c("n", "detect"),
    ages = c(0, 0.5) * 365
  )
  expect_true(all(dim(custom) == c(n_days * n_particles, 4)))
})

test_that("make 2D works as expected", {
  n_particles <- 3
  n_days <- 30
  params <- get_parameters(stochastic = TRUE, n_days = n_days) |> set_equilibrium(init_EIR = 10)
  sim <- run_simulation(params, n_particles = n_particles)
  expect_true(dim(sim)[1] == 30)
  expect_true(dim(sim)[3] == 3)

  sim_flat <- make_2d(sim)
  expect_true(length(dim(sim_flat)) == 2)
  expect_true(nrow(sim_flat) == n_particles * n_days)
})


