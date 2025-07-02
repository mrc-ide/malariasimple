test_that("n_days must be > 0", {
  expect_error(get_parameters(n_days = -1))
})

##age_rendering
test_that("Default values give correct age_rendering outputs",{
  params <- list(age_vector = c(0,0.25,0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8.5,10,20,30,40,60,80)*365)
  age_render <- age_rendering(params = params,
                              prevalence_rendering_max_ages = NULL,
                              prevalence_rendering_min_ages = NULL,
                              clin_inc_rendering_max_ages = NULL,
                              clin_inc_rendering_min_ages = NULL)
  expect_equal(age_render$min_age_prev, 6)
  expect_equal(age_render$max_age_prev, 14)
  expect_equal(age_render$min_age_inc, 1)
  expect_equal(age_render$max_age_inc, 20)
})

##ft
test_that("Zero ft provides zero treated patients",{
  age_vector <- c(0,5,50) * 365
  params <- get_parameters(daily_ft = 0, age_vector = age_vector) |> set_equilibrium(init_EIR = 100)
  sim <- run_simulation(params)
  expect_equal(sum(sim[,"T_count"]),0)
})

test_that("Varying ft in get_parameters gives varying T_count",{
  ft <- seq(0.1,0.9,length.out = 200)
  age_vector <- c(0,5,50) * 365
  params <- get_parameters(age_vector = age_vector,
                           n_days = 100,
                           daily_ft = ft) |>
    set_equilibrium(init_EIR = 20)

  sim <- run_simulation(params)
  expect_lt(sim[1,"T_count"], sim[100,"T_count"])
})
