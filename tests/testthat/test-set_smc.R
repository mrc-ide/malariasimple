test_that("set_drug_params produces a sensible output for SP_AQ", {
  drug_params <- set_drug_params("SP_AQ")
  expect_equal(length(drug_params),4)
  expect_true(is.numeric(drug_params))
})

test_that("get_smc_usage_mat returns zero usage from zero coverage input",{
  smc_usage_mat_random <- get_smc_usage_mat(
    days = c(10, 100, 200),
    coverages = c(0, 0, 0),
    n_days = 300,
    distribution_type = "random"
  )
  smc_usage_mat_cohort <- get_smc_usage_mat(
    days = c(10, 100, 200),
    coverages = c(0, 0, 0),
    n_days = 300,
    distribution_type = "cohort"
  )

  expect_true(all(smc_usage_mat_cohort[,4] == 1))
  expect_true(all(smc_usage_mat_random[,4] == 1))
  expect_true(all(smc_usage_mat_cohort[,1:3] == 0))
  expect_true(all(smc_usage_mat_random[,1:3] == 0))
})

test_that("alpha_smc is zero when drug_efficacy is zero",{
  params <- get_parameters(n_days = 200) |>
    set_smc(days = c(1,50,99), coverages = 0.5, drug_efficacy = 0)
  expect_true(all(params$alpha_smc == 0))
})

test_that("rel_c_days works at edge cases", {
  days <- c(1,50,99)
  cov <- 0.5
  params_1 <- get_parameters(n_days = 200) |>
    set_smc(days = days, coverages = cov, drug_rel_c = 1)
  expect_true(all(params_1$rel_c_days == 1))

  params_0 <- get_parameters(n_days = 200) |>
    set_smc(days = days, coverages = cov, drug_rel_c = 0)
  x <- 1 - cov / params_0$max_smc_cov #max_smc_cov*x + (1-max_smc_cov)*1 == cov

  expect_true(all(params_0$rel_c_days[days+1] == x))
  expect_true(all(params_0$rel_c_days %in% c(1,x)))
})

test_that("SMC clearance is working as expected", {
  set.seed(1)
  params <- get_parameters(n_days = 100,
                           prevalence_rendering_max_ages = 5*365,
                           prevalence_rendering_min_ages = 0.25*365,
                           stochastic = TRUE,
                           human_pop = 10000) |>
    set_smc(days = c(40,50),
            coverages = c(1,1),
            min_age = 0.25*365,
            max_age = 5*365,
            drug_efficacy = 1,
            smc_clearance_lag = 5) |>
    set_equilibrium(init_EIR = 10)
  out <- run_simulation(params)
  expect_equal(as.numeric(out[47,"n_detect_91.25_1825"]), 0)
})

test_that("Nonsense inputs produce errors",{
  params <- get_parameters()
  expect_error(
    set_smc(
    params = params,
    coverages = c(1, 1),
    min_age = 0.25 * 365,
    max_age = 5 * 365,
    drug_efficacy = 1,
    smc_clearance_lag = 5
  ))
})

test_that("Random SMC distribution is more effective than correlated",{
  params <- get_parameters(n_days = 100)
  corr_params <- set_smc(params,
                         days = c(20, 50),
                         min_age = 2*365,
                         max_age = 10*365,
                         coverages = c(0.5, 0.5),
                         distribution_type = "correlated") |>
    set_equilibrium(init_EIR = 50)
  corr_sim <- run_simulation(corr_params) |> as.data.frame()
  rand_params <- set_smc(params,
                         days = c(20, 50),
                         coverages = c(0.5, 0.5),
                         min_age = 2*365,
                         max_age = 10*365,
                         distribution_type = "random") |>
    set_equilibrium(init_EIR = 50)
  rand_sim <- run_simulation(rand_params) |> as.data.frame()

  expect_equal(corr_sim[10, "n_detect_730_3650"], rand_sim[10, "n_detect_730_3650"])
  expect_true(corr_sim[60, "n_detect_730_3650"] > rand_sim[60, "n_detect_730_3650"])
})

test_that("set_smc sets coverage properly for a single SMC dosage", {
  smc_days <- 20
  smc_cov <- 0.7
  params <- get_parameters(n = 50) |>
    set_smc(coverages = smc_cov,
            days = smc_days) |>
    set_equilibrium(init_EIR = 10)
  expect_equal(params$max_smc_cov, smc_cov)
})

test_that("Nonsense parameters produce errors", {
  params <- get_parameters(n_days = 100)
  expect_error(set_smc(params, days = 5))
  expect_error(set_smc(params, coverages = 0.5))
  expect_error(set_smc(params, days = 6,
                       coverages = 1.5))
  expect_error(set_smc(params, days = 0,
                       coverages = 1.5))
  expect_error(set_smc(params,
                       days = 5,
                       coverages = 0.5,
                       drug = "party"))
})
