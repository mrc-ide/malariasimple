test_that("get_itn_usage_mat produces zero usage from zero coverage", {
  itn_usage_mat <- get_itn_usage_mat(days = c(1,50,99), coverages = c(0,0,0),
                                     retention = 100, n_days = 200)
  expect_true(all(itn_usage_mat[,4] == 1))
  expect_true(all(itn_usage_mat[,1:3] == 0))
})

test_that("daily outputs are all the same length (continuous_distribution == TRUE)", {
  cov_days <- 500
  n_days <- 300
  itn_coverage <- 0.2 + 0.15*(sin(2 * pi * (1:cov_days / 365)) + 1)
  params <- get_parameters(n_days = n_days) |>
    set_bednets(continuous_distribution = TRUE,
            daily_continuous_cov = itn_coverage) |>
    set_equilibrium(init_EIR = 10)

  expect_equal(length(params$itn_decay_daily),n_days+1)
  expect_equal(length(params$itn_eff_cov_daily),n_days+1)

  cov_days <- 300
  n_days <- 500
  itn_coverage <- 0.2 + 0.15*(sin(2 * pi * (1:cov_days / 365)) + 1)
  expect_error(suppressMessages(get_parameters(n_days = n_days) |>
                                  set_bednets(continuous_distribution = TRUE,
                                          daily_continuous_cov = itn_coverage)))
})

test_that("daily outputs are all the same length (continuous_distribution == FALSE)", {
  n_days <- 300
  itn_days <- c(4,n_days)
  itn_cov <- c(0.2,0.4)
  params <- get_parameters(n_days = n_days) |>
    set_bednets(continuous_distribution = FALSE,
            days = itn_days,
            coverages = itn_cov) |>
    set_equilibrium(init_EIR = 10)

  expect_equal(length(params$itn_decay_daily),n_days+1)
  expect_equal(length(params$itn_eff_cov_daily),n_days+1)

  n_days <- 300
  itn_days <- c(4,500)
  itn_cov <- c(0.2,0.4)

  expect_error(suppressMessages(get_parameters(n_days = n_days) |>
                                  set_bednets(continuous_distribution = FALSE,
                                          days = itn_days,
                                          coverages = itn_cov)))
})

test_that("Zero valued itn coverage does not result in NaN values",{
  params <- get_parameters(n_days = 100) |>
    set_bednets(days = 50, coverages = 0) |>
    set_equilibrium(init_EIR = 5)
  expect_equal(sum(params$itn_eff_cov_daily), 0)
})

