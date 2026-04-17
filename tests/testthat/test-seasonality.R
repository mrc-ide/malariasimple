
test_that("get_theta2 is producing fourier functions correctly", {
  g0 <- 0.2845931
  g <- c(-0.2992386, -0.0313528, 0.1659284)
  h <- c(-0.35471881, 0.32595719, -0.07911448)
  theta2 <- get_seasonal_forcing(t = 1:100,
                                 g0 = g0,
                                 g = g,
                                 h = h)
  expect_equal(round(theta2[10],5), 0.40319)
  expect_true(length(theta2) == 100)
})

test_that("set_seasonality returns a numeric output of n_days",{
  n_days <- 100
  g0 <- 0.2845931
  g <- c(-0.2992386, -0.0313528, 0.1659284)
  h <- c(-0.35471881, 0.32595719, -0.07911448)
  params <- get_parameters(n_days = n_days) |>
    set_seasonality(g0 = g0,
                    g = g,
                    h = h)

  expect_equal(length(params$daily_rain_input),n_days+1)
  expect_equal(sum(is.na(params$daily_rain_input)),0)
})

test_that("Error is produced if trying to set equilibrium before seasonality",{
  g0 <- 0.2845931
  g <- c(-0.2992386, -0.0313528, 0.1659284)
  h <- c(-0.35471881, 0.32595719, -0.07911448)
  expect_warning(
    params <- get_parameters() |>
    set_equilibrium(init_EIR = 5) |>
    set_seasonality(g0 = g0,
                    g = g,
                    h = h)
  )
})

test_that("set_rainfall_manual gives an appropriate length output",{
  n_days <- 100
  rain_days <- 200
  cc_df <- data.frame(day = 1:rain_days,
                            cc = 3*sin((1:rain_days * 2*pi) / 365) + 3)
  params <- get_parameters(n_days = n_days) |>
    set_rainfall_manual(cc_df$cc)

  expect_equal(length(params$daily_rain_input), n_days + 1)
})

test_that("set_rainfall_manual produces errors and warnings as expected",{
  n_days <- 100
  rain_days <- 200
  cc_df <- data.frame(day = 1:rain_days,
                      cc = 3*sin((1:rain_days * 2*pi) / 365) + 3)
  params <- get_parameters(n_days = n_days)
  expect_warning(
      set_equilibrium(params, init_EIR = 10) |>
      set_rainfall_manual(cc_df$cc)
  )

  cc_list <- as.list(cc_df$cc)

  expect_error(
    set_rainfall_manual(params, cc_list)
  )

  expect_error(
    set_rainfall_manual(params, 5)
  )

  expect_error(
    set_rainfall_manual(params, cc_df$cc - 5)
  )


})


