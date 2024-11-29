
test_that("get_theta2 is producing fourier functions correctly", {
  ff_params <- list(g0 = 0.2845931,
                    g = c(-0.2992386, -0.0313528, 0.1659284),
                    h = c(-0.35471881, 0.32595719, -0.07911448),
                    theta_c = 0.2845931)
  theta2 <- get_seasonal_forcing(ff_params,n_days = 100)
  expect_equal(round(theta2[10],5), 0.40319)
  expect_true(length(theta2) == 100)
})

test_that("set_seasonality returns a numeric output of n_days",{
  ff_params <- list(g0 = 0.2845931,
                    g = c(-0.2992386, -0.0313528, 0.1659284),
                    h = c(-0.35471881, 0.32595719, -0.07911448),
                    theta_c = 0.2845931)
  params <- get_parameters(n_days = 100) |>
    set_seasonality(season_params = ff_params)

  expect_equal(length(params$daily_rain_input),100)
  expect_equal(sum(is.na(params$daily_rain_input)),0)
})

