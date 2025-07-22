#' @title Add seasonal patterns
#' @description Adds seasonal forcing to the model by approximating average seasonal rainfall as a Fourier series.
#' @param params Other malariasimple parameters
#' @param g0 Mean baseline coefficient
#' @param g Cosine coefficients
#' @param h Sine coefficients
#' @param floor Minimum value of seasonal forcing
#' @return Updates the input parameter list to include seasonal parameters
#' @examples
#' # Define seasonality Fourier coefficients
#' g0 = 0.28
#' g = c(-0.3, -0.03, 0.17)
#' h = c(-0.35, 0.33, -0.08)
#' params <- get_parameters() |>
#'           set_seasonality(g0=g0,
#'           g=g,
#'           h=h,
#'           floor = 0.005) |>
#'           set_equilibrium(init_EIR = 5)
#'@export
set_seasonality <- function(params, g0, g, h, floor = 0.001){
  if (params$equilibrium_set == 1) warning(message("Equilbrium must be set last"))
  rain_input <- get_seasonal_forcing(t=1:params$n_days, g0, g, h, floor = floor)
  params$daily_rain_input <- c(1,rain_input)
  params$seasonality_set <- 1
  return(params)
}

#' @title Add manual rainfall forcing
#' @description Takes daily rain-forcing time series input and adds it to the parameter set
#' @param params Other malariasimple parameters
#' @param cc_ts Carrying capacity time series. Vector of daily carrying capacity. Length must equal or exceed params$n_days
#' @return Updates the input parameter list to include seasonal parameters
#' @examples
#' n_days = 1000
#' t <- 1:n_days
#' rainfall_ts <- sin((t*2*pi)/365) + 1.1
#' params <- get_parameters() |>
#'           set_rainfall_manual(rainfall_ts) |>
#'           set_equilibrium(init_EIR = 5)
#'@export
#'
set_rainfall_manual <- function(params, cc_ts){
  if (params$equilibrium_set == 1) warning(message("Equilbrium must be set last"))
  if (!is.vector(cc_ts) | is.list(cc_ts)) stop(message("cc_ts must be a vector"))
  if (min(cc_ts) < 0) stop(message("cc_ts contains negative values. Vector carrying capacity cannot be negative."))
  if (length(cc_ts) < params$n_days) stop(message("cc_ts must be at least as long as n_days"))
  if (0 %in% cc_ts) warning(message("cc_ts contains zero values. This may cause the model to behave strangely."))
  params$daily_rain_input <- c(1,cc_ts[1:params$n_days])
  params$seasonality_set <- 1
  return(params)
}


#' @title Get daily seasonal forcing
#' @description Convert Fourier coefficients into a smooth vector of daily seasonal forcing. Used within set_seasonality function.
#' @param t Day-of-year
#' @param g0 Mean baseline coefficient
#' @param g Cosine coefficients
#' @param h Sine coefficients
#' @param floor Minimum permitted value of output
#' @export
get_seasonal_forcing <- function(t, g0, g, h, floor = 0.001) {
  if (floor <= 0) stop(message("floor must be greater than 0"))
  result <- g0
  for (i in seq_along(g)) {
    result <- result +
      g[i] * cos(2 * pi * t * (i) / 365) +
      h[i] * sin(2 * pi * t * (i) / 365)
  }

  #Estimate annual mean
  result_annual <- g0
  for (i in seq_along(g)) {
    result_annual <- result_annual +
      g[i] * cos(2 * pi * 1:365 * (i) / 365) +
      h[i] * sin(2 * pi * 1:365 * (i) / 365)
  }

  result_norm <- result / mean(result_annual)
  result_norm[result_norm < floor] <- floor
  return(result_norm)
}



