#' @title Add seasonal patterns
#' @description Adds seasonal forcing to the model using a Fourier function of rainfall.
#' @param params Other malariasimple parameters
#' @param admin_unit Name of admin 2 region
#' @param country Name of admin 2 country
#' @param season_params List of Fourier parameters for manual seasonal input
#' @param floor Minimum value of seasonal forcing
#' @return Updates the input parameter list to include seasonal parameters
#' @examples
#' season_params <- list(g0 = 0.28,
#' g = c(0.2, -0.07, -0.001),
#' h = c(0.2, -0.07, -0.1),
#' theta_c = 1)
#' params <- get_parameters() |>
#'           set_seasonality(season_params = season_params, floor = 0.005) |>
#'           set_equilibrium(init_EIR = 5)
#'
#'
set_seasonality <- function(params, season_params, floor = 0.001){
  #if(!is.null(params$daily_rain_input)) warning("Seasonality profile has replaced daily_rain input")
  if (floor <= 0) stop(message("floor must be greater than 0"))
  if (params$equilibrium_set == 1) warning(message("Equilbrium must be set last"))

  theta2 <- get_seasonal_forcing(season_params, params$n_days, floor = floor)
  params$daily_rain_input <- theta2
  params$seasonality_set <- 1
  return(params)
}

#' @title Get seasonal forcing
#' @description Convert Fourier function parameters into a smooth vector of daily seasonal forcing. Used within set_seasonality function
#' @param season_params List of Fourier function parameters describing seasonality
#' @param n_days Numerical value of output vector length
#' @param floor Minimum permitted value of output
#'
#' @examples
#' season_params <- list(g0 = 0.28,
#'                           g = c(-0.30, -0.03, 0.17),
#'                           h = c(-0.35, 0.32, -0.08),
#'                           theta_c = 0.28)
#' get_seasonal_forcing(season_params,
#'                       n_days = 1000,
#'                       floor = 0.005)
#'
get_seasonal_forcing <- function(season_params, n_days, floor = 0.001){
  #Ensure parameters have been correctly described
  if (length(season_params) != 4 |
      sum(names(season_params) %in% c("g0","g","h","theta_c")) != 4 |
      length(season_params$g) != 3 |
      length(season_params$h) != 3){
    stop(message("season_params improperly described. Check ?get_seasonal_forcing."))
  }
  ssa0 <- season_params$g0 |> as.numeric()
  ssa1 <- season_params$g[1] |> as.numeric()
  ssa2 <- season_params$g[2] |> as.numeric()
  ssa3 <- season_params$g[3] |> as.numeric()
  ssb1 <- season_params$h[1] |> as.numeric()
  ssb2 <- season_params$h[2] |> as.numeric()
  ssb3 <- season_params$h[3] |> as.numeric()
  theta_c <- season_params$theta_c |> as.numeric()
  time <- 1:n_days
  theta2 <- (ssa0+
               ssa1*cos(2*pi*time/365)+
               ssa2*cos(2*2*pi*time/365)+
               ssa3*cos(3*2*pi*time/365)+
               ssb1*sin(2*pi*time/365)+
               ssb2*sin(2*2*pi*time/365)+
               ssb3*sin(3*2*pi*time/365) ) /theta_c
  theta2[theta2 < 0.001] <- floor #Set a floor of 0.001
  return(theta2)
}

