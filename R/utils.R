#' @title Timing of peak carrying capacity
#' @description Returns the day-of-year on which carrying capacity is maximum
#' @param season_params List of Fourier function parameters describing seasonality
#' @examples
#' season_params <- list(g0 = 0.28,
#'                           g = c(-0.30, -0.03, 0.17),
#'                           h = c(-0.35, 0.32, -0.08),
#'                           theta_c = 0.28)
#' peak_cc <- get_peak_cc(season_params)

get_peak_cc <- function(season_params){
  seasonal_forcing <- get_seasonal_forcing(season_params, n_days = 365)
  peak_day <- which.max(seasonal_forcing)
  return(peak_day)
}

#' @title Generate malariasimple deterministic model
#' @examples
#' gen <- gen_model()
#' params <- get_parameters() |> set_equilibrium(init_EIR = 5)
#' out <- run_simulation(gen,params)
#'

gen_model <- function(){
  model_path <- system.file("malariasimple_deterministic.R", package = "malariasimple")
  gen <- odin.dust::odin_dust(model_path)
}


