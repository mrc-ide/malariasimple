#' @title Timing of peak carrying capacity
#' @description Returns the day-of-year on which carrying capacity is maximum
#' @param season_params List of Fourier function parameters describing seasonality
#' @examples
#' season_params <- list(g0 = 0.28,
#'                           g = c(-0.30, -0.03, 0.17),
#'                           h = c(-0.35, 0.32, -0.08),
#'                           theta_c = 0.28)
#' peak_cc <- get_peak_cc(season_params)
#' @export

get_peak_cc <- function(season_params){
  seasonal_forcing <- get_seasonal_forcing(season_params, n_days = 365)
  peak_day <- which.max(seasonal_forcing)
  return(peak_day)
}

##Intervention coverage
get_decay_mat <- function(days,n_days,gamman_itn = NULL,scale_smc = NULL,shape_smc = NULL,intervention = NULL){
  n_dists <- length(days) #Number of smc distribution events
  decay_mat <- matrix(0,nrow = n_days, ncol = n_dists+1)

  for(i in 1:n_dists){
    dist_day <- days[i]
    t <- 1:(n_days - dist_day + 1)
    if(intervention == "SMC"){
      decay_mat[dist_day:n_days,i] <- exp(-((t / scale_smc) ^ shape_smc))
    } else if(intervention == "ITN"){
      decay_mat[dist_day:n_days,i] <- exp(-t/gamman_itn)
    }
  }
  return(decay_mat)
}

get_daily_decay <- function(usage_mat, decay_mat){
  n_dists <- ncol(usage_mat) - 1
  prop_cat <- usage_mat[,1:n_dists] / rowSums(usage_mat[,1:n_dists]) #Proportion of pop. in each distribution category
  prop_cat[is.nan(prop_cat)] <- 0 #Dividing by zero causes NaNs
  mean_decay <- rowSums(prop_cat * decay_mat[,1:n_dists])
  return(mean_decay)
}
#Convert intervention usage matrix into daily overall coverage time series
get_daily_cov <- function(usage_mat){
  n_dists <- ncol(usage_mat) - 1
  if(n_dists > 1){
    out <- rowSums(usage_mat[,1:n_dists])
  } else {
    out <- usage_mat[ ,1]
  }
  return(out)
}

