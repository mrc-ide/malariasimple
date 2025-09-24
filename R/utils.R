#' @title Timing of peak carrying capacity
#' @description Returns the day-of-year on which vector carrying capacity is maximum
#' @param g0 Mean baseline coefficient
#' @param g Cosine coefficients
#' @param h Sine coefficients
#' @examples
#' #Seasonality parameters
#' g0 = 0.28
#' g = c(-0.3, -0.03, 0.17)
#' h = c(-0.35, 0.32, -0.07)
#'
#' peak_cc <- get_peak_cc(g0, g, h)
#' @export

get_peak_cc <- function(g0, g, h){
  seasonal_forcing <- get_seasonal_forcing(t = 1:365, g0, g, h)
  peak_day <- which.max(seasonal_forcing)
  return(peak_day)
}

#' @title Decay matrix for ITN or SMC
#' @description Produces an ij matrix representing the decay in efficacy of insecticide-treated nets (ITN) or
#' prophylactic protection (SMC) on day j for individuals in receipt of distribution event i.
#' @param days Days on which distribution events occur
#' @param n_days Length of simulation (days)
#' @param gamman_itn Half-life of ITN insecticide (days)
#' @param scale_smc Scale parameter of Weibull distribution defining decay of prophylactic protection of SMC
#' @param shape_smc Shape parameter of Weibull distribution defining decay of prophylactic protection of SMC
#' @param intervention "ITN" or "SMC"
#'
#' @export
##Intervention coverage
get_decay_mat <- function(days,n_days,gamman_itn = NULL,scale_smc = NULL,shape_smc = NULL,intervention = NULL){
  n_dists <- length(days) #Number of ITN/SMC distribution events
  decay_mat <- matrix(0,nrow = n_days, ncol = n_dists+1)

  for(i in 1:n_dists){
    dist_day <- days[i] #Days on which distribution events occur
    t <- 1:(n_days - dist_day + 1)
    if(intervention == "SMC"){
      decay_mat[dist_day:n_days,i] <- exp(-((t / scale_smc) ^ shape_smc))
    } else if(intervention == "ITN"){
      #decay_rate <- log(2) / gamman_itn[i] #This is what it should be?
      decay_rate <- 1 / gamman_itn[i] #This matches malariasimulation? https://github.com/mrc-ide/malariasimulation/blob/master/R/vector_control.R bednet_decay()
      decay_mat[dist_day:n_days,i] <- exp(-t*decay_rate)
    }
  }
  return(decay_mat)
}

#' @title Daily Effective Decay
#' @description Estimates the mean average intervention decay across all distribution events
#' @param usage_mat Matrix of usage values for each distribution event and simulation day.
#' @param decay_mat Matrix of decay values for each distribution event and simulation day.

#' @export
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

