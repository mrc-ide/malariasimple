#' Set ITN parameters
#' @export
#' @param params malariasimple parameters
#' @param continuous_distribution Is ITN distribution continuous? If FALSE, distribution is assumed to occur in discrete events.
#' @param daily_continuous_cov Vector of daily ITN coverage (required when continuous_distribution = TRUE). A single value is also accepted
#' @param days Vector of days on which ITN distribution events occur (required when continuous_distribution = FALSE)
#' @param coverages Vector detailing the proportion of the population receiving an ITN during each intervention (required when continuous_distribution = FALSE)
#' @param gamman Half-life of ITN insecticide (days)
#' @param retention Average number of days a net is kept for
#' @param dn0 Probability of mosquito dying upon an encounter with ITN (max)
#' @param rn Probability of repeating behaviour with ITN (max)
#' @param rnm Probability of repeating behaviour with ITN (min)
#'
#' @returns Updates the input parameter list to include ITN parameters
#'
#' @examples
#' n_days <- 500
#' #Discrete distribution scenario
#' discrete_itn_params <- get_parameters(n_days = n_days) |>
#'   set_bednets(days = c(50,100,200),
#'           coverages = c(0.2,0.5,0.1)) |>
#'   set_equilibrium(init_EIR = 10)
#'
#' #Continuous distribution scenario
#' coverageserage <- 0.2 + 0.15*(sin(2 * pi * (1:n_days / 365)) + 1)
#' discrete_itn_params <- get_parameters(n_days = n_days) |>
#'   set_bednets(continuous_distribution = TRUE,
#'           daily_continuous_cov = coverageserage) |>
#'   set_equilibrium(init_EIR = 10)
set_bednets <- function(params,
                    continuous_distribution = FALSE, #Assume discrete distributions or continuous
                    daily_continuous_cov = NULL,
                    days = NULL,
                    coverages = NULL,
                    gamman =   2.64 * 365,
                    retention =   3 * 365,
                    dn0 = 0.41,
                    rn = 0.56,
                    rnm = 0.24){
  #----------------------------------- Sanity Checks -----------------------------------------------
  if (params$equilibrium_set == 1) warning(message("Equilbrium must be set last"))

  if (!continuous_distribution){
    if (max(coverages) > 1) stop(message("coverages cannot exceed 1"))
    if(is.null(days)) stop(message("days must be specific when continuous_distribution == FALSE"))
    if(is.null(coverages)) stop(message("coverages must be specific when continuous_distribution == FALSE"))
    if(max(days) > params$n_days){
      days <- days[days >= n_days]
      coverages <- coverages[days >= n_days]
    }
    if (length(days) != length(coverages)){
      if(length(coverages) != 1){
        stop(message("days and coverages must be equal in length"))
      } else {
        coverages <- rep(coverages,length(days))
      }}
  }
  if (continuous_distribution){
    if(is.null(daily_continuous_cov)) stop(message("daily_continuous_cov must be specified when continuous_distribution == TRUE"))
    if(length(daily_continuous_cov) < params$n_days) stop(message("n_days cannot exceed length(daily_continuous_cov)"))
    if(max(daily_continuous_cov) > 1) stop(message("daily_continuous_cov cannot exceed 1"))
    daily_continuous_cov <- daily_continuous_cov[1:params$n_days]
  }
  #----------------------------------- Set Parameters -----------------------------------------------
  params$dn0 = dn0
  params$rn = rn
  params$rnm = rnm

  if(continuous_distribution){
    params <- itn_continuous_distribution_params(params,
                                                 daily_continuous_cov = daily_continuous_cov,
                                                 gamman = gamman,
                                                 retention = retention)
  }

  if(!continuous_distribution){
    params <- itn_discrete_distribution_params(params = params,
                                               days = days,
                                               coverages = coverages,
                                               gamman = gamman,
                                               retention = retention)
  }
  params$itn_set <- 1
  params$equilibrium_set <- 0
  return(params)
}

#' @title Produce matrix of population ITN usage
#' @description Produces an ij matrix defining the proportion of the population using an ITN from distribution i on day j. Final column represents population with no ITN.
#' @param days Vector of days on which ITN distribution events occur (required when continuous_distribution = FALSE)
#' @param coverages Vector detailing the proportion of the population receiving an ITN during each intervention (required when continuous_distribution = FALSE)
#' @param retention Average number of days a net is kept for
#' @param n_days Length of simulation
#' @export
get_itn_usage_mat <- function(days,coverages,retention,n_days){
  n_dists <- length(coverages) #Number of ITN distribution events

  ##Matrix of the population prop. using a net from distribution event i (col number) on day t (row number)
  #Final column is those with no bednet
  itn_mat <- matrix(nrow = n_days, ncol = n_dists+1)
  itn_mat[1,] <- c(rep(0,n_dists),1) #Initialise no bednets
  for(i in 2:n_days){
    #People discard nets
    itn_mat[i,1:n_dists] <- itn_mat[(i-1),(1:n_dists)]*exp(-1/retention)
    itn_mat[i,(n_dists+1)] <- 1 - sum(itn_mat[i,(1:n_dists)])

    #New nets are distributed on days
    if(i %in% days){
      dist_index <- which(days == i) #Which distribution event
      new_nets <- coverages[dist_index]
      itn_mat[i,dist_index] <- new_nets
      itn_mat[i,-dist_index] <- itn_mat[i,-dist_index] - new_nets*itn_mat[i,-dist_index]
    }
  }
  return(itn_mat)
}

itn_discrete_distribution_params <- function(params,
                                             days,
                                             coverages,
                                             gamman, #ITN half life
                                             retention #Average number of days a net is kept for
){
  #Matrix of the proportion of individuals currently using a net from each distribution event
  usage_mat <- get_itn_usage_mat(days = days, coverages = coverages,
                                 retention = retention, n_days = params$n_days)

  #Overall itn coverage
  daily_coverages <- get_daily_cov(usage_mat)

  #The level of insecticide decay of the average net from each itn distribution event
  decay_mat <- get_decay_mat(days = days, gamman_itn = gamman,
                             n_days = params$n_days, intervention = "ITN")

  #Decay level averaged over all itns currently in use
  if (length(days) == 1){
    daily_itn_decay <- usage_mat[ ,1] * decay_mat[ ,1]
  } else {
    daily_itn_decay <- get_daily_decay(usage_mat, decay_mat)
  }

  #Update parameters
  params$itn_decay_daily <- c(0,daily_itn_decay)
  params$max_itn_cov <- max(daily_coverages)
  if(params$max_itn_cov == 0){
    params$itn_eff_cov_daily <- rep(0,params$n_days+1)
  } else {
    params$itn_eff_cov_daily <- c(0,daily_coverages / params$max_itn_cov)
  }
  return(params)
}

itn_continuous_distribution_params <- function(params,daily_continuous_cov, gamman, retention){
  params$max_itn_cov <- max(daily_continuous_cov)
  params$itn_eff_cov_daily <- c(0,daily_continuous_cov/max(daily_continuous_cov))
  decay_rate <- log(2) / gamman
  mean_itn_decay <- mean(exp(-((1:retention)) * decay_rate))
  params$itn_decay_daily <- rep(mean_itn_decay, params$n_days + 1)
  return(params)
}

