#' @title Set IRS parameters
#' @export
#' @param params malariasimple parameters
#' @param days Days on which IRS is implemented
#' @param coverages Population coverage for each day
#' @param distribution_type "random" or "cohort" distribution of IRS
#' @param ls_theta Initial proportion of mosquitoes dying following entering a hut
#' @param ls_gamma Defines how ls changes over time
#' @param ks_theta Initial proportion of mosquitoes successfully feeding
#' @param ks_gamma Defines how ks changes over time
#' @param ms_theta Initial proportion of mosquitoes repelled
#' @param ms_gamma Defines how ms changes over time
set_spraying <- function(params,
                    days,
                    coverages,
                    distribution_type = "random",
                    ls_theta = 2.025,
                    ls_gamma = -0.009,
                    ks_theta = -2.222,
                    ks_gamma = -1.232,
                    ms_theta = -1.232,
                    ms_gamma = -0.009) {
  params$irs_set <- 1
  n_days <- params$n_days
  irs_usage_mat <- get_irs_usage_mat(days, coverages, n_days, distribution_type)
  n_dists <- length(days)
  ls_mat <- matrix(NA, nrow = n_days, ncol = n_dists + 1) #Proportion of mosquitoes dying from IRS following entering a hut
  ks_mat <- matrix(NA, nrow = n_days, ncol = n_dists + 1) #Proportion of mosquitoes successfully feeding
  ms_mat <- matrix(NA, nrow = n_days, ncol = n_dists + 1) #Proportion of mosquitoes deterred
  for (i in 1:n_dists) {
    dist_day <- days[i] #Days on which distribution events occur
    t_since <- 1:(n_days - dist_day + 1)
    ls_mat[dist_day:n_days, i] <- spraying_decay(t_since, ls_theta, ls_gamma)
    ks_mat[dist_day:n_days, i] <- params$k0 * spraying_decay(t_since, ks_theta, ks_gamma)
    ms_mat[dist_day:n_days, i] <- spraying_decay(t_since, ms_theta, ms_gamma)
  }
  ones <- matrix(1, nrow = n_days, ncol = n_dists + 1)
  js_mat <- ones - ls_mat - ks_mat #Proportion of mosquitoes that enter the hut but are repelled without being killed or feeding
  ms_comp_mat <- ones - ms_mat
  ls_prime_mat <- ls_mat * ms_comp_mat
  ks_prime_mat <- ks_mat * ms_comp_mat
  js_prime_mat <- js_mat * ms_comp_mat + ms_mat

  rs_mat <- prob_spraying_repels(ls_prime_mat, ks_prime_mat, js_prime_mat, params$k0)
  rs_mat[is.na(rs_mat)] <- 0

  ss_mat <- prob_survives_spraying(ks_prime_mat, params$k0)
  ss_mat[is.na(ss_mat)] <- 1

  #Decay level averaged over all itns currently in use
  if (n_dists == 1){
    ss <- ss_mat[ ,1]
    rs <- rs_mat[ ,1]
  } else {
    ss <- get_daily_irs(irs_usage_mat, ss_mat)
    rs <- get_daily_irs(irs_usage_mat, rs_mat)
  }
  ss[is.nan(ss)] <- 1
  rs[is.nan(rs)] <- 0

  daily_irs_cov <- get_daily_cov(irs_usage_mat)
  params$max_irs_cov <- max(daily_irs_cov)

  irs_eff_cov_daily <- daily_irs_cov / params$max_irs_cov
  params$irs_eff_cov_daily <- c(0,irs_eff_cov_daily)
  params$ss <- c(1,ss)
  params$rs <- c(0,rs)
  return(params)
}

#' @title Produces matrix of IRS usage
#' @export
#' @description Produces an ij matrix defining the proportion of the population using an IRS from distribution i on day j. Final column represents population with no IRS.
#' @param days Days on which IRS is implemented
#' @param coverages Population coverage for each day
#' @param n_days Length of simulation in days
#' @param distribution_type "random" or "cohort" distribution of IRS
get_irs_usage_mat <- function(days,coverages,n_days,distribution_type){
  n_dists <- length(coverages) #Number of irs distribution events
  ##Matrix of the population prop. protected by irs from distribution event i (col number) on day t (row number)
  irs_mat <- matrix(nrow = n_days, ncol = n_dists+1)
  irs_mat[1:(days[1]-1),] <- matrix(rep(c(rep(0,n_dists),1), (days[1]-1)), nrow = (days[1]-1), byrow = TRUE)  #Initialise no irs
  days_complete <- c(days, n_days)
  for(i in seq_along(coverages)){
    day <- days[i]
    new_recipients <- coverages[i]
    previous_row <- irs_mat[day-1,]
    new_row <- previous_row
    new_row[i] <- new_recipients
    if(distribution_type == "random"){
      new_row[-i] <- new_row[-i] - new_recipients*new_row[-i]
    } else if(distribution_type == "cohort"){
      current_cov <- 1-previous_row[n_dists+1]
      if(current_cov != 0){
        new_row[-i] <- new_row[-i] * (max(current_cov - new_recipients,0) / current_cov)
      }
      new_row[n_dists+1] <- 1-sum(new_row[1:n_dists])
    }
    days_in_distribution <- day:(days_complete[i+1]-1)
    irs_mat[days_in_distribution,] <- matrix(rep(new_row, length(days_in_distribution)), nrow = length(days_in_distribution), byrow = TRUE)
  }
  irs_mat[n_days,] <- irs_mat[(n_days-1),]
  return(irs_mat)
}

get_mean_irs_mats <- function(usage_mat, irs_mat){
  multiply <- usage_mat * irs_mat
  means <- rowSums(multiply)
  return(means)
}

spraying_decay <- function(t, theta, gamma) {
  1 / (1 + exp(-(theta + gamma * t)))
}

prob_spraying_repels <- function(ls_prime, ks_prime, js_prime, k0) {
  (1 - ks_prime / k0) * (js_prime / (ls_prime + js_prime))
}

prob_survives_spraying <- function(ks_prime, k0) {
  ks_prime / k0
}

get_daily_irs <- function(irs_usage_mat, parameter_mat){
  n_dists <- ncol(irs_usage_mat) - 1
  prop_cat <- irs_usage_mat[,1:n_dists] / rowSums(irs_usage_mat[,1:n_dists]) #Proportion of pop. in each distribution category
  mean_decay <- rowSums(prop_cat * parameter_mat[,1:n_dists])
  return(mean_decay)
}
