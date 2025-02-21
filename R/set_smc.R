#' Set SMC parameters
#'
#' @param params malariasimple parameters
#' @param drug Specify default drug parameters. See ?set_drug_params for details
#' @param smc_cov Numeric vector detailing proportion of target age population receiving SMC for each distribution event. A single value is also accepted
#' @param smc_days Numeric vector detailing days on which SMC is administered
#' @param smc_lower_age Minimum age eligible for SMC
#' @param smc_upper_age Maximum age eligible for SMC
#' @param smc_clearance_lag No. of days between receipt of SMC and infectious clearance
#' @param drug_efficacy Proportion of infections successfully cleared
#' @param drug_rel_c Relative infectivity of treated individuals during infection clearing period
#' @param shape_smc Shape of Weibull distribution describing prophylactic waning
#' @param scale_smc Scale of Weibull distribution describing prophylactic waning
#' @param distribution_type SMC distributions can either be delived randomly across the eligible population each time "random", or within a consistent cohort "cohort".
#' @returns Updates the input parameter list to include SMC parameters
#'
#' @examples
#'
#' smc_params <- get_parameters(n_days = 300) |>
#'           set_smc(smc_days = c(100,150,200),
#'                   smc_cov = 0.4,
#'                   distribution_type = "cohort") |>
#'           set_equilibrium(init_EIR = 10)
#' @export
set_smc <- function(params,
                    drug = "SP_AQ",
                    smc_cov = NULL,
                    smc_days = NULL,
                    smc_lower_age = 0.25*365,
                    smc_upper_age = 5*365,
                    smc_clearance_lag = 5,
                    drug_efficacy = NULL,
                    drug_rel_c = NULL,
                    shape_smc = NULL,
                    scale_smc = NULL,
                    distribution_type = "random"){

  #----------------------------------- Sanity Checks -----------------------------------------------
  if(is.null(smc_days)){stop(message("smc_days must provided"))}
  if(is.null(smc_cov)){stop(message("smc_cov must be provided"))}

  if(!smc_lower_age %in% params$age_vector){stop(message("smc_lower_age must correspond to an age category boundary. Consider adjusting age_vector"))}
  if(!smc_upper_age %in% params$age_vector){stop(message("smc_upper_age must correspond to an age category boundary. Consider adjusting age_vector"))}

  if(max(smc_cov) > 1 | min(smc_cov) < 0) stop(message("smc_cov values must be between 0 and 1"))
  if(min(smc_days) < 1) stop("smc_days must be greater than 1")
  if(max(smc_days) > (params$n_days - smc_clearance_lag)){
    warning(message("SMC_days which exceed (n_days - smc_clearance_lag) have been removed."))
    smc_days <- smc_days[smc_days <= (params$n_days - smc_clearance_lag)]
  }
  if(length(smc_days) == 0) stop("No SMC distribution events are included")
  if(!(length(smc_cov) == 1 | length(smc_cov) == length(smc_days))) stop(message("length(smc_cov) must be either 1 or length(smc_days)"))
  smc_days <- round(smc_days) #It's easiest if SMC days are all integer values.
  smc_days <- smc_days + 1 #To account for day 0.
  #----------------------------------- Set Parameters -----------------------------------------------
  drug_params <- set_drug_params(drug)
  if(is.null(drug_efficacy)){drug_efficacy <- drug_params[1]}
  if(is.null(drug_rel_c)){drug_rel_c <- drug_params[2]}
  if(is.null(shape_smc)) {shape_smc <- drug_params[3]}
  if(is.null(scale_smc)) {scale_smc <- drug_params[4]}
  if(length(smc_cov) == 1) smc_cov <- rep(smc_cov,length(smc_days))

  if(drug_efficacy > 1 | drug_efficacy < 0) stop(message("drug_efficacy must be between 0 and 1"))
  if(drug_rel_c > 1 | drug_rel_c < 0) stop(message("drug_rel_c must be between 0 and 1"))
  params$smc_age <- params$age_vector[params$age_vector >= smc_lower_age & params$age_vector < smc_upper_age] #age_cats to receive SMC

  #----------------------------------- Prophylaxis Calculations ------------------------------------------------
  usage_mat <- get_smc_usage_mat(smc_days,smc_cov,params$n_days,distribution_type)
  daily_smc_cov <- get_daily_cov(usage_mat)

  decay_mat <- get_decay_mat(days = smc_days,n_days=params$n_days,scale_smc=scale_smc,shape_smc=shape_smc,intervention="SMC")
  daily_smc_decay <- get_daily_decay(usage_mat, decay_mat)
  params$P_smc_daily <- c(0,daily_smc_decay) #Prophylactic effect
  params$max_smc_cov <- max(daily_smc_cov) #Size of SMC compartment

  #------------------------------------- Infection Clearance ----------------------------------------------------
  #Schedule infection clearance
  # params$alpha_smc <- rep(0,(params$n_ts+1))
  # params$alpha_smc[(smc_days + smc_clearance_lag)*params$tsd] <- drug_efficacy
  # params$alpha_smc[params$alpha_smc != 0] <- params$alpha_smc[params$alpha_smc != 0] * (smc_cov / max(daily_smc_cov))

  eff_drug_efficacy <- (smc_cov / max(daily_smc_cov)) * drug_efficacy
  smc_days_mat <- cbind(
    time = c(0,smc_days + smc_clearance_lag, smc_days + smc_clearance_lag + params$dt),
    alpha_smc = c(0,eff_drug_efficacy,rep(0,length(smc_days))))
  sorted_smc_days_mat <- smc_days_mat[order(smc_days_mat[, 'time']), ]
  params$alpha_smc_times <- sorted_smc_days_mat[,"time"]
  params$alpha_smc_set <- sorted_smc_days_mat[,"alpha_smc"]

  #------------------------------------- Reduced Infectivity ---------------------------------------------------
  #The in the time between smc treatment and infection clearance - existing infections are less infectious by a factor of #drug_rel_c
  rel_c_days <- rep(1,params$n_days + 1)
  rel_c_days[as.vector(outer(smc_days, 0:(smc_clearance_lag-1), "+"))] <- drug_rel_c
  params$rel_c_days <- rel_c_days #Time vector. Days where SMC influences infectivity are set to SMC_rel_c. Else 1.

  #--------------------------------------------------------------------------------------------------------------
  params$smc_set <- 1
  return(params)
}


get_smc_usage_mat <- function(smc_days,smc_cov,n_days,distribution_type){
  n_dists <- length(smc_cov) #Number of SMC distribution events

  ##Matrix of the population prop. protected by smc from distribution event i (col number) on day t (row number)
  smc_mat <- matrix(nrow = n_days, ncol = n_dists+1)
  smc_mat[1,] <- c(rep(0,n_dists),1) #Initialise no smc
  days_since_last_adm <- 0
  for(i in 2:n_days){
    smc_mat[i,1:n_dists] <- smc_mat[(i-1),(1:n_dists)]
    smc_mat[i,(n_dists+1)] <- 1 - sum(smc_mat[i,(1:n_dists)])
    days_since_last_adm <- days_since_last_adm + 1

    if(distribution_type == "random"){
      ##After 100 days without an SMC administration, the system resets
      if(days_since_last_adm == 150){
        smc_mat[i,] <- c(rep(0,n_dists),1)
      }}

    #New SMC are distributed on smc_days
    if(i %in% smc_days){
      dist_index <- which(smc_days == i) #Which distribution event
      new_recipients <- smc_cov[dist_index] #Proportion which is now a 'new' recipient
      smc_mat[i,dist_index] <- new_recipients
      if(distribution_type == "random"){
        smc_mat[i,-dist_index] <- smc_mat[i,-dist_index] - new_recipients*smc_mat[i,-dist_index]
      } else if(distribution_type == "cohort"){
        smc_mat[i,-dist_index] <- 0
        smc_mat[i,ncol(smc_mat)] <- 1-new_recipients
      } else{
        stop(message("SMC distribution type must be either 'random' or 'cohort'"))
      }
      days_since_last_adm <- 0
    }
  }
  return(smc_mat)
}

##See https://github.com/mrc-ide/malariasimulation/blob/master/R/drug_parameters.R for references
#' @title Set drug parameter defaults
#' @description Sets default parameters for commonly used SMC drugs. Values from SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014).
#
#' @param drug A character string specifying the drug type. One of:
#' \itemize{
#'   \item "SP_AQ": Sulfadoxine-pyrimethamine with Amodiaquine [0.9, 0.32, 4.3, 38.1]
#'   \item "AL": Artemether-Lumefantrine [0.95,0. 05094, 11.3, 10.6]
#'   \item "DHA_PQP": Dihydroartemisinin-Piperaquine [0.95, 0.09434, 4.4, 28.1]
#' }
#'
#' @return A vector of drug parameters for the specified drug in the format
#' [drug_efficacy, drug_rel_c, prophylaxis_shape, prophylaxis_scale]
#'
#' @examples
#' set_drug_params(drug = "SP_AQ")

#' @export
set_drug_params <- function(drug){
  #SMC drug parameters [drug_efficacy, drug_rel_c, drug_prophylaxis_shape, drug_prophylaxis_scale]
  if(drug == "SP_AQ"){
    drug_params <- c(0.9, 0.32, 4.3, 38.1) #sulphadoxine-pyrimethamine and amodiaquine
  } else if(drug == "AL"){
    drug_params <- c(.95, 0.05094, 11.3, 10.6) #artemether-lumefantrine
  } else if(drug == "DHA_PQP"){
    drug_params <- c(.95, 0.09434, 4.4, 28.1) #dihydroartemisinin-piperaquine
  } else{
    stop(message("SMC drug '",drug,"' not recognised. Available pre-sets are `SP_AQ`, `AL`, and `DHA_PQP`.\n",
                 "Custom drug allocation may be applied by defining parameters directly"))
  }
  names(drug_params) <- c("drug_efficacy", "drug_rel_c", "drug_prophylaxis_shape", "drug_prophylaxis_scale")
  return(drug_params)
}
