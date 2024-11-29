#' @title Run deterministic malariasimple
#' @description Generates and runs the deterministic malariasimple model
#' @param params List of model parameters generated using get_parameters() and other helper functions
#' @examples
#' params <- get_parameters() |>
#'           set_equilibrium(init_EIR = 10)
#' gen <- gen_model()
#' simulation_output <- run_simulation(gen,params)


run_simulation <- function(gen, params, full_output = FALSE) {
  ## Set up model
  mod <- gen$new(
    pars = params,
    time = 1,
    n_particles = 1,
    n_threads = 1,
    seed = 1L
  )

  mod_info <- mod$info()
  out <- array(NA, dim = c(mod_info$len, 1, params$n_ts)) #Initialize output array. Second dimension is set to 1 because model is deterministic

  ## Run model
  for (t in seq_len(params$n_ts)) {
    out[, , t] <- mod$run(t)
  }

  dim(out) <- c(dim(out)[1], dim(out)[3]) #Make it 2D
  out <- t(out)
  out <- out[!duplicated(out[,1]),] #Return only one output per day
  out <- out[-1,] #Discard the first row

  ## Get column names
  colname_length <- sum(sapply(mod_info$index, length))
  index <- vector("character", colname_length)
  pos <- 1

  for (i in seq_along(mod_info$index)) {
    len <- length(mod_info$index[[i]])
    index[pos:(pos + len - 1)] <- names(mod_info$index)[i]
    pos <- pos + len
  }

  #Improve naming of user-defined outputs
  index[index == "n_ud_prev"] <- paste("n",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_age,sep="_")
  index[index == "n_ud_detect_prev"] <- paste("n_detect",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_ages,sep="_")
  index[index == "n_ud_inc"] <- paste("n_clin_inc",params$age_vector[params$min_age_inc], params$clin_inc_rendering_max_ages,sep="_")
  colnames(out) <- index

  if(full_output == FALSE){
    selected_cols <- c("day","EIR_mean","natural_deaths","mu_mosq",
                       grep("_count$", colnames(out), value = TRUE),
                       "ica_mean","icm_mean","ib_mean", "id_mean",
                       grep("^n_", colnames(out), value = TRUE),
                       "EL","LL","PL","Sv","Pv","Iv","mv")

    out[,"natural_deaths"] <- out[,"natural_deaths"] * params$tsd
    out <- out[,selected_cols]
  }
  return(out)
}
