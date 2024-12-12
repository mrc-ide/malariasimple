#' @title Run deterministic malariasimple
#' @description Generates and runs the deterministic malariasimple model
#' @param params List of model parameters generated using get_parameters() and other helper functions
#' @param full_output Boolean variable stating whether model should output model outputs, or a selected few.
#' @examples
#' params <- get_parameters() |>
#'           set_equilibrium(init_EIR = 10)
#' simulation_output <- run_simulation(params)
#' @export

run_simulation <- function(params, full_output = FALSE){
  # model_path = "inst/odin/malariasimple_deterministic.R"
  # gen <- odin2::odin(model_path)
  gen <- malariasimple_deterministic

  sys <- dust2::dust_system_create(gen(), params, n_particles = 1, dt = 1/params$tsd)

  dust2::dust_system_set_state_initial(sys)
  time <- 0:params$n_days
  out <- dust2::dust_system_simulate(sys, time)
  out <- out[,-1] |> t()
  out <- cbind(time = 1:params$n_days, out)

  ## Get column names
  colname_ls <- dust2::dust_unpack_index(sys)
  colname_length <- sum(sapply(colname_ls, length))
  index <- vector("character", colname_length)
  pos <- 1

  for (i in seq_along(colname_ls)) {
    len <- length(colname_ls[[i]])
    index[pos:(pos + len - 1)] <- names(colname_ls)[i]
    pos <- pos + len
  }

  #Improve naming of user-defined outputs
  index[index == "n_ud_prev"] <- paste("n",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_age,sep="_")
  index[index == "n_ud_detect_prev"] <- paste("n_detect",params$age_vector[params$min_age_prev], params$prevalence_rendering_max_ages,sep="_")
  index[index == "n_ud_inc"] <- paste("n_clin_inc",params$age_vector[params$min_age_inc], params$clin_inc_rendering_max_ages,sep="_")
  colnames(out) <- c("time",index)

  if(full_output == FALSE){
    selected_cols <- c("time","EIR_mean","natural_deaths","mu_mosq",
                       grep("_count$", colnames(out), value = TRUE),
                       "ica_mean","icm_mean","ib_mean", "id_mean",
                       grep("^n_", colnames(out), value = TRUE),
                       "EL","LL","PL","Sv","Pv","Iv","mv")

    out[,"natural_deaths"] <- out[,"natural_deaths"] * params$tsd
    out <- out[,selected_cols]
  }
  return(out)
}
