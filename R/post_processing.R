#' @title Get more specific outputs
#' @description Allows the user to filter to more specific human outputs. Aggregated values for any human compartment as well as
#' detected cases and clinical incidence may be requested for specific age, biting heterogeneity or intervention groups.
#' @param params malariasimple input parameters
#' @param sim Output from malariasimple::run_simulation(params, full_output = TRUE)
#' @param output_variable Output variable of interest.
#' @param ages Age groups of interest (defined by the lower bracket)
#' @param biting_groups Biting groups of interest
#' @param int_groups Intervention groups of interest. 1 = No intervention, 2 = Bednets, 3 = SMC, 4 = Bednets and SMC
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt

get_custom <- function(params,
                       sim,
                       output_variable,
                       ages = NULL,
                       biting_groups = NULL,
                       int_groups = NULL) {
  grid <- expand.grid(
    min_age = params$age_vector,
    biting_group = 1:params$biting_groups,
    int_group = 1:params$num_int
  )
  grid$idx <- 1:nrow(grid)
  if (!is.null(ages))
    grid <- grid[grid$min_age %in% ages, ]
  if (!is.null(biting_groups))
    grid <- grid[grid$biting_group %in% biting_groups, ]
  if (!is.null(int_groups))
    grid <- grid[grid$int_group %in% int_groups, ]

  colnames <- colnames(sim)
  if(length(dim(sim)) == 3){
    n_particles <- dim(sim)[3]
    output <- sim[, colnames %in% c(output_variable),]
    output_filtered <- output[, grid$idx, ]
    output_summed <- apply(output_filtered, c(1,3), sum) |>
      as.data.frame() |>
      mutate(time = sim[,"time",1])
    colnames(output_summed) <- c(1:n_particles, "time")
    out_df <- output_summed |> melt(id.vars = "time")
    colnames(out_df) <- c("time", "particle", output_variable)
  } else {
    output <- sim[, colnames == output_variable]
    output_filtered <- output[, grid$idx, ]
    out_df <- data.frame(sim[, "time"], rowSums(output_filtered))
    colnames(out_df) <- c("time", output_variable)
  }

  return(out_df)
}

#' @title Get more specific outputs
#' @description Allows the user to filter to more specific human outputs. Aggregated values for any human compartment as well as
#' detected cases and clinical incidence may be requested for specific age, biting heterogeneity or intervention groups.
#' @param params malariasimple input parameters
#' @param sim Output from malariasimple::run_simulation(params, full_output = TRUE)
#' @param output_variables Output variable of interest.
#' @param ages Age groups of interest (defined by the lower bracket)
#' @param biting_groups Biting groups of interest
#' @param int_groups Intervention groups of interest. 1 = No intervention, 2 = Bednets, 3 = SMC, 4 = Bednets and SMC
#' @export

get_custom_output <- function(params,
                              sim,
                              output_variables,
                              ages = NULL,
                              biting_groups = NULL,
                              int_groups = NULL) {
  var_list <- vector(mode = "list", length = length(output_variables))
  for(i in 1:length(output_variables)){
    var_list[[i]] <- get_custom(
      params = params,
      sim = sim,
      output_variable = output_variables[i],
      ages = ages,
      biting_groups = biting_groups,
      int_groups = int_groups
    )
  }
  df <- Reduce(function(x,y) merge(x,y), var_list)
  return(df)
}


#' @title Converts 3D to long 2D
#' @description Takes a 3D array simulation output with multiple particles, and melts the third dimension to produce a 2D dataframe.
#' @param array_output Any 3D array
#' @examples
#' params <- get_parameters(stochastic = TRUE) |>
#'             set_equilibrium(init_EIR = 5)
#' long_output <- run_simulation(params, n_particles = 3) |> make_2d()
#' @export
make_2d <- function(array_output){
  n_particles <- dim(array_output)[3]
  n_days <- dim(array_output)[1]
  out_2d <- do.call(rbind, lapply(1:n_particles, function(i) array_output[,,i])) |>
    as.data.frame()
  out_2d$particle <- rep(1:n_particles, each = n_days)
  return(out_2d)
}
