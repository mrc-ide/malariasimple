#' Parameter draws (P. falciparum)
#'
#' Draws from the joint posterior distribution of the *Plasmodium falciparum*
#' transmission model fit, used to propagate parameter uncertainty via the
#' `parameter_draws` argument of [get_parameters()].
#'
#' @format ## `parameter_draws_df`
#' A data frame with 21,000 rows (1000 draws x 21 parameters) and 5 columns:
#' \describe{
#'   \item{draw}{Integer draw identifier (1 to 1000).}
#'   \item{malsim_name}{Parameter name as used in \pkg{malariasimulation}.}
#'   \item{malsim_val}{Drawn parameter value (malariasimulation parameterisation).}
#'   \item{simple_name}{Corresponding parameter name as used in malariasimple.}
#'   \item{simple_val}{Drawn parameter value (malariasimple parameterisation).}
#' }
#'
#' @source <https://www.nature.com/articles/ncomms4136>
"parameter_draws_df"
