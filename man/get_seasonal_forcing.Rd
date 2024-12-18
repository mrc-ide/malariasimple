% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasonality.R
\name{get_seasonal_forcing}
\alias{get_seasonal_forcing}
\title{Get seasonal forcing}
\usage{
get_seasonal_forcing(season_params, n_days, floor = 0.001)
}
\arguments{
\item{season_params}{List of Fourier function parameters describing seasonality}

\item{n_days}{Numerical value of output vector length}

\item{floor}{Minimum permitted value of output}
}
\description{
Convert Fourier function parameters into a smooth vector of daily seasonal forcing. Used within set_seasonality function
}
\examples{
season_params <- list(g0 = 0.28,
                          g = c(-0.30, -0.03, 0.17),
                          h = c(-0.35, 0.32, -0.08),
                          theta_c = 0.28)
get_seasonal_forcing(season_params,
                      n_days = 1000,
                      floor = 0.005)
}
