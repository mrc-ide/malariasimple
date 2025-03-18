
# malariasimple

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/malariasimple/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/malariasimple/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/malariasimple/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/malariasimple)
<!-- badges: end -->

A fast, time-discrete compartmental approximation of
[malariasimulation](https://mrc-ide.github.io/malariasimulation/).

## Installation

You can install the development version of malariasimple from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("mrc-ide/malariasimple")
```

## Usage

Performing a model simulation in malariasimple is a two-step process:

``` r
library(malariasimple)
#1. Define model parameters
params <- get_parameters() |>
          set_equilibrium(init_EIR = 10)

#2. Run model
out <- run_simulation(params)
```

## Code Organisation

The model itself is written using
[odin2](https://mrc-ide.github.io/odin.dust/) and is found in two
scripts:

- `inst/odin/malariasimple_deterministic.R`
- `inst/odin/malariasimple_stochastic.R`

Customisations occur in the parameter set-up stage, which is facilitated
by the helper functions found in `R/`.

## License

[MIT](https://choosealicense.com/licenses/mit/)
