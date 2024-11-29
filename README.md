
<!-- README.md is generated from README.Rmd. Please edit that file -->

# malariasimple

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/malariasimple/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/malariasimple/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/malariasimple/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/malariasimple)
<!-- badges: end -->

The goal of malariasimple is to …

## Installation

You can install the development version of malariasimple from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("mrc-ide/malariasimple")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
 #library(malariasimple)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:
=======
A fast, time-discrete compartmental approximation of [malariasimulation](https://mrc-ide.github.io/malariasimulation/). 

# Installation
------

# Usage
Performing a model simulation in malariasimple is a three-step process:

```{r, eval=FALSE}
library(malariasimple)

#1. Generate model
gen <- gen_model()

#2. Define model parameters
params <- get_parameters() |>
          set_equilibrium(init_EIR = 10)

#3. Run model
out <- run_simulation(gen,params)
```

# Code Organisation
The model itself is defined in a single script (`inst/malariasimple_deterministic.R`) written with [odin.dust](https://mrc-ide.github.io/odin.dust/).

Customisations occur in the parameter set-up stage, which is facilitated by the helper functions found in `R/`.

# License
[MIT](https://choosealicense.com/licenses/mit/)

