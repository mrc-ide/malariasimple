# Set ITN parameters

Set ITN parameters

## Usage

``` r
set_bednets(
  params,
  continuous_distribution = FALSE,
  daily_continuous_cov = NULL,
  days = NULL,
  coverages = NULL,
  gamman = 2.64 * 365,
  retention = 3 * 365,
  dn0 = 0.41,
  rn = 0.56,
  rnm = 0.24,
  distribution_type = "random"
)
```

## Arguments

- params:

  malariasimple parameters

- continuous_distribution:

  Is ITN distribution continuous? If FALSE, distribution is assumed to
  occur in discrete events.

- daily_continuous_cov:

  Vector of daily ITN coverage (required when continuous_distribution =
  TRUE). A single value is also accepted

- days:

  Vector of days on which ITN distribution events occur (required when
  continuous_distribution = FALSE). Analogous to 'timesteps' argument in
  malariasimulation

- coverages:

  Vector detailing the proportion of the population receiving an ITN
  during each intervention (required when continuous_distribution =
  FALSE)

- gamman:

  Mean lifetime ITN insecticide efficacy (days).

- retention:

  Average number of days a net is kept for

- dn0:

  Probability of mosquito dying upon an encounter with ITN (max)

- rn:

  Probability of repeating behaviour with ITN (max)

- rnm:

  Probability of repeating behaviour with ITN (min)

- distribution_type:

  Either 'random' or 'correlated'

## Value

Updates the input parameter list to include ITN parameters

## Examples

``` r
n_days <- 500
#Discrete distribution scenario
discrete_itn_params <- get_parameters(n_days = n_days) |>
  set_bednets(days = c(50,100,200),
          coverages = c(0.2,0.5,0.1)) |>
  set_equilibrium(init_EIR = 10)

#Continuous distribution scenario
continuous_cov <- 0.2 + 0.15*(sin(2 * pi * (1:n_days / 365)) + 1)
discrete_itn_params <- get_parameters(n_days = n_days) |>
  set_bednets(continuous_distribution = TRUE,
          daily_continuous_cov = continuous_cov) |>
  set_equilibrium(init_EIR = 10)
```
