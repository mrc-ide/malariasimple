# Produce matrix of population ITN usage

Produces an ij matrix defining the proportion of the population using an
ITN from distribution i on day j. Final column represents population
with no ITN.

## Usage

``` r
get_itn_usage_mat(days, coverages, retention, n_days, distribution_type)
```

## Arguments

- days:

  Vector of days on which ITN distribution events occur (required when
  continuous_distribution = FALSE)

- coverages:

  Vector detailing the proportion of the population receiving an ITN
  during each intervention (required when continuous_distribution =
  FALSE)

- retention:

  Average number of days a net is kept for

- n_days:

  Length of simulation

- distribution_type:

  Either 'random' or 'correlated'
