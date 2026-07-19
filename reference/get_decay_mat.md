# Decay matrix for ITN or SMC

Produces an ij matrix representing the decay in efficacy of
insecticide-treated nets (ITN) or prophylactic protection (SMC) on day j
for individuals in receipt of distribution event i.

## Usage

``` r
get_decay_mat(
  days,
  n_days,
  gamman_itn = NULL,
  scale_smc = NULL,
  shape_smc = NULL,
  intervention = NULL
)
```

## Arguments

- days:

  Days on which distribution events occur

- n_days:

  Length of simulation (days)

- gamman_itn:

  Half-life of ITN insecticide (days)

- scale_smc:

  Scale parameter of Weibull distribution defining decay of prophylactic
  protection of SMC

- shape_smc:

  Shape parameter of Weibull distribution defining decay of prophylactic
  protection of SMC

- intervention:

  "ITN" or "SMC"
