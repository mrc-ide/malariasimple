# Set SMC parameters

Set SMC parameters

## Usage

``` r
set_smc(
  params,
  drug = "SP_AQ",
  coverages = NULL,
  days = NULL,
  min_age = 0.25 * 365,
  max_age = 5 * 365,
  smc_clearance_lag = 5,
  drug_efficacy = NULL,
  drug_rel_c = NULL,
  shape_smc = NULL,
  scale_smc = NULL,
  distribution_type = "random"
)
```

## Arguments

- params:

  malariasimple parameters

- drug:

  Specify default drug parameters. See ?set_drug_params for details

- coverages:

  Numeric vector detailing proportion of target age population receiving
  SMC for each distribution event. A single value is also accepted

- days:

  Numeric vector detailing days on which SMC is administered

- min_age:

  Minimum age eligible for SMC

- max_age:

  Maximum age eligible for SMC

- smc_clearance_lag:

  No. of days between receipt of SMC and infectious clearance

- drug_efficacy:

  Proportion of infections successfully cleared

- drug_rel_c:

  Relative infectivity of treated individuals during infection clearing
  period. 1 indicates no difference. 0 indicates no forward infection.

- shape_smc:

  Shape of Weibull distribution describing prophylactic waning

- scale_smc:

  Scale of Weibull distribution describing prophylactic waning

- distribution_type:

  SMC distributions can either be delivered randomly across the eligible
  population each time "random", or given to the same children each time
  "correlated".

## Value

Updates the input parameter list to include SMC parameters

## Examples

``` r

smc_params <- get_parameters(n_days = 300) |>
          set_smc(days = c(100,150,200),
                  coverages = 0.4,
                  distribution_type = "correlated") |>
          set_equilibrium(init_EIR = 10)
```
