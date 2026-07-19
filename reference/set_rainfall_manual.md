# Add manual rainfall forcing

Takes daily rain-forcing time series input and adds it to the parameter
set

## Usage

``` r
set_rainfall_manual(params, cc_ts)
```

## Arguments

- params:

  Other malariasimple parameters

- cc_ts:

  Carrying capacity time series. Vector of daily carrying capacity.
  Length must equal or exceed params\$n_days

## Value

Updates the input parameter list to include seasonal parameters

## Examples

``` r
n_days = 1000
t <- 1:n_days
rainfall_ts <- sin((t*2*pi)/365) + 1.1
params <- get_parameters() |>
          set_rainfall_manual(rainfall_ts) |>
          set_equilibrium(init_EIR = 5)
```
