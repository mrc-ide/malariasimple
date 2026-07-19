# Add seasonal patterns

Adds seasonal forcing to the model by approximating average seasonal
rainfall as a Fourier series.

## Usage

``` r
set_seasonality(params, g0, g, h, floor = 0.001)
```

## Arguments

- params:

  Other malariasimple parameters

- g0:

  Mean baseline coefficient

- g:

  Cosine coefficients

- h:

  Sine coefficients

- floor:

  Minimum value of seasonal forcing

## Value

Updates the input parameter list to include seasonal parameters

## Examples

``` r
# Define seasonality Fourier coefficients
g0 <- 0.28
g <- c(-0.3, -0.03, 0.17)
h <- c(-0.35, 0.33, -0.08)
params <- get_parameters() |>
          set_seasonality(g0=g0,
          g=g,
          h=h,
          floor = 0.005) |>
          set_equilibrium(init_EIR = 5)
```
