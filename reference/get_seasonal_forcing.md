# Get daily seasonal forcing

Convert Fourier coefficients into a smooth vector of daily seasonal
forcing. Used within set_seasonality function.

## Usage

``` r
get_seasonal_forcing(t, g0, g, h, floor = 0.001)
```

## Arguments

- t:

  Day-of-year

- g0:

  Mean baseline coefficient

- g:

  Cosine coefficients

- h:

  Sine coefficients

- floor:

  Minimum permitted value of output
