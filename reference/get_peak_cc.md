# Timing of peak carrying capacity

Returns the day-of-year on which vector carrying capacity is maximum

## Usage

``` r
get_peak_cc(g0, g, h)
```

## Arguments

- g0:

  Mean baseline coefficient

- g:

  Cosine coefficients

- h:

  Sine coefficients

## Examples

``` r
#Seasonality parameters
g0 <- 0.28
g <- c(-0.3, -0.03, 0.17)
h <- c(-0.35, 0.32, -0.07)

peak_cc <- get_peak_cc(g0, g, h)
```
