# Converts 3D to long 2D

Takes a 3D array simulation output with multiple particles, and melts
the third dimension to produce a 2D dataframe.

## Usage

``` r
make_2d(array_output)
```

## Arguments

- array_output:

  Any 3D array

## Examples

``` r
params <- get_parameters(stochastic = TRUE) |>
            set_equilibrium(init_EIR = 5)
long_output <- run_simulation(params, n_particles = 3) |> make_2d()
```
