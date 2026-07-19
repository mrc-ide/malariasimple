# Run deterministic malariasimple

Generates and runs the deterministic malariasimple model

## Usage

``` r
run_simulation(params, n_particles = 1, full_output = FALSE)
```

## Arguments

- params:

  List of model parameters generated using get_parameters() and other
  helper functions

- n_particles:

  Number of simulations to perform

- full_output:

  Boolean variable stating whether model should output model outputs, or
  a selected few.

## Examples

``` r
params <- get_parameters() |>
          set_equilibrium(init_EIR = 10)
simulation_output <- run_simulation(params)
```
