# Set equilibrium

Updates parameter list to include equilibrium values for a given EIR
(entomological inoculation rate) for the malariasimple parameter list.
This function also includes some 'finishing touches' for the smooth
running of the simulation and so is essential that this function is
performed last.

## Usage

``` r
set_equilibrium(params, init_EIR)
```

## Arguments

- params:

  List of parameters

- init_EIR:

  Value of EIR at equilbrium

## Examples

``` r
params <- get_parameters() |>
            set_equilibrium(init_EIR = 6)
```
