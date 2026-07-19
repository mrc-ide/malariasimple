# Get more specific outputs

Allows the user to filter to more specific human outputs. Aggregated
values for any human compartment as well as detected cases and clinical
incidence may be requested for specific age, biting heterogeneity or
intervention groups.

## Usage

``` r
get_custom_output(
  params,
  sim,
  output_variables,
  ages = NULL,
  biting_groups = NULL,
  int_groups = NULL
)
```

## Arguments

- params:

  malariasimple input parameters

- sim:

  Output from malariasimple::run_simulation(params, full_output = TRUE)

- output_variables:

  Output variable of interest.

- ages:

  Age groups of interest (defined by the lower bracket)

- biting_groups:

  Biting groups of interest

- int_groups:

  Intervention groups of interest, given as indices into the model's
  intervention strata (1 = No intervention, 2 = Bednets, 3 = SMC, 4 =
  Bednets and SMC). These indices reflect the interventions currently
  supported by the model (ITNs and SMC); adding a new intervention type
  would require extending the intervention structure in
  \`set_equilibrium()\` and the odin model, after which the index
  meanings would change accordingly.

## Examples

``` r
  params <- get_parameters(
    n_days = 50,
    age_vector = c(0, 0.25, 0.5, 1, 3, 5, 10, 50) * 365) |>
    set_equilibrium(init_EIR = 10)
  sim <- run_simulation(params, full_output = TRUE)
  custom_output <- get_custom_output(
    params,
    sim,
    output_variables = c("n", "detect"),
    ages = c(0, 0.5) * 365)
```
