# Parameter draws (P. falciparum)

Draws from the joint posterior distribution of the \*Plasmodium
falciparum\* transmission model fit, used to propagate parameter
uncertainty via the \`parameter_draws\` argument of
\[get_parameters()\].

## Usage

``` r
parameter_draws_df
```

## Format

\## \`parameter_draws_df\` A data frame with 21,000 rows (1000 draws x
21 parameters) and 5 columns:

- draw:

  Integer draw identifier (1 to 1000).

- malsim_name:

  Parameter name as used in malariasimulation.

- malsim_val:

  Drawn parameter value (malariasimulation parameterisation).

- simple_name:

  Corresponding parameter name as used in malariasimple.

- simple_val:

  Drawn parameter value (malariasimple parameterisation).

## Source

\<https://www.nature.com/articles/ncomms4136\>
