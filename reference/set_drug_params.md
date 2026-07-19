# Set drug parameter defaults

Sets default parameters for commonly used SMC drugs. Values from SI of
Commun. 5:5606 doi: 10.1038/ncomms6606 (2014).

## Usage

``` r
set_drug_params(drug)
```

## Arguments

- drug:

  A character string specifying the drug type. One of:

  - "SP_AQ": Sulfadoxine-pyrimethamine with Amodiaquine \[0.9, 0.32,
    4.3, 38.1\]

  - "AL": Artemether-Lumefantrine \[0.95,0. 05094, 11.3, 10.6\]

  - "DHA_PQP": Dihydroartemisinin-Piperaquine \[0.95, 0.09434, 4.4,
    28.1\]

## Value

A vector of drug parameters for the specified drug in the format
\[drug_efficacy, drug_rel_c, prophylaxis_shape, prophylaxis_scale\]

## Examples

``` r
set_drug_params(drug = "SP_AQ")
#>          drug_efficacy             drug_rel_c drug_prophylaxis_shape 
#>                   0.90                   0.32                   4.30 
#> drug_prophylaxis_scale 
#>                  38.10 
```
