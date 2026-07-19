# Create parameters for malariasimple

Helper function to provide defaults for most necessary values required
for running the malariasimple model

## Usage

``` r
get_parameters(
  stochastic = FALSE,
  parameter_draws = "median",
  n_days = 100,
  human_pop = 1e+05,
  tsd = 3,
  age_vector = c(0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8.5, 10, 20, 30, 40,
    60, 80) * 365,
  biting_groups = 4,
  lag_rates = 10,
  eta = 1/(21 * 365),
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100 * 365,
  rA = 1/195,
  rT = 0.2,
  rD = 0.2,
  rU = 1/110.299,
  rP = 1/15,
  dE = 12,
  delayGam = 12.5,
  cD = 0.0676909,
  cT = 0.322 * cD,
  cU = 0.006203,
  gamma1 = 1.82425,
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  eip = 10,
  foraging_time = 0.69,
  gonotrophic_cycle = 2.31,
  mum = 0.132,
  Q0 = 0.92,
  chi = 0.86,
  phi_bednets = 0.85,
  phi_indoors = 0.9,
  muEL = 0.0338,
  muLL = 0.0348,
  muPL = 0.249,
  dEL = 6.64,
  dLL = 3.72,
  dPL = 0.643,
  gammaL = 13.25,
  betaL = 21.2,
  daily_ft = 0,
  clin_inc_rendering_min_ages = NULL,
  clin_inc_rendering_max_ages = NULL,
  prevalence_rendering_min_ages = NULL,
  prevalence_rendering_max_ages = NULL,
  ...
)
```

## Arguments

- stochastic:

  Boolean variable. Set to false for deterministic simulation

- parameter_draws:

  Default 'median'. If given a value 1-1000 model parameters are taken
  with a single draw from the fitted joint posterior

- n_days:

  Number of days for which the simulation will run.

- human_pop:

  Size of human population (count)

- tsd:

  Number of time-steps per day. Fewer is faster, more better
  approximates the continuous solution.

- age_vector:

  Lower bound of each age category. See function for default.

- biting_groups:

  Number of biting heterogeneity groups.

- lag_rates:

  Number of sub-compartments within FOI and FOIv which approximate
  delay-differential equation. Higher values are a closer approximation,
  but computationally more expensive.

- eta:

  Death rate for exponential population distribution, i.e. 1/Mean
  Population Age

- rho:

  Age-dependent biting parameter

- a0:

  Age-dependent biting parameter

- sigma2:

  Variance of the log heterogeneity in biting rates

- max_age:

  Upper boundary of oldest age category

- rA:

  Rate of leaving asymptomatic infection

- rT:

  Rate of leaving treatment

- rD:

  Rate of leaving clinical disease

- rU:

  Rate of recovering from subpatent infection

- rP:

  Rate of leaving prophylaxis

- dE:

  Latent period of human infection

- delayGam:

  Lag from parasites to infectious gametocytes

- cD:

  Untreated disease contribution to infectiousness

- cT:

  Treated disease contribution to infectiousness

- cU:

  Subpatent disease contribution to infectiousness

- gamma1:

  Parameter for infectiousness of state A

- d1:

  Minimum probability due to maximum immunity

- dID:

  Inverse of decay rate

- ID0:

  Scale parameter

- kD:

  Shape parameter

- uD:

  Duration in which immunity is not boosted

- aD:

  Scale parameter relating age to immunity

- fD0:

  Time-scale at which immunity changes with age

- gammaD:

  Shape parameter relating age to immunity

- alphaA:

  PCR detection probability parameters state A

- alphaU:

  PCR detection probability parameters state U

- b0:

  Maximum probability due to no immunity

- b1:

  Maximum relative reduction due to immunity

- dB:

  Inverse of decay rate

- IB0:

  Scale parameter

- kB:

  Shape parameter

- uB:

  Duration in which pre-erythrocytic immunity is not boosted

- phi0:

  Maximum probability due to no immunity

- phi1:

  Maximum relative reduction due to immunity

- dCA:

  Inverse of decay rate

- IC0:

  Scale parameter

- kC:

  Shape parameter

- uCA:

  Duration in which clinical immunity is not boosted

- PM:

  New-born immunity relative to mothers

- dCM:

  Inverse of decay rate of maternal immunity

- eip:

  Extrinsic incubation period (days). Either a scalar value (constant
  EIP), or a vector of daily values.

- foraging_time:

  Duration of host seeking, assumed to be constant between species

- gonotrophic_cycle:

  Duration of mosquito resting after feed

- mum:

  Daily mortality of adult mosquitoes

- Q0:

  Anthrophagy probability

- chi:

  Endophily probability

- phi_bednets:

  Percentage of bites indoors and in bed

- phi_indoors:

  Percentage of bites indoors

- muEL:

  Per capita daily mortality rate of early stage larvae (low density)

- muLL:

  Per capita daily mortality rate of late stage larvae (low density)

- muPL:

  Per capita daily mortality rate of pupae

- dEL:

  Development time of early stage larvae

- dLL:

  Development time of late stage larvae

- dPL:

  Development time of pupae

- gammaL:

  Relative effect of density dependence on late instars relative to
  early instars

- betaL:

  Number of eggs laid per day per mosquito

- daily_ft:

  Daily vector of percentage of population that gets treated. A scalar
  value is permitted and assumes constant ft.

- clin_inc_rendering_min_ages:

  Vector of values (or singe value) of lower age boundaries for clinical
  incidence output (days)

- clin_inc_rendering_max_ages:

  Vector of values (or singe value) of upper age boundaries for clinical
  incidence output (days)

- prevalence_rendering_min_ages:

  Vector of values (or singe value) of lower age boundaries for
  prevalence output (days)

- prevalence_rendering_max_ages:

  Vector of values (or singe value) of upper age boundaries for
  prevalence output (days)

- ...:

  Additional arguments

## Examples

``` r
#Output default model for clinical incidence in range [0-10] and 2+
params <- get_parameters(n_days = 500,
                         clin_inc_rendering_min_ages = c(0,2*365),
                         clin_inc_rendering_max_ages = c(10*365,Inf)) |>
          set_equilibrium(init_EIR = 10)
```
