# malariasimple 0.1.0

## New features

* Added a stochastic version of the model, selectable via `get_parameters(stochastic = TRUE)`.
* Added the ability to propagate parameter uncertainty by drawing from the joint
  posterior distribution via the `parameter_draws` argument of `get_parameters()`
  (defaults to `"median"`, which reproduces the previous central estimates).
  See the new `parameter_draws_df` dataset.
* Added a fully correlated option for ITN distribution via
  `set_bednets(distribution_type = "correlated")` (default remains `"random"`).
* Added seasonal malaria chemoprevention (SMC) via `set_smc()`.
* Added support for a time-varying extrinsic incubation period (`eip`).
* Added custom output post-processing helpers (`get_custom_output()`, `make_2d()`).
* `run_simulation()` now errors clearly if `set_equilibrium()` has not been called.

## Bug fixes

* Corrected the ITN insecticide decay rate to use the mean retention time
  (rate `1 / gamman`), matching `malariasimulation`'s `bednet_decay()`.
* Fixed a double-counting error in ITN decay when only a single distribution
  event was specified.

## Breaking changes

These changes alter model behaviour and/or the public interface. Existing
scripts will need updating, and numerical results may differ from previous
versions.

* Several `get_parameters()` arguments were renamed:
  * `het_brackets` -> `biting_groups`
  * `ft` -> `daily_ft`
  * `delayMos` -> `eip`
  * `mu0` -> `mum`
  * `tau1` / `tau2` -> `foraging_time` / `gonotrophic_cycle`
  * `bites_Bed` / `bites_Indoors` -> `phi_bednets` / `phi_indoors`

  Old argument names are no longer recognised and are silently ignored.
* Some default parameter values changed: `tsd` (4 -> 3), number of biting
  heterogeneity groups (5 -> 4), `phi_bednets` (0.89 -> 0.85),
  `phi_indoors` (0.97 -> 0.9). A new `human_pop` argument (default 100000)
  now scales the initial human state sizes.
* `set_seasonality()` now takes Fourier coefficients `g0`, `g` and `h`
  directly instead of a `season_params` object. The `get_seasonal_forcing()`
  and `get_peak_cc()` signatures changed accordingly.
