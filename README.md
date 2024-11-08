# malariasimple

A fast, time-discrete compartmental approximation of [malariasimulation](https://mrc-ide.github.io/malariasimulation/). 

# Installation
------

# Usage
Performing a model simulation in malariasimple is a three-step process:

```{r, eval=FALSE}
library(malariasimple)

#1. Generate model
gen <- gen_model()

#2. Define model parameters
params <- get_parameters() |>
          set_equilibrium(init_EIR = 10)

#3. Run model
out <- run_simulation(gen,params)
```

# Code Organisation
The model itself is defined in a single script (`inst/malariasimple_deterministic.R`) written with [odin.dust](https://mrc-ide.github.io/odin.dust/).

Customisations occur in the parameter set-up stage, which is facilitated by the helper functions found in `R/`.

# License
[MIT](https://choosealicense.com/licenses/mit/)






