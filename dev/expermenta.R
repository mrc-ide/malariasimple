year <- 365
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 50
days <- c(183, 730)
retention <- 5*365
dn0 <- c(.2, .9)
rn <- c(.4, .8)
rnm <- c(.2, .5)
gamman <- c(1000, 2000)

dn0 <- c(.533, .533)
rn <- c(.56, .56)
rnm <- c(.24, .24)
gamman <- c(2.64,2.64)*365
coverages = c(.5, .5)
itn_params <- malariasimple::get_parameters(human_pop = human_population,
                                            n_days = sim_length,
                                            bites_Bed = 0.85,
                                            bites_Indoors = 0.97) |>
  malariasimple::set_bednets(days = days,
                         coverages = coverages,
                         retention = retention,
                         dn0 = dn0,
                         rn = rn,
                         rnm = rnm) |>
  malariasimple::set_equilibrium(init_EIR = starting_EIR)
itn_output <- malariasimple::run_simulation(itn_params) |> as.data.frame()
library(ggplot2)
ggplot() +
  #geom_line(aes(x = 1:1461, y = itn_params$r_itn_daily)) +
  geom_line(aes(x = 1:1461, y = itn_params$s_itn_daily))

test <- itn_params$r_itn_daily
plot(test)
plot(itn_params$r_itn_daily)

ggplot() +
  geom_line(data = itn_output, aes(x=time, y = n_detect_730_3650 / n_730_3650, col = "malariasimple")) +
  geom_vline(xintercept = days, lty = 2, col = "grey") +
  theme_bw()

##malariasimulation parameters
malsim_params <- malariasimulation::get_parameters(
  list(human_population = human_population)
) |> malariasimulation::set_equilibrium(init_EIR = starting_EIR)
malsim_itn_params <- malariasimulation::set_bednets(
    malsim_params,
    timesteps = days,
    coverages = coverages,
    retention = retention,
    dn0 = matrix(dn0, nrow = length(dn0), ncol = 1),
    rn = matrix(rn, nrow = length(dn0), ncol = 1),
    rnm = matrix(rnm, nrow = length(dn0), ncol = 1),
    gamman = gamman
  )
bednet_output <- malariasimulation::run_simulation_with_repetitions(timesteps = sim_length,
                                                                    repetitions = 5,
                                                                    overrides = malsim_itn_params,
                                                                    parallel = TRUE)
itn_output2 <- readRDS("C:/Users/Debbie/Downloads/itn_test.rds")
ggplot() +
  geom_line(data = bednet_output, aes(x=timestep, y = n_detect_lm_730_3650 / n_age_730_3650, col = "malariasimulation", group = repetition)) +
  geom_line(data = itn_output, aes(x=time, y = n_detect_730_3650 / n_730_3650, col = "malariasimple")) +
  geom_line(data = itn_output2, aes(x=time, y = n_detect_730_3650 / n_730_3650, col = "malariasimple - original")) +
  geom_vline(xintercept = days, lty = 2, col = "grey") +
  theme_bw()
ggplot() +
  #geom_line(data = bednet_output, aes(x=timestep, y = n_detect_lm_730_3650 / n_age_730_3650, col = "malariasimulation", group = repetition)) +
  geom_line(aes(x=1:1461, y = itn_params$s_itn_daily, col = "malariasimple")) +
  geom_line(data = itn_output2, aes(x=time, y = s_itn_out, col = "malariasimple - original"), lty = 2) +
  geom_vline(xintercept = days, lty = 2, col = "grey") +
  theme_bw()


