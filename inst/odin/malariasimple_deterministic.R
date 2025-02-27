## DECLARE TIME STEPS
n_days <- parameter()
n_ts <- parameter()

dim(daily_rain_input) <- n_days +1
daily_rain_input <- parameter()

dim(days) <- n_days + 1
days <- parameter()

dim(iterations) <- n_ts + 1
iterations <- parameter()


## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- parameter()
nh <- parameter()
ft <- parameter()

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------
# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

#--------- DEFINE PARAMETERS FOR HUMAN TRANSITIONS ---------
eta <- parameter()
age_rate <- parameter()
dim(age_rate) <- na
het_wt <- parameter()
dim(het_wt) <- nh
rA <- parameter()
rT <- parameter()
rD <- parameter()
rU <- parameter()
rP <- parameter()

#----------- INITIALISE HUMAN COMPARTMENTS -----------------
# S - SUSCEPTIBLE
init_S <- parameter()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

# T- SUCCESSFULLY TREATED
init_T <- parameter()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

# D - CLEAR DISEASE
init_D <- parameter()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

# A - ASYMPTOMATIC DISEASE
init_A <- parameter()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

# U - SUBPATENT DISEASE
init_U <- parameter()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

# P - PROPHYLAXIS
init_P <- parameter()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

#-------------------- DEFINE HUMAN TRANSITION EVENTS -------------------------
#### BIRTHS ####
dim(births) <- c(1,nh,num_int)
births[1,,] <- H*dt*cov[k]*eta*het_wt[j]


#### TRANSITIONS FROM S #####
# Rates of transition
dim(ST_rate) <- c(na,nh,num_int)
dim(SD_rate) <- c(na,nh,num_int)
dim(SA_rate) <- c(na,nh,num_int)

ST_rate[,,] <- ft*phi[i,j,k]*FOI_smc[i,j,k]
SD_rate[,,] <- (1-ft)*phi[i,j,k]*FOI_smc[i,j,k]
SA_rate[,,] <- (1-phi[i,j,k])*FOI_smc[i,j,k]


# Number of humans transitioning per timestep
dim(ST_trans) <- c(na,nh,num_int)
dim(SD_trans) <- c(na,nh,num_int)
dim(SA_trans) <- c(na,nh,num_int)
dim(S_death) <- c(na,nh,num_int)
dim(S_age) <- c(na,nh,num_int)

ST_trans[,,] <- if(S[i,j,k]*dt*ST_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*ST_rate[i,j,k]
SD_trans[,,] <- if(S[i,j,k]*dt*SD_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SD_rate[i,j,k]
SA_trans[,,] <- if(S[i,j,k]*dt*SA_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SA_rate[i,j,k]
S_death[,,] <- if(S[i,j,k]*dt*eta < 0) 0 else S[i,j,k]*dt*eta
S_age[,,] <- if(S[i,j,k]*dt*age_rate[i] < 0) 0 else S[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM T #####
dim(TP_trans) <- c(na,nh,num_int)
dim(T_death) <- c(na,nh,num_int)
dim(T_age) <- c(na,nh,num_int)

TP_trans[,,] <- if(T[i,j,k]*dt*rT < 0) 0 else T[i,j,k]*dt*rT
T_death[,,] <- if(T[i,j,k]*dt*eta < 0) 0 else T[i,j,k]*dt*eta
T_age[,,] <- if(T[i,j,k]*dt*age_rate[i] < 0) 0 else T[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM D #####
dim(DA_trans) <- c(na,nh,num_int)
dim(DS_trans) <- c(na,nh,num_int)
dim(D_death) <- c(na,nh,num_int)
dim(D_age) <- c(na,nh,num_int)

DA_trans[,,] <- if(D[i,j,k]*dt*rD < 0) 0 else D[i,j,k]*dt*rD
DS_trans[,,] <- if(D[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else D[i,j,k]*alpha_smc_array[i,j,k]
D_death[,,] <- if(D[i,j,k]*dt*eta < 0) 0 else D[i,j,k]*dt*eta
D_age[,,] <- if(D[i,j,k]*dt*age_rate[i] < 0) 0 else D[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM A #####
#Rates
dim(AT_rate) <- c(na,nh,num_int)
dim(AD_rate) <- c(na,nh,num_int)

AT_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * ft
AD_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * (1-ft)

# Number of humans transitioning per timestep
dim(AU_trans) <- c(na,nh,num_int)
dim(AT_trans) <- c(na,nh,num_int)
dim(AD_trans) <- c(na,nh,num_int)
dim(AS_trans) <- c(na,nh,num_int) #SMC Infection Clearing
dim(A_age) <- c(na,nh,num_int)
dim(A_death) <- c(na,nh,num_int)

AU_trans[,,] <- if(A[i,j,k]*dt*rA < 0) 0 else A[i,j,k]*dt*rA
AT_trans[,,] <- if(A[i,j,k]*dt*AT_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AT_rate[i,j,k]
AD_trans[,,] <- if(A[i,j,k]*dt*AD_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AD_rate[i,j,k]
AS_trans[,,] <- if(A[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else A[i,j,k]*alpha_smc_array[i,j,k] #smc Infection Clearing
A_death[,,] <- if(A[i,j,k]*dt*eta < 0) 0 else A[i,j,k]*dt*eta
A_age[,,] <- if(A[i,j,k]*dt*age_rate[i] < 0) 0 else A[i,j,k]*dt*age_rate[i]

#### TRANSITIONS FROM U #####
#Rates
dim(UA_rate) <- c(na,nh,num_int)
dim(UD_rate) <- c(na,nh,num_int)
dim(UT_rate) <- c(na,nh,num_int)

UA_rate[,,] <- (1-phi[i,j,k]) * FOI_smc[i,j,k]
UD_rate[,,] <- phi[i,j,k] * (1-ft) * FOI_smc[i,j,k]
UT_rate[,,] <- phi[i,j,k] * ft * FOI_smc[i,j,k]

# Number of humans transitioning per timestep
dim(UA_trans) <- c(na,nh,num_int)
dim(UD_trans) <- c(na,nh,num_int)
dim(UT_trans) <- c(na,nh,num_int)
dim(US_trans) <- c(na,nh,num_int)
dim(US_trans_SMC) <- c(na,nh,num_int)
dim(U_age) <- c(na,nh,num_int)
dim(U_death) <- c(na,nh,num_int)

US_trans[,,] <- if(U[i,j,k]*dt*rU < 0) 0 else U[i,j,k]*dt*rU
US_trans_SMC[,,] <- if(U[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else U[i,j,k]*alpha_smc_array[i,j,k]
UA_trans[,,] <- if(U[i,j,k]*dt*UA_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UA_rate[i,j,k]
UD_trans[,,] <- if(U[i,j,k]*dt*UD_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UD_rate[i,j,k]
UT_trans[,,] <- if(U[i,j,k]*dt*UT_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UT_rate[i,j,k]
U_age[,,] <- if(U[i,j,k]*dt*age_rate[i] < 0) 0 else U[i,j,k]*dt*age_rate[i]
U_death[,,] <- if(U[i,j,k]*dt*eta < 0) 0 else U[i,j,k]*dt*eta

#### TRANSITIONS FROM P #####
dim(PS_trans) <- c(na,nh,num_int)
dim(P_death) <- c(na,nh,num_int)
dim(P_age) <- c(na,nh,num_int)

PS_trans[,,] <-if(P[i,j,k]*dt*rP < 0) 0 else P[i,j,k]*dt*rP
P_death[,,] <- if(P[i,j,k]*dt*eta < 0) 0 else P[i,j,k]*dt*eta
P_age[,,] <- if(P[i,j,k]*dt*age_rate[i] < 0) 0 else P[i,j,k]*dt*age_rate[i]

#------------------------------ TRANSITIONS ----------------------------------
update(S[1,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_SMC[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + births[1,j,k]
update(S[2:na,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_SMC[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + S_age[i-1,j,k]

update(T[1,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k]
update(T[2:na,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k] + T_age[i-1,j,k]

update(D[1,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k]
update(D[2:na,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k] + D_age[i-1,j,k]

update(A[1,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k]
update(A[2:na,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k] + A_age[i-1,j,k]

update(U[1,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_SMC[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k]
update(U[2:na,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_SMC[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k] + U_age[i-1,j,k]

update(P[1,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k]
update(P[2:na,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k] + P_age[i-1,j,k]

##------------------------------------------------------------------------------
##########################################
## KEY INFECTION PARAMETERS (FOI & EIR) ##
##########################################
##------------------------------------------------------------------------------

#----------FORCE OF INFECTION-----------------
# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- parameter()
lag_rates <- parameter(type = "integer")

dim(FOI_eq) <- c(na,nh)
FOI_eq <- parameter()

dim(FOI_XL) <- c(na,nh,num_int,lag_rates)
initial(FOI_XL[,,,]) <- FOI_eq[i,j]

update(FOI_XL[,,,1]) <- FOI_XL[i,j,k,1] + dt*((lag_rates/dE)*FOI_lag[i,j,k] - (lag_rates/dE)*FOI_XL[i,j,k,1])
update(FOI_XL[,,,2:lag_rates]) <- FOI_XL[i,j,k,l] + dt*((lag_rates/dE)*FOI_XL[i,j,k,l-1] - (lag_rates/dE)*FOI_XL[i,j,k,l])

dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- FOI_XL[i,j,k,lag_rates]


#------------ EIR -----------------
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age <- parameter()
dim(rel_foi) <- nh
rel_foi <- parameter()

dim(EIR_species1) <- c(na,nh,num_int)
EIR_species1[,,] <- av_human1[k] * rel_foi[j] * foi_age[i] * Iv1/omega

dim(EIR_species2) <- c(na,nh,num_int)
EIR_species2[,,] <- av_human2[k] * rel_foi[j] * foi_age[i] * Iv2/omega

dim(EIR) <- c(na,nh,num_int)
EIR[,,] <- EIR_species1[i,j,k] + EIR_species2[i,j,k]

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

#----------------- DEFINE IMMUNITY PARAMETERS -------------------
dCM <- parameter()
uCA <- parameter()
dCA <- parameter()
dB <- parameter()
uB <- parameter()
dID <- parameter()
uD <- parameter()
x_I <- parameter()
dim(x_I) <- na
age20l <- parameter(type = "integer")
age20u <- parameter(type = "integer")
age_20_factor <- parameter()
PM <- parameter()

#---------- UPDATE IMMUNITY STATES FOR EACH HUMAN COMPARTMENTS ---------

# ICM - maternally acquired immunity
init_ICM <- parameter()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)
dim(init_ICM_pre) <- c(nh,num_int)
init_ICM_pre[1:nh,1:num_int] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))

update(ICM[1, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i])
update(ICM[2:na, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i])

# ICA - exposure driven immunity
init_ICA <- parameter()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

update(ICA[1, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i])
update(ICA[2:na, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i])

dim(IC) <- c(na,nh,num_int)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# IB - infection blocking immunity
init_IB <- parameter()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

update(IB[1, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i])
update(IB[2:na, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i])

# ID - detection immunity
init_ID <- parameter()
dim(init_ID) <- c(na,nh,num_int)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,num_int)

update(ID[1, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i])
update(ID[2:na, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i])

#----------- PROBABILITIES INFERRED FROM IMMUNITY --------------------
# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- parameter()
phi1 <- parameter()
IC0 <- parameter()
kC <- parameter()
dim(phi) <- c(na,nh,num_int)
phi[1:na,1:nh,1:num_int] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- parameter()
b1 <- parameter()
kB <- parameter()
IB0 <- parameter()
dim(b) <- c(na,nh,num_int)
b[1:na, 1:nh, 1:num_int] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)


# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- parameter()
fD0 <- parameter()
gammaD <- parameter()
d1 <- parameter()
ID0 <- parameter()
kD <- parameter()
dim(age_vector) <- na
age_vector <- parameter()

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age_vector[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Pv - Latently infected (Pre-Infectious?) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- parameter()
init_Pv <- parameter()
init_Iv <- parameter()
initial(Sv) <- init_Sv * mv0
initial(Pv) <- init_Pv * mv0
initial(Iv) <- init_Iv * mv0
initial(Sv1) <- init_Sv * mv0
initial(Pv1) <- init_Ev * mv0
initial(Iv1) <- init_Iv * mv0
initial(Sv2) <- init_Sv * mv0 * 0.01 ## This is a bit of a fudge but initial low carrying capacity for stephensi will bring it even lower during equilibration time.
initial(Pv2) <- init_Ev * mv0 * 0.01 ## Still to do: equilibrium solution for two (or N) mosquitoes and updating this to reflect that.
initial(Iv2) <- init_Iv * mv0 * 0.01

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- parameter()
cD <- parameter()
cT <- parameter()
gamma1 <- parameter()
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
lag_ratesMos <- parameter(type = "integer")

FOIv_eq <- parameter()
dim(FOIv_species1) <- lag_ratesMos
initial(FOIv_species1[]) <- FOIv_eq*delayGam/lag_ratesMos
dim(FOIv_species2) <- lag_ratesMos
initial(FOIv_species2[]) <- 0.01 * FOIv_eq*delayGam/lag_ratesMos

update(FOIv_species1[1]) <- FOIv_species1[1] + dt*(lag_FOIv_species1 - (lag_ratesMos/delayGam)*FOIv_species1[1])
update(FOIv_species1[2:lag_ratesMos]) <- FOIv_species1[i] + dt*((lag_ratesMos/delayGam)*FOIv_species1[i-1] - (lag_ratesMos/delayGam)*FOIv_species1[i])
update(FOIv_species2[1]) <- FOIv_species2[1] + dt*(lag_FOIv_species2 - (lag_ratesMos/delayGam)*FOIv_species2[1])
update(FOIv_species2[2:lag_ratesMos]) <- FOIv_species2[i] + dt*((lag_ratesMos/delayGam)*FOIv_species2[i-1] - (lag_ratesMos/delayGam)*FOIv_species2[i])

dim(FOIvijk_species1) <- c(na,nh,num_int)
dim(FOIvijk_species2) <- c(na,nh,num_int)
omega <- parameter()
FOIvijk_species1[1:na, 1:nh, 1:num_int] <- ((cT*smc_rel_c_mask[i,j,k]*T[i,j,k] + cD*smc_rel_c_mask[i,j,k]*D[i,j,k] + cA[i,j,k]*smc_rel_c_mask[i,j,k]*A[i,j,k] + cU*smc_rel_c_mask[i,j,k]*U[i,j,k])/H) * rel_foi[j] * av_mosq1[k]*foi_age[i]/omega ## For discrete human compartments
FOIvijk_species2[1:na, 1:nh, 1:num_int] <- ((cT*smc_rel_c_mask[i,j,k]*T[i,j,k] + cD*smc_rel_c_mask[i,j,k]*D[i,j,k] + cA[i,j,k]*smc_rel_c_mask[i,j,k]*A[i,j,k] + cU*smc_rel_c_mask[i,j,k]*U[i,j,k])/H) * rel_foi[j] * av_mosq2[k]*foi_age[i]/omega ## For discrete human compartments

lag_FOIv_species1=sum(FOIvijk_species1)
lag_FOIv_species2=sum(FOIvijk_species2)

ince1 <- FOIv_species1[lag_ratesMos] * lag_ratesMos/delayGam * Sv1
initial(ince_delay1[]) <- FOIv_eq*init_Sv*mv0*delayMos_use/lag_ratesMos
dim(ince_delay1) <- lag_ratesMos

update(ince_delay1[1]) <- ince_delay1[1] + dt*(ince1 - (lag_ratesMos/delayMos_use)*ince_delay1[1])
update(ince_delay1[2:lag_ratesMos]) <- ince_delay1[i] + dt*((lag_ratesMos/delayMos_use)*ince_delay1[i-1] - (lag_ratesMos/delayMos_use)*ince_delay1[i])
incv1 <- ince_delay1[lag_ratesMos]*lag_ratesMos/delayMos_use * surv1

ince2 <- FOIv_species2[lag_ratesMos] * lag_ratesMos/delayGam * Sv2
initial(ince_delay2[]) <- FOIv_eq*init_Sv*mv0*delayMos_use/lag_ratesMos
dim(ince_delay2) <- lag_ratesMos

update(ince_delay2[1]) <- ince_delay2[1] + dt*(ince2 - (lag_ratesMos/delayMos_use)*ince_delay2[1])
update(ince_delay2[2:lag_ratesMos]) <- ince_delay2[i] + dt*((lag_ratesMos/delayMos_use)*ince_delay2[i-1] - (lag_ratesMos/delayMos_use)*ince_delay2[i])
incv2 <- ince_delay2[lag_ratesMos]*lag_ratesMos/delayMos_use * surv2

## Outputs for checking (new two species model)
# initial(ince1_out) <- 0
# update(ince1_out) <- ince1
# initial(incv1_out) <- 0
# update(incv1_out) <- incv1
# initial(ince2_out) <- 0
# update(ince2_out) <- ince2
# initial(incv2_out) <- 0
# update(incv2_out) <- incv2

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- parameter()
delayMos <- parameter()
delayMos_use <- delayMos
# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos_use)

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa1 <- 0.5*PL1/dPL
betaa2 <- 0.5*PL2/dPL

update(Sv1) <- if(Sv1 + dt*(-ince1 - mu1*Sv1 + betaa1) < 0) 0 else Sv1 + dt*(-ince1 - mu1*Sv1 + betaa1)
update(Pv1) <- if(Pv1 + dt*(ince1 - incv1 - mu1*Pv1) < 0) 0 else Pv1 + dt*(ince1 - incv1 - mu1*Pv1)
update(Iv1) <- if(Iv1 + dt*(incv1 - mu1*Iv1) < 0) 0 else Iv1 + dt*(incv1 - mu1*Iv1)

update(Sv2) <- if(Sv2 + dt*(-ince2 - mu2*Sv2 + betaa2) < 0) 0 else Sv2 + dt*(-ince2 - mu2*Sv2 + betaa2)
update(Pv2) <- if(Pv2 + dt*(ince2 - incv2 - mu2*Pv2) < 0) 0 else Pv2 + dt*(ince2 - incv2 - mu2*Pv2)
update(Iv2) <- if(Iv2 + dt*(incv2 - mu2*Iv2) < 0) 0 else Iv2 + dt*(incv2 - mu2*Iv2)

# Total mosquito population
initial(mv1) <- 0
update(mv1) <- Sv1 + Pv1 + Iv1
initial(mv2) <- 0
update(mv2) <- Sv2 + Pv2 + Iv2
mv <- mv1 + mv2 ## CHECK: Does this need to be an initial/update type of process instead? Are we going to be out-by-one doing it this way?

##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dLL <- parameter()
dPL <- parameter()
dEL <- parameter()
muLL <- parameter()
muEL <- parameter()
muPL <- parameter()
gammaL <- parameter()

# fitted entomological parameters:
mv0 <- parameter()
mu0 <- parameter()
tau1 <- parameter()
tau2 <- parameter()
betaL <- parameter()

# Entomological variables:
mu0_use <- mu0
p10 <- exp(-mu0_use * tau1)  # probability of surviving one feeding cycle
p2 <- exp(-mu0_use * tau2)  # probability of surviving one resting cycle

eov1 <- betaL/mu1*(exp(mu1/fv1)-1)
beta_larval1 <- eov1*mu1*exp(-mu1/fv1)/(1-exp(-mu1/fv1)) # Number of eggs laid per day
eov2 <- betaL/mu2*(exp(mu2/fv2)-1)
beta_larval2 <- eov2*mu2*exp(-mu2/fv2)/(1-exp(-mu2/fv2)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)

lambda_species1 <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval1*muLL*dEL/(2*muEL*mu0_use*dLL*(1+dPL*muPL)))
lambda_species2 <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval2*muLL*dEL/(2*muEL*mu0_use*dLL*(1+dPL*muPL)))

### CHECK: THIS IS WHAT WAS THERE PREVIOUSLY, AND IT HAS mv0 IN THE NUMERATOR, WHEREAS MINE DOESN'T (AND POSS ICDMM DIDN'T)
# K0 <- 2*mv0*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)
K0_species1 <- if(as.integer(step) == 0) 2*density_vec_sp1[as.integer(1)]*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda_species1+1)/(lambda_species1/(muLL*dEL)-1/(muLL*dLL)-1) else
  2*density_vec_sp1[as.integer(step)]*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda_species1+1)/(lambda_species1/(muLL*dEL)-1/(muLL*dLL)-1)
K0_species2 <- if(as.integer(step) == 0) 2*density_vec_sp2[as.integer(1)]*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda_species2+1)/(lambda_species2/(muLL*dEL)-1/(muLL*dLL)-1) else
  2*density_vec_sp2[as.integer(step)]*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda_species2+1)/(lambda_species2/(muLL*dEL)-1/(muLL*dLL)-1)
## density_vec is defined above and is a monotonically increasing function designed to mimic and simulate invasion and establishment of An. stephensi

# parameters for species 1 density and species 2 density (latter increasing in abundance over time - NOTE this is distinct and on top of the custom seasonality)
time_length <- parameter()
density_vec_sp1[] <- parameter()
dim(density_vec_sp1) <- time_length # CHECK AND MAKE THIS SOM FUNCTION OF ALREADY EXISTING INPUTS
density_vec_sp2[] <- parameter()
dim(density_vec_sp2) <- time_length # CHECK AND MAKE THIS SOM FUNCTION OF ALREADY EXISTING INPUTS

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
theta_species1 <- interpolate(days, daily_rain_input, "linear")
KL_species1 <- K0_species1 * theta_species1

## Custom seasonality function for species 2 (which is An.stephensi for us)
custom_seasonality[] <- parameter()
time_length <- parameter(type = "integer")
dim(custom_seasonality) <- time_length # CHECK: NEED TO CHANGE THIS AS I SUSPECT THIS CAN BE SOME FUNCTION OF ALREADY INPUTTED PARAMETERS
theta_species2 <- if(as.integer(step) == 0) custom_seasonality[as.integer(1)] else custom_seasonality[as.integer(step)]
KL_species2 <- K0_species2 * theta_species2
initial(theta_species2_out) <- 0
update(theta_species2_out) <- theta_species2

fv1 <- 1/( tau1/(1-zbar1) + tau2 ) # mosquito feeding rate for species 1 (zbar is a derived intervention param.)
fv2 <- 1/( tau1/(1-zbar2) + tau2 ) # mosquito feeding rate species 2 (zbar is a derived intervention param.)
mu1 <- -fv1*log(p1_1*p2) # mosquito death rate - species 1
mu2 <- -fv2*log(p1_2*p2) # mosquito death rate - species 2

# finding equilibrium and initial values for EL, LL & PL
init_PL <- parameter()
initial(PL) <- init_PL
init_LL <- parameter()
initial(LL) <- init_LL
init_EL <- parameter()
initial(EL) <- init_EL

## This is a bit of a fudge but initial low carrying capacity for stephensi will bring it even lower during equilibration time.
## Still to do: equilibrium solution for two (or N) mosquito species and updating this to reflect that.
initial(PL1) <- init_PL
initial(LL1) <- init_LL
initial(EL1) <- init_EL
initial(PL2) <- init_PL * 0.01
initial(LL2) <- init_LL * 0.01
initial(EL2) <- init_EL * 0.01

## For Species 1
# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
update(EL1) <- if(EL1 + dt*(beta_larval1*mv1-muEL*(1+(EL1+LL1)/KL_species1)*EL1 - EL1/dEL) < 0) 0 else (EL1 + dt*(beta_larval1*mv1-muEL*(1+(EL1+LL1)/KL_species1)*EL1 - EL1/dEL))
# egg hatching - den. dep. mortality - maturing larvae
update(LL1) <- if(LL1 + dt*(EL1/dEL - muLL*(1+gammaL*(EL1 + LL1)/KL_species1)*LL1 - LL1/dLL) < 0) 0 else (LL1 + dt*(EL1/dEL - muLL*(1+gammaL*(EL1 + LL1)/KL_species1)*LL1 - LL1/dLL))
# pupae - mortality - fully developed pupae
update(PL1) <- if(PL1 + dt*(LL1/dLL - muPL*PL1 - PL1/dPL) < 0) 0 else (PL1 + dt*(LL1/dLL - muPL*PL1 - PL1/dPL))

## For Species 2
# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
update(EL2) <- if(EL2 + dt*(beta_larval2*mv2-muEL*(1+(EL2+LL2)/KL_species2)*EL2 - EL2/dEL) < 0) 0 else (EL2 + dt*(beta_larval2*mv2-muEL*(1+(EL2+LL2)/KL_species2)*EL2 - EL2/dEL))
# egg hatching - den. dep. mortality - maturing larvae
update(LL2) <- if(LL2 + dt*(EL2/dEL - muLL*(1+gammaL*(EL2 + LL2)/KL_species2)*LL2 - LL2/dLL) < 0) 0 else (LL2 + dt*(EL2/dEL - muLL*(1+gammaL*(EL2 + LL2)/KL_species2)*LL2 - LL2/dLL))
# pupae - mortality - fully developed pupae
update(PL2) <- if(PL2 + dt*(LL2/dLL - muPL*PL2 - PL2/dPL) < 0) 0 else (PL2 + dt*(LL2/dLL - muPL*PL2 - PL2/dPL))

##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

################## SMC #######################
#See supplementary materials of [Thompson, 2022] - https://doi.org/10.1016/S2214-109X(22)00416-8
max_smc_cov <- parameter()

##Parameters relevant to clearance of existing infections by SMC treatment
dim(alpha_smc) <- n_ts + 1
alpha_smc <- parameter()
alpha_smc_timestep <- interpolate(iterations, alpha_smc, "constant")

dim(alpha_smc_array) <-  c(na,nh,num_int)
alpha_smc_array[,,] <- smc_mask[i,j,k] * alpha_smc_timestep #Multiply by proportion of individuals in SMC compartment currently receiving SMC

##Parameters relating to prophylaxis effect of SMC
dim(P_smc_daily) <- n_days + 1
P_smc_daily <- parameter()
P_smc <- interpolate(days, P_smc_daily, "linear")

dim(smc_mask) <- c(na,nh,num_int)
smc_mask <- parameter()
dim(FOI_smc) <- c(na,nh,num_int)
FOI_smc[,,] <- FOI[i,j,k] * (1 - (smc_mask[i,j,k] * P_smc))

##Parameters relating to reduced infectivity of existing infections upon SMC dose
dim(rel_c_days) <- n_days + 1
rel_c_days <- parameter()

rel_c <- interpolate(days, rel_c_days, "constant")
dim(smc_rel_c_mask) <- c(na,nh,num_int)
smc_rel_c_mask[,,] <- 1 - (smc_mask[i,j,k] * (1-rel_c)) #Value = SMC_rel_c for those scheduled to receive SMC on rel_c_days. Else = 1.


################## ITN #######################
# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
max_itn_cov <- parameter()

Q0_species1 <- parameter() # species 1 proportion of anthropophagy
Q0_species2 <- parameter() # species 2 proportion of anthropophagy
bites_Bed_species1 <- parameter() # endophagy in bed
bites_Bed_species2 <- parameter() # endophagy in bed
chi_species1 <- parameter() # proportion of vector endophily
chi_species2 <- parameter() # proportion of vector endophily

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN
# d - probability of dying after hitting ITN
# s - probability of successful feed after hitting ITN

# The maximum (and then minimum) r and d values for ITN on day 0 before they decay
rn_1 <- parameter()
rn_2 <- parameter()
dn0_1 <- parameter()
dn0_2 <- parameter()
rnm_1 <- parameter()
rnm_2 <- parameter()

dim(itn_decay_daily) <- n_days + 1 #Bednet insecticide decay
dim(itn_eff_cov_daily) <- n_days + 1 #Proportion of individuals in the ITN compartment who currently have ITNs.= daily_ITN / max(daily_ITN)
itn_decay_daily <- parameter()
itn_eff_cov_daily <- parameter()

itn_eff_cov <- interpolate(days, itn_eff_cov_daily, "linear")
itn_decay <- interpolate(days, itn_decay_daily, "linear")

d_itn_1 <- dn0_1*itn_decay*itn_eff_cov
r_itn_1 <- (rnm_1 + (rn_2 - rnm_1)*itn_decay)*itn_eff_cov #Bednet repellency does not decay to zero.
s_itn_1 <- 1 - d_itn_1 - r_itn_1

d_itn_2 <- dn0_2*itn_decay*itn_eff_cov
r_itn_2 <- (rnm_2 + (rn_2 - rnm_2)*itn_decay)*itn_eff_cov #Bednet repellency does not decay to zero.
s_itn_2 <- 1 - d_itn_2 - r_itn_2

################## GENERAL INTERVENTION PARAMETERS #######################
num_int <- parameter()

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-max_itn_cov)*(1-max_smc_cov)  # {No intervention}
cov_[2] <- max_itn_cov*(1-max_smc_cov) # 	   {ITN only}
cov_[3] <- (1-max_itn_cov)*max_smc_cov	#      {SMC only}
cov_[4] <- max_itn_cov*max_smc_cov #	   {Both itn and SMC}
cov[] <- cov_[i]
dim(cov) <- num_int

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - bites_Bed + bites_Bed*s_itn
w_[3] <- 1
w_[4] <- 1 - bites_Bed + bites_Bed*s_itn
w[] <- w_[i]
dim(w) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- bites_Bed*r_itn
z_[3] <- 0
z_[4] <- bites_Bed*r_itn
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- sum(zhi)
wh <- sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.


##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------
dim(clin_inc) <- c(na,nh,num_int)
clin_inc[,,] <- ST_trans[i,j,k] + SD_trans[i,j,k] + AT_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k]  + UT_trans[i,j,k]

prev_dim <- parameter()
dim(min_age_prev) <- prev_dim
dim(max_age_prev) <- prev_dim
min_age_prev <- parameter(type = "integer")
max_age_prev <- parameter(type = "integer")

dim(n_prev) <- prev_dim

n_prev[1:prev_dim] <- sum(S[min_age_prev[i]:max_age_prev[i],,]) + sum(T[min_age_prev[i]:max_age_prev[i],,]) + sum(D[min_age_prev[i]:max_age_prev[i],,]) +
  sum(A[min_age_prev[i]:max_age_prev[i],,]) + sum(U[min_age_prev[i]:max_age_prev[i],,]) + sum(P[min_age_prev[i]:max_age_prev[i],,])

dim(detect_prev_full) <- c(na,nh,num_int)
detect_prev_full[,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]

dim(detect_prev) <- prev_dim
detect_prev[1:prev_dim] <- sum(detect_prev_full[min_age_prev[i]:max_age_prev[i],,])

##Proportion of population in user defined prevalence age groups
dim(n_ud_prev) <- prev_dim
initial(n_ud_prev[]) <- min_age_prev[i] #user defined prevalence
update(n_ud_prev[]) <- n_prev[i]

##Proportion of population in user defined prevalence age groups with detectable malaria
dim(n_ud_detect_prev) <- prev_dim
initial(n_ud_detect_prev[]) <- min_age_prev[i]
update(n_ud_detect_prev[]) <- detect_prev[i]

# #Incidence
inc_dim <- parameter()
dim(min_age_inc) <- inc_dim
dim(max_age_inc) <- inc_dim
min_age_inc <- parameter(type = "integer")
max_age_inc <- parameter(type = "integer")

dim(n_ud_inc) <- inc_dim
initial(n_ud_inc[]) <- min_age_inc[i]
update(n_ud_inc[]) <- sum(clin_inc[min_age_inc[i]:max_age_inc[i],,]) / dt

Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph

initial(S_count) <- 0
initial(T_count) <- 0
initial(D_count) <- 0
initial(A_count) <- 0
initial(U_count) <- 0
initial(P_count) <- 0

update(S_count) <- Sh
update(T_count) <- Th
update(D_count) <- Dh
update(A_count) <- Ah
update(U_count) <- Uh
update(P_count) <- Ph

dim(all) <- c(na,nh,num_int)
all[,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]
##-------------------- IMMUNITY OUTPUTS ----------------------------------
dim(icm_pop) <- c(na,nh,num_int)
icm_pop[,,] <- all[i,j,k] * ICM[i,j,k]
initial(icm_mean) <- 0
update(icm_mean) <- sum(icm_pop[,,]) / H

dim(ica_pop) <- c(na,nh,num_int)
ica_pop[,,] <- all[i,j,k] * ICA[i,j,k]
initial(ica_mean) <- 0
update(ica_mean) <- sum(ica_pop[,,]) / H

dim(id_pop) <- c(na,nh,num_int)
id_pop[,,] <- all[i,j,k] * ID[i,j,k]
initial(id_mean) <- 0
update(id_mean) <- sum(id_pop[,,]) / H

dim(ib_pop) <- c(na,nh,num_int)
ib_pop[,,] <- all[i,j,k] * IB[i,j,k]
initial(ib_mean) <- 0
update(ib_mean) <- sum(ib_pop[,,]) / H

dim(ic_pop) <- c(na,nh,num_int)
ic_pop[,,] <- all[i,j,k] * IC[i,j,k]
initial(ic_mean) <- 0
update(ic_mean) <- sum(ic_pop[,,]) / H

##-------------------- OTHER OUTPUTS ----------------------------------
dim(all_deaths) <- c(na,nh,num_int)
all_deaths[,,] <- S_death[i,j,k] + T_death[i,j,k] + D_death[i,j,k] + A_death[i,j,k] + U_death[i,j,k] + P_death[i,j,k]
initial(natural_deaths) <- 0
update(natural_deaths) <- sum(all_deaths[,,])

dim(EIR_pop) <- c(na,nh,num_int)
EIR_pop[,,] <- EIR[i,j,k] * all[i,j,k]
initial(EIR_mean) <- 0
update(EIR_mean) <- sum(EIR_pop[,,]) / H

initial(mu_mosq) <- 0
update(mu_mosq) <- mu
