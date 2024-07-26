library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(lubridate)
library(R2jags)
library(doParallel)

ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  
setwd("/Users/dmorina/Downloads/")
df <- read_excel("df_monthly_med.xlsx", guess_max = 10000)
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")

df <- df %>% distinct()

df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals), names_sep = "_")

### Select time period with complete information (2014-2021)
df2 <- df2[df2$date>="2014-01-01" & df2$date <= "2021-06-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)


df_dates <- readRDS("df_dates.RDS")
df_dates$study_period <- as.Date(df_dates$study_period, format="%Y-%m-%d")
df_dates <- df_dates %>% 
  group_by(month = lubridate::floor_date(study_period, "month")) %>%
  summarize(SAR_SeaWatch123=sum(SAR_SeaWatch123), SAR_SeaWatch4=sum(SAR_SeaWatch4),
            SAR_Lifeline=sum(SAR_Lifeline), SAR_SeaEye_TheSeaEye=sum(SAR_SeaEye_TheSeaEye),
            SAR_SeaEye_TheSeefuchs=sum(SAR_SeaEye_TheSeefuchs), SAR_SeaEye_AlanKurdi=sum(SAR_SeaEye_AlanKurdi),
            SAR_OpenArms_Astral=sum(SAR_OpenArms_Astral), SAR_OpenArms_GolfoAzzuro=sum(SAR_OpenArms_GolfoAzzuro),
            SAR_OpenArms_OpenArms=sum(SAR_OpenArms_OpenArms), SAR_Mediterranea=sum(SAR_Mediterranea),
            SAR_SMH=sum(SAR_SMH), SAR_LouiseMichel=sum(SAR_LouiseMichel), SAR_RefugeeRescue=sum(SAR_RefugeeRescue),
            SAR_MOAS=sum(SAR_MOAS), SAR_JugendRettet=sum(SAR_JugendRettet), SAR_MSFandSOS=sum(SAR_MSFandSOS),
            SAR_SavetheChildren=sum(SAR_SavetheChildren), SAR_MSF_BourbonArgos=sum(SAR_MSF_BourbonArgos),
            SAR_MSF_Dignity1=sum(SAR_MSF_Dignity1), SAR_MSF_VosPrudence=sum(SAR_MSF_VosPrudence),
            SAR_Resqship=sum(SAR_Resqship), SAR_Lifeboat=sum(SAR_Lifeboat), FRONTEX_op_POSEIDON=sum(FRONTEX_op_POSEIDON),
            FRONTEX_op_HERA=sum(FRONTEX_op_HERA), FRONTEX_op_MINERVA=sum(FRONTEX_op_MINERVA), 
            FRONTEX_op_TRITON=sum(FRONTEX_op_TRITON), FRONTEX_op_THEMIS=sum(FRONTEX_op_THEMIS),
            EUNAVFOR_SOPHIA=sum(EUNAVFOR_SOPHIA), EUNAVFOR_IRINI=sum(EUNAVFOR_IRINI), MARE_NOSTRUM=sum(MARE_NOSTRUM),
            COASTGUARD_LIBYA=sum(COASTGUARD_LIBYA), Extension_SAR_LIBYA=sum(Extension_SAR_LIBYA))

df_dates$SAR_Total <- rowSums(df_dates[, 2:23])
# df_dates$TotalDays <- days_in_month(df_dates$month)
# df_dates$SAR_Total <- df_dates$SAR_Total/df_dates$TotalDays
df_dates$FrontiersControl <- rowSums(df_dates[, 24:33])
df_dates <- df_dates[as.Date(df_dates$month)>="2014-01-01" & as.Date(df_dates$month)<="2021-06-01", ]

meteo <- read.csv(file = "df_monthly_deaths_share.csv", sep = ";")
meteo$date <- as.Date(meteo$date, format = "%Y-%m-%d")
meteo <- meteo[!duplicated(meteo$date), ]
meteo <- meteo[meteo$date>="2014-01-01" & meteo$date<="2021-06-01", ]



### MODEL AMB q_b1 * At[t]

# JAGS to estimate parameters
ar1_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0)/(1-phi), (sigma2_w/(1-phi^2))^-1)
  for (t in 2:N) {
    x[t] ~ dnorm((delta_b0) + phi*x[t-1], (sigma2_w/(1-phi^2))^-1)
  }
  for (t in 1:N)
  {
    logit(q[t]) <- q_b0+q_b1*At[t]
    logit(omega[t]) <- w_b0+w_b1*At[t]
    mu[1, t] <- (delta_b0)/(1-phi) 
    mu[2, t] <- q[t]*(delta_b0)/(1-phi)
    tau.m[1, t] <- (sigma2_w/(1-phi^2))^-1
    tau.m[2, t] <- (1/(q[t]^2))*(sigma2_w/(1-phi^2))^-1
    y[t] ~ dnormmix(mu[, t], tau.m[, t], p.m[, t])
    p.m[1, t] <- 1-omega[t]
    p.m[2, t] <- omega[t]
    component_chosen[t] ~ dcat(p.m[, t])
  }
  # priors
  q_b0 ~ dnorm(0, 0.01)
  q_b1 ~ dnorm(0, 0.01)
  w_b0 ~ dnorm(0, 0.01)
  w_b1 ~ dnorm(0, 0.01)
  delta_b0 ~ dnorm(0, 1/1000)
  phi ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR), At=df2$At_CMR)

params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_At <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar1_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))

save(list="cmed_route_At", file="Resultats/cmed_route_covs_At.RData")



### MODEL AMB q_b1 * At[t] i delta_b1 * FC

# JAGS to estimate parameters
ar1_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0+delta_b1*FC[1])/(1-phi), (sigma2_w/(1-phi^2))^-1)
  for (t in 2:N) {
    x[t] ~ dnorm((delta_b0+delta_b1*FC[t]) + phi*x[t-1], (sigma2_w/(1-phi^2))^-1)
  }
  for (t in 1:N)
  {
    logit(q[t]) <- q_b0+q_b1*At[t]
    logit(omega[t]) <- w_b0+w_b1*At[t]
    mu[1, t] <- (delta_b0+delta_b1*FC[t])/(1-phi) 
    mu[2, t] <- q[t]*(delta_b0+delta_b1*FC[t])/(1-phi)
    tau.m[1, t] <- (sigma2_w/(1-phi^2))^-1
    tau.m[2, t] <- (1/(q[t]^2))*(sigma2_w/(1-phi^2))^-1
    y[t] ~ dnormmix(mu[, t], tau.m[, t], p.m[, t])
    p.m[1, t] <- 1-omega[t]
    p.m[2, t] <- omega[t]
    component_chosen[t] ~ dcat(p.m[, t])
  }
  # priors
  q_b0 ~ dnorm(0, 0.01)
  q_b1 ~ dnorm(0, 0.01)
  w_b0 ~ dnorm(0, 0.01)
  w_b1 ~ dnorm(0, 0.01)
  delta_b0 ~ dnorm(0, 1/1000)
  delta_b1 ~ dnorm(0, 1/1000)
  phi ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}
### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR), FC=log(df_dates$FrontiersControl), At=df2$At_CMR)

params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "delta_b1", "phi", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_At_FC <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar1_model,
                                           n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))

save(list="cmed_route_At_FC", file="Resultats/cmed_route_covs_At_FC.RData")



### MODEL AMB q_b1 * At[t] + q_b2 * tempMa i delta_b1 * FC

ar1_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0+delta_b1*FC[1])/(1-phi), (sigma2_w/(1-phi^2))^-1)
  for (t in 2:N) {
    x[t] ~ dnorm((delta_b0+delta_b1*FC[t]) + phi*x[t-1], (sigma2_w/(1-phi^2))^-1)
  }
  for (t in 1:N)
  {
    logit(q[t]) <- q_b0+q_b1*At[t]+q_b2*tempMa[t]
    logit(omega[t]) <- w_b0+w_b1*At[t]+w_b2*tempMa[t]
    mu[1, t] <- (delta_b0+delta_b1*FC[t])/(1-phi) 
    mu[2, t] <- q[t]*(delta_b0+delta_b1*FC[t])/(1-phi)
    tau.m[1, t] <- (sigma2_w/(1-phi^2))^-1
    tau.m[2, t] <- (1/(q[t]^2))*(sigma2_w/(1-phi^2))^-1
    y[t] ~ dnormmix(mu[, t], tau.m[, t], p.m[, t])
    p.m[1, t] <- 1-omega[t]
    p.m[2, t] <- omega[t]
    component_chosen[t] ~ dcat(p.m[, t])
  }
  # priors
  q_b0 ~ dnorm(0, 0.01)
  q_b1 ~ dnorm(0, 0.01)
  q_b2 ~ dnorm(0, 0.01)
  w_b0 ~ dnorm(0, 0.01)
  w_b1 ~ dnorm(0, 0.01)
  w_b2 ~ dnorm(0, 0.01)
  delta_b0 ~ dnorm(0, 1/1000)
  delta_b1 ~ dnorm(0, 1/1000)
  phi ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR),  FC=log(df_dates$FrontiersControl), At=df2$At_CMR,
                 tempMa=meteo$temperature_malt)

params <- c("w_b0", "w_b1", "w_b2", "q_b0", "q_b1", "q_b2", "delta_b0", "delta_b1", "phi", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_At_tMa <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar1_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))

save(list="cmed_route_At_tMa", file="Resultats/cmed_route_covs_At_tMa.RData")


### MODEL AMB q_b1 * tempMa + q_b2 * At i delta_b1 * FC
# model igual que l'anterior però valors girats i canvia el resultat

ar1_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0+delta_b1*FC[1])/(1-phi), (sigma2_w/(1-phi^2))^-1)
  for (t in 2:N) {
    x[t] ~ dnorm((delta_b0+delta_b1*FC[t]) + phi*x[t-1], (sigma2_w/(1-phi^2))^-1)
  }
  for (t in 1:N)
  {
    logit(q[t]) <- q_b0+q_b1*tempMa[t]+q_b2*At[t]
    logit(omega[t]) <- w_b0+w_b1*tempMa[t]+w_b2*At[t]
    mu[1, t] <- (delta_b0+delta_b1*FC[t])/(1-phi) 
    mu[2, t] <- q[t]*(delta_b0+delta_b1*FC[t])/(1-phi)
    tau.m[1, t] <- (sigma2_w/(1-phi^2))^-1
    tau.m[2, t] <- (1/(q[t]^2))*(sigma2_w/(1-phi^2))^-1
    y[t] ~ dnormmix(mu[, t], tau.m[, t], p.m[, t])
    p.m[1, t] <- 1-omega[t]
    p.m[2, t] <- omega[t]
    component_chosen[t] ~ dcat(p.m[, t])
  }
  # priors
  q_b0 ~ dnorm(0, 0.01)
  q_b1 ~ dnorm(0, 0.01)
  q_b2 ~ dnorm(0, 0.01)
  w_b0 ~ dnorm(0, 0.01)
  w_b1 ~ dnorm(0, 0.01)
  w_b2 ~ dnorm(0, 0.01)
  delta_b0 ~ dnorm(0, 1/1000)
  delta_b1 ~ dnorm(0, 1/1000)
  phi ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR),  FC=log(df_dates$FrontiersControl), At=df2$At_CMR,
                 tempMa=meteo$temperature_malt)

params <- c("w_b0", "w_b1", "w_b2", "q_b0", "q_b1", "q_b2", "delta_b0", "delta_b1", "phi", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_tMa_At <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar1_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))

save(list="cmed_route_tMa_At", file="Resultats/cmed_route_covs_tMa_At.RData")








