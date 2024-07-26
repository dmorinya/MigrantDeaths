library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(lubridate)
library(R2jags)
library(doParallel)

df <- read_excel("Data/df_monthly_med.xlsx", guess_max = 10000)
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")

df <- df %>% distinct()

df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals), names_sep = "_")

### Select time period with complete information (2014-2021)
df2 <- df2[df2$date>="2014-01-01" & df2$date < "2021-12-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)

# JAGS to estimate parameters
ar3_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[2] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[3] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  for (t in 4:N) {
    x[t] ~ dnorm((delta_b0) + phi1*x[t-1] + phi2*x[t-2] + phi3*x[t-3], (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  }
  for (t in 1:N)
  {
    logit(q[t]) <- q_b0+q_b1*At[t]
    logit(omega[t]) <- w_b0+w_b1*At[t]
    mu[1, t] <- (delta_b0)/(1-phi1-phi2-phi3)
    mu[2, t] <- q[t]*(delta_b0)/(1-phi1-phi2-phi3)
    tau.m[1, t] <- (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1
    tau.m[2, t] <- (1/(q[t]^2))*(sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1
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
  phi1 ~ dunif(-1, 1)
  phi2 ~ dunif(-1, 1)
  phi3 ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  

df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR), At=df2$At_CMR)

params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1", "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar3_model, n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="cmed_route", file="Results/cmed_route_covsAR3.RData")

### WESTERN MEDITERRANEAN
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  

df2$Incidence_WMR <- df2$deaths_WMR/(df2$deaths_WMR+df2$survivors_WMR+df2$arrivals_WMR)*100
df2$Incidence_WMR[is.nan(df2$Incidence_WMR)] <- 0.01
data_WMR <- list(y=df2$Incidence_WMR, N=length(df2$Incidence_WMR), At=df2$At_WMR)

params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1", "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_WMR", "params"))
system.time(wmed_route <- jags.parallel(data = data_WMR, inits = NULL, parameters.to.save = params, model.file = ar3_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="wmed_route", file="Results/wmed_route_covsAR3.RData")

### EASTERN MEDITERRANEAN
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  

df2$Incidence_EMR <- df2$deaths_EMR/(df2$deaths_EMR+df2$survivors_EMR+df2$arrivals_EMR)*100
df2$Incidence_EMR[is.nan(df2$Incidence_EMR)] <- 0.01
data_EMR <- list(y=df2$Incidence_EMR, N=length(df2$Incidence_EMR), At=df2$At_EMR)

params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1", "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_EMR", "params"))
system.time(emed_route <- jags.parallel(data = data_EMR, inits = NULL, parameters.to.save = params, model.file = ar3_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="emed_route", file="Results/emed_route_covsAR3.RData")
