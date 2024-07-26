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
ar2_model <- function()
{
  # likelihood
  x[1] ~ dnorm((delta_b0)/(1-phi1-phi2), (sigma2_w/(1-phi1^2-phi2^2))^-1)
  x[2] ~ dnorm((delta_b0)/(1-phi1-phi2), (sigma2_w/(1-phi1^2-phi2^2))^-1)
  for (t in 3:N) {
    x[t] ~ dnorm((delta_b0) + phi1*x[t-1] + phi2*x[t-2], (sigma2_w/(1-phi1^2-phi2^2))^-1)
  }
  delta_b0 ~ dnorm(0, 1/1000)
  phi1 ~ dunif(-1, 1)
  phi2 ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
data_CMR <- list(x=df2$Incidence_CMR, N=length(df2$Incidence_CMR))

params <- c("delta_b0", "phi1", "phi2", "sigma2_w")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar2_model,
                                                n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="cmed_route", file="Results/cmed_route_std_ar2.RData")

### WESTERN MEDITERRANEAN
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  

df2$Incidence_WMR <- df2$deaths_WMR/(df2$deaths_WMR+df2$survivors_WMR+df2$arrivals_WMR)*100
df2$Incidence_WMR[is.nan(df2$Incidence_WMR)] <- 0.01
data_WMR <- list(x=df2$Incidence_WMR, N=length(df2$Incidence_WMR), At=df2$At_WMR)

params <- c("delta_b0", "phi1", "phi2", "sigma2_w")
clusterExport(cl, list("data_WMR", "params"))
system.time(wmed_route <- jags.parallel(data = data_WMR, inits = NULL, parameters.to.save = params, model.file = ar2_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="wmed_route", file="Results/wmed_route_std_ar2.RData")

### EASTERN MEDITERRANEAN
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  

df2$Incidence_EMR <- df2$deaths_EMR/(df2$deaths_EMR+df2$survivors_EMR+df2$arrivals_EMR)*100
df2$Incidence_EMR[is.nan(df2$Incidence_EMR)] <- 0.01
data_EMR <- list(x=df2$Incidence_EMR, N=length(df2$Incidence_EMR), At=df2$At_EMR)

params <- c("delta_b0", "phi1", "phi2", "sigma2_w")
clusterExport(cl, list("data_EMR", "params"))
system.time(emed_route <- jags.parallel(data = data_EMR, inits = NULL, parameters.to.save = params, model.file = ar2_model,
                                        n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="emed_route", file="Results/emed_route_std_ar2.RData")
