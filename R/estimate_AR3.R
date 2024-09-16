library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(lubridate)
library(R2jags)
library(doParallel)

source("R/get_data.R")
df <- df_monthly_med %>% distinct()
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")

df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals), names_sep = "_")

### Select time period with complete information (2014-2023)
df2 <- df2[df2$date>="2014-01-01" & df2$date <= "2023-12-01", ]

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
data_CMR <- list(x=df2$Incidence_CMR, N=length(df2$Incidence_CMR))

params <- c("delta_b0", "phi1", "phi2", "phi3", "sigma2_w")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params, model.file = ar3_model,
                                                n.chains = 5, n.iter = 50000, n.burnin = 10000, n.thin = 1000, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="cmed_route", file="Results/cmed_route_std_ar3.RData")