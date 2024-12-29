library(dplyr)
library(R2jags)
library(doParallel)

load("../Data/df_update_12_2024.Rdata")

###
### Select time period with complete information (depending on the mortality rate definition)
df2 <- df[df$date_month >="2016-01-01" & df$date_month  <= "2023-11-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)

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
ncores <- detectCores() - 1
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores)

data_CMR <- list(y=df2$mortality_rate*100, N=length(df2$mortality_rate),
                 At=df2$hhi_deaths_month)

params <- c("w_b0", "q_b0", "q_b1", "w_b1", "phi", "delta_b0", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params,
                                        model.file = ar1_model, n.chains = 5, n.iter = 500000,
                                        n.burnin = 20000, n.thin = 10000, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="cmed_route", file="../Results/cmed_route_AR1.RData")
