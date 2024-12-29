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
  delta_b0 ~ dnorm(0, 1/1000)
  phi ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}

### CENTRAL MEDITERRANEAN
ncores <- detectCores() - 1
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores)

data_CMR <- list(x=df2$mortality_rate*100, N=length(df2$mortality_rate))
params <- c("delta_b0", "phi", "sigma2_w")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route <- jags.parallel(data = data_CMR, inits = NULL, parameters.to.save = params,
                                        model.file = ar1_model, n.chains = 5, n.iter = 500000,
                                        n.burnin = 20000, n.thin = 10000, DIC = T, jags.module = "mix"))
stopCluster(cl)
save(list="cmed_route", file="../Results/cmed_route_std_ar1.RData")
