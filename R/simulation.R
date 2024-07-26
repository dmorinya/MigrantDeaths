library(R2jags)
library(dplyr)
library(tidyr)
library(doParallel)

simulate_data <- function(N, q_b0, q_b1, w_b0, w_b1, delta_b0, phi, sigma2_w, At) 
{
  q <- exp(q_b0+q_b1*At)/(1+exp(q_b0+q_b1*At))
  w <- exp(w_b0+w_b1*At)/(1+exp(w_b0+w_b1*At))
  x <- vector() # Unobserved process
  z <- rbinom(N, 1, w) # Underreporting indicator
  x[1] <- delta_b0/(1-phi)
  for (i in 2:N) {
    x[i] <- delta_b0 + phi*x[i-1] + rnorm(1, 0, sqrt(sigma2_w))
  }
  y <- ifelse(z == 1, q*x, x) # Observed process 
  return(y)
}

# JAGS to estimate parameters
jags_GBV_sim <- function()
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

###################################################
############ Parameter variation ##############
###################################################

# Define a list of parameter values
gs <- list(q_b0 = c(-1, 0, 1),
           q_b1 = c(-1, 0, 1),
           w_b0 = c(-1, 0, 1),
           w_b1 = c(-1, 0, 1),
           delta_b0 = c(50),
           phi = c(-0.5, 0.5),
           sigma2_w = c(10)) %>% 
  expand.grid()

# Fractionalization index
At <- runif(89)
  
simulation <- function(i)
{  
  print(paste0("simulation_", i))
  for (j in 1:100)
  {
    print(paste0("repetition_", j))
    y <- simulate_data(89, gs$q_b0[i], gs$q_b1[i], gs$w_b0[i], gs$w_b1[i], gs$delta_b0[i], gs$phi[i], gs$sigma2_w[i], At) 
    datos <- list(y=y, N=length(y), At=At)
    params <- c("q_b0", "q_b1", "w_b0", "w_b1", "delta_b0", "phi", "sigma2_w")
    init_values <- function(){
      list(q_b0=rnorm(1), q_b1=rnorm(1), w_b0=rnorm(1), w_b1=rnorm(1), delta_b0=rnorm(1, 10), phi=runif(1, -1, 1), 
           tau=0.1, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
    }

    fit_GBV_sim <- jags.parallel(data = datos, inits = init_values, parameters.to.save = params, 
                                 model.file = jags_GBV_sim, n.chains = 5, n.iter = 5000, 
                                 n.burnin = 1000, n.thin = 10, DIC = T, jags.module = "mix")
    lm1_mcmc <- as.data.frame(summary(as.mcmc(fit_GBV_sim))$quantiles)
    
    #Add real parameters
    lm1_mcmc$real_value <- NULL
    lm1_mcmc[1, "real_value"] <- gs[i, 5]
    lm1_mcmc[3, "real_value"] <- gs[i, 6]
    lm1_mcmc[4, "real_value"] <- gs[i, 1]
    lm1_mcmc[5, "real_value"] <- gs[i, 2]
    lm1_mcmc[6, "real_value"] <- gs[i, 7]
    lm1_mcmc[7, "real_value"] <- gs[i, 3]
    lm1_mcmc[8, "real_value"] <- gs[i, 4]
    
    sim_name <- paste0("sim_n_", i)
    save(lm1_mcmc, file = paste0("Results/Simulation/", sim_name, "_", j, ".RData"))
  }
}

system.time(foreach(i=1:81, combine=rbind, .packages=c("R2jags", "doParallel")) %do% simulation(i))
