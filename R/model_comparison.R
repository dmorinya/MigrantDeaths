#### Models comparison
library(R2jags)

### Standard models (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_std_ar1.RData") #DIC=708.1981
load("Results/cmed_route_std_ar2.RData") #DIC=709.3282
load("Results/cmed_route_std_ar3.RData") #DIC=705.859

### Underreported full models (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR1.RData") #DIC=595.5421
load("Results/cmed_route_covsAR2.RData") #DIC=591.7625
load("Results/cmed_route_covsAR3.RData") #DIC=594.2973

### Underreported model 2 (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR1_2.RData") #DIC=594.8627
load("Results/cmed_route_covsAR2_2.RData") #DIC=595.6263
load("Results/cmed_route_covsAR3_2.RData") #DIC=594.5411

### Underreported model 3 (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR1_3.RData") #DIC=592.0279
load("Results/cmed_route_covsAR2_3.RData") #DIC=592.2409
load("Results/cmed_route_covsAR3_3.RData") #DIC=591.5445

### Underreported model 4 (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR1_4.RData") #DIC=593.1802
load("Results/cmed_route_covsAR2_4.RData") #DIC=593.189
load("Results/cmed_route_covsAR3_4.RData") #DIC=593.4032

#### Best models estimates (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR3_3.RData") 
cmed_route_mcmc <- coda::as.mcmc(cmed_route)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)
median(cmed_route_mcmc_df$delta_b0); quantile(cmed_route_mcmc_df$delta_b0, 0.025); quantile(cmed_route_mcmc_df$delta_b0, 0.975)
median(cmed_route_mcmc_df$phi1); quantile(cmed_route_mcmc_df$phi1, 0.025); quantile(cmed_route_mcmc_df$phi1, 0.975)
median(cmed_route_mcmc_df$phi2); quantile(cmed_route_mcmc_df$phi2, 0.025); quantile(cmed_route_mcmc_df$phi2, 0.975)
median(cmed_route_mcmc_df$phi3); quantile(cmed_route_mcmc_df$phi3, 0.025); quantile(cmed_route_mcmc_df$phi3, 0.975)
median(cmed_route_mcmc_df$q_b0); quantile(cmed_route_mcmc_df$q_b0, 0.025); quantile(cmed_route_mcmc_df$q_b0, 0.975)
median(cmed_route_mcmc_df$w_b0); quantile(cmed_route_mcmc_df$w_b0, 0.025); quantile(cmed_route_mcmc_df$w_b0, 0.975)
median(cmed_route_mcmc_df$w_b1); quantile(cmed_route_mcmc_df$w_b1, 0.025); quantile(cmed_route_mcmc_df$w_b1, 0.975)