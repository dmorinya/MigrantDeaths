#### Models comparison
library(R2jags)

### Standard models (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_std_ar1.RData") #DIC=536.4834
load("Results/cmed_route_std_ar2.RData") #DIC=538.9838
load("Results/cmed_route_std_ar3.RData") #DIC=531.8396

### Underreported models (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR1.RData") #DIC=466.518
load("Results/cmed_route_covsAR2.RData") #DIC=466.1173
load("Results/cmed_route_covsAR3.RData") #DIC=466.5408

### Standard models (EASTERN MEDITERRANEAN ROUTE)
load("Results/emed_route_std_ar1.RData") #DIC=68.12661
load("Results/emed_route_std_ar2.RData") #DIC=70.37025
load("Results/emed_route_std_ar3.RData") #DIC=72.14745

### Underreported models (EASTERN MEDITERRANEAN ROUTE)
load("Results/emed_route_covsAR1.RData") #DIC= -161.517
load("Results/emed_route_covsAR2.RData") #DIC=-161.1473
load("Results/emed_route_covsAR3.RData") #DIC=-160.48

### Standard models (WESTERN MEDITERRANEAN ROUTE)
load("Results/wmed_route_std_ar1.RData") #DIC=409.6808
load("Results/wmed_route_std_ar2.RData") #DIC=409.0943
load("Results/wmed_route_std_ar3.RData") #DIC=410.5748

### Underreported models (WESTERN MEDITERRANEAN ROUTE)
load("Results/wmed_route_covsAR1.RData") #DIC=263.6892
load("Results/wmed_route_covsAR2.RData") #DIC=263.3959
load("Results/wmed_route_covsAR3.RData") #DIC=263.7694

#### Best models estimates (CENTRAL MEDITERRANEAN ROUTE)
load("Results/cmed_route_covsAR2.RData") 
cmed_route_mcmc <- coda::as.mcmc(cmed_route)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)
median(cmed_route_mcmc_df$delta_b0); quantile(cmed_route_mcmc_df$delta_b0, 0.025); quantile(cmed_route_mcmc_df$delta_b0, 0.975)
median(cmed_route_mcmc_df$phi1); quantile(cmed_route_mcmc_df$phi1, 0.025); quantile(cmed_route_mcmc_df$phi1, 0.975)
median(cmed_route_mcmc_df$phi2); quantile(cmed_route_mcmc_df$phi2, 0.025); quantile(cmed_route_mcmc_df$phi2, 0.975)
median(cmed_route_mcmc_df$q_b0); quantile(cmed_route_mcmc_df$q_b0, 0.025); quantile(cmed_route_mcmc_df$q_b0, 0.975)
median(cmed_route_mcmc_df$q_b1); quantile(cmed_route_mcmc_df$q_b1, 0.025); quantile(cmed_route_mcmc_df$q_b1, 0.975)
median(cmed_route_mcmc_df$w_b0); quantile(cmed_route_mcmc_df$w_b0, 0.025); quantile(cmed_route_mcmc_df$w_b0, 0.975)
median(cmed_route_mcmc_df$w_b1); quantile(cmed_route_mcmc_df$w_b1, 0.025); quantile(cmed_route_mcmc_df$w_b1, 0.975)

#### Best models estimates (EASTERN MEDITERRANEAN ROUTE)
load("Results/emed_route_covsAR1.RData") 
emed_route_mcmc <- coda::as.mcmc(emed_route)
emed_route_mcmc_df <- codatools::coda_df(emed_route_mcmc)
median(emed_route_mcmc_df$delta_b0); quantile(emed_route_mcmc_df$delta_b0, 0.025); quantile(emed_route_mcmc_df$delta_b0, 0.975)
median(emed_route_mcmc_df$phi); quantile(emed_route_mcmc_df$phi, 0.025); quantile(emed_route_mcmc_df$phi, 0.975)
median(emed_route_mcmc_df$q_b0); quantile(emed_route_mcmc_df$q_b0, 0.025); quantile(emed_route_mcmc_df$q_b0, 0.975)
median(emed_route_mcmc_df$q_b1); quantile(emed_route_mcmc_df$q_b1, 0.025); quantile(emed_route_mcmc_df$q_b1, 0.975)
median(emed_route_mcmc_df$w_b0); quantile(emed_route_mcmc_df$w_b0, 0.025); quantile(emed_route_mcmc_df$w_b0, 0.975)
median(emed_route_mcmc_df$w_b1); quantile(emed_route_mcmc_df$w_b1, 0.025); quantile(emed_route_mcmc_df$w_b1, 0.975)

#### Best models estimates (WESTERN MEDITERRANEAN ROUTE)
load("Results/wmed_route_covsAR2.RData") 
wmed_route_mcmc <- coda::as.mcmc(wmed_route)
wmed_route_mcmc_df <- codatools::coda_df(wmed_route_mcmc)
median(wmed_route_mcmc_df$delta_b0); quantile(wmed_route_mcmc_df$delta_b0, 0.025); quantile(wmed_route_mcmc_df$delta_b0, 0.975)
median(wmed_route_mcmc_df$phi1); quantile(wmed_route_mcmc_df$phi1, 0.025); quantile(wmed_route_mcmc_df$phi1, 0.975)
median(wmed_route_mcmc_df$phi2); quantile(wmed_route_mcmc_df$phi2, 0.025); quantile(wmed_route_mcmc_df$phi2, 0.975)
median(wmed_route_mcmc_df$q_b0); quantile(wmed_route_mcmc_df$q_b0, 0.025); quantile(wmed_route_mcmc_df$q_b0, 0.975)
median(wmed_route_mcmc_df$q_b1); quantile(wmed_route_mcmc_df$q_b1, 0.025); quantile(wmed_route_mcmc_df$q_b1, 0.975)
median(wmed_route_mcmc_df$w_b0); quantile(wmed_route_mcmc_df$w_b0, 0.025); quantile(wmed_route_mcmc_df$w_b0, 0.975)
median(wmed_route_mcmc_df$w_b1); quantile(wmed_route_mcmc_df$w_b1, 0.025); quantile(wmed_route_mcmc_df$w_b1, 0.975)
