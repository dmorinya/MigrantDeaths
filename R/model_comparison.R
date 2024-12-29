#### Models comparison
library(R2jags)

### Standard models (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_std_ar0.RData") #DIC=478.1066
load("Resultsx100/cmed_route_std_ar1.RData") #DIC=477.9009
load("Resultsx100/cmed_route_std_ar2.RData") #DIC=482.2191
load("Resultsx100/cmed_route_std_ar3.RData") #DIC=473.0417

### Underreported full models (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_AR0.RData") #DIC=413.6247
load("Resultsx100/cmed_route_AR1.RData") #DIC=413.0294
load("Resultsx100/cmed_route_AR2.RData") #DIC=412.6529
load("Resultsx100/cmed_route_AR3.RData") #DIC=413.7146

### Underreported model 2 (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_AR0_2.RData") #DIC=412.591
load("Resultsx100/cmed_route_AR1_2.RData") #DIC=412.558
load("Resultsx100/cmed_route_AR2_2.RData") #DIC=413.0414
load("Resultsx100/cmed_route_AR3_2.RData") #DIC=413.462

### Underreported model 3 (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_AR0_3.RData") #DIC=413.3053
load("Resultsx100/cmed_route_AR1_3.RData") #DIC=413.063
load("Resultsx100/cmed_route_AR2_3.RData") #DIC=413.3024
load("Resultsx100/cmed_route_AR3_3.RData") #DIC=413.4159

### Underreported model 4 (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_AR0_4.RData") #DIC=414.102
load("Resultsx100/cmed_route_AR1_4.RData") #DIC=412.9352
load("Resultsx100/cmed_route_AR2_4.RData") #DIC=413.0852
load("Resultsx100/cmed_route_AR3_4.RData") #DIC=413.7725

#### Best models estimates (CENTRAL MEDITERRANEAN ROUTE)
load("Resultsx100/cmed_route_AR1_2.RData")
cmed_route_mcmc <- coda::as.mcmc(cmed_route)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)
median(cmed_route_mcmc_df$delta_b0); quantile(cmed_route_mcmc_df$delta_b0, 0.025); quantile(cmed_route_mcmc_df$delta_b0, 0.975)
median(cmed_route_mcmc_df$phi); quantile(cmed_route_mcmc_df$phi, 0.025); quantile(cmed_route_mcmc_df$phi, 0.975)
median(cmed_route_mcmc_df$q_b0); quantile(cmed_route_mcmc_df$q_b0, 0.025); quantile(cmed_route_mcmc_df$q_b0, 0.975)
median(cmed_route_mcmc_df$q_b1); quantile(cmed_route_mcmc_df$q_b1, 0.025); quantile(cmed_route_mcmc_df$q_b1, 0.975)
median(cmed_route_mcmc_df$w_b0); quantile(cmed_route_mcmc_df$w_b0, 0.025); quantile(cmed_route_mcmc_df$w_b0, 0.975)
