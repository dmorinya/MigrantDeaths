library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(lubridate)
library(R2jags)
library(stringr)
library(rlist)
library(matrixStats)
library(ggdist)
library(doParallel)
library(forecast)
library(knitr)
library(MisRepARMA)
df <- read_excel("Data/df_monthly_med.xlsx", guess_max = 10000)
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")
df <- df %>% distinct()
df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals),
              names_sep = "_")
df2 <- df2[df2$date>="2014-01-01" & df2$date <= "2021-12-01", ]
df2 <- df2 %>% replace(is.na(.), 0)
df2 <- select(df2, -deaths_WAR, - survivors_WAR, - At_WAR, -arrivals_WAR)
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+
                                       df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
df2$Incidence_EMR <- df2$deaths_EMR/(df2$deaths_EMR+df2$survivors_EMR+
                                       df2$arrivals_EMR)*100
df2$Incidence_EMR[is.nan(df2$Incidence_EMR)] <- 0.01
df2$Incidence_WMR <- df2$deaths_WMR/(df2$deaths_WMR+df2$survivors_WMR+
                                       df2$arrivals_WMR)*100
df2$Incidence_WMR[is.nan(df2$Incidence_WMR)] <- 0.01
pacf(df2$Incidence_CMR)
pacf(df2$Incidence_EMR)

pacf(df2$Incidence_EMR)
# MODEL CLÀSSIC AR(1)
model_AR1_CMR <- Arima(df2$Incidence_CMR, order = c(1, 0, 0))
summary(model_AR1_CMR); AIC_AR1_CMR <- AIC(model_AR1_CMR)
BIC_AR1_CMR <- BIC(model_AR1_CMR)
# MODEL CLÀSSIC AR(3)
model_AR3_CMR <- Arima(df2$Incidence_CMR, order = c(3, 0, 0))
summary(model_AR3_CMR); AIC_AR3_CMR <- AIC(model_AR3_CMR)
BIC_AR3_CMR <- BIC(model_AR3_CMR)
model_AR3_EMR <- Arima(df2$Incidence_EMR, order = c(3, 0, 0))
summary(model_AR3_EMR); AIC_AR3_EMR <- AIC(model_AR3_EMR)
BIC_AR3_EMR <- BIC(model_AR3_EMR)
model_AR3_WMR <- Arima(df2$Incidence_WMR, order = c(3, 0, 0))
summary(model_AR3_WMR); AIC_AR3_WMR <- AIC(model_AR3_WMR)
BIC_AR3_WMR <- BIC(model_AR3_WMR)
# MODEL MISREP ARMA
AR1_CMR <- fitMisRepARMA(df2$Incidence_CMR, tol=1e-16, B=1000, p_AR=1, q_MA=0)
summary(AR1_CMR)
AR3_CMR <- fitMisRepARMA(df2$Incidence_CMR, tol=1e-16, B=1000, p_AR=3, q_MA=0)
summary(AR3_CMR)
AR1_EMR <- fitMisRepARMA(df2$Incidence_EMR, tol=1e-16, B=1000, p_AR=1, q_MA=0)
summary(AR1_EMR)
AR3_EMR <- fitMisRepARMA(df2$Incidence_EMR, tol=1e-16, B=1000, p_AR=3, q_MA=0)
summary(AR3_EMR)
AR1_WMR <- fitMisRepARMA(df2$Incidence_WMR, tol=1e-16, B=1000, p_AR=1, q_MA=0)
summary(AR1_WMR)
AR3_WMR <- fitMisRepARMA(df2$Incidence_WMR, tol=1e-16, B=1000, p_AR=3, q_MA=0)
summary(AR3_WMR)
# RECONSTRUCCIÓ
REC1_CMR <- reconstruct(AR1_CMR)
REC3_CMR <- reconstruct(AR3_CMR)
REC1_EMR <- reconstruct(AR1_EMR)
REC3_EMR <- reconstruct(AR3_EMR)
REC1_WMR <- reconstruct(AR1_WMR)
REC3_WMR <- reconstruct(AR3_WMR)

plot(df2$date, df2$Incidence_CMR, type="l",col="blue4", lty=1, ylim = c(0,35),
     xlab = "Any", ylab = "Incidència",main = "Reconstrucció amb AR(3)")
lines(df2$date, REC3_CMR, col="gold1", lty=2, lwd=2)
legend("topleft", legend=c("Sèrie registrada", "Sèrie ajustant subnotificació"),
       col=c("blue4", "gold1"), lty=c(2,3), cex=0.8)
plot(df2$date, df2$Incidence_EMR, type="l",col="blue4", lty=1, ylim = c(0,10),
     xlab = "Any", ylab = "Incidència",main = "Reconstrucció amb AR(3)")
lines(df2$date, REC3_EMR, col="red4", lty=3, lwd=2)
legend("topleft", legend=c("Sèrie registrada", "Sèrie ajustant subnotificació"),
       col=c("blue4", "red4"), lty=c(2,3), cex=0.8)
plot(df2$date, df2$Incidence_WMR, type="l",col="blue4", lty=1, ylim = c(0,20),
     xlab = "Any", ylab = "Incidència",main = "Reconstrucció amb AR(3)")
lines(df2$date, REC3_WMR, col="red4", lty=3, lwd=2)
legend("topleft", legend=c("Sèrie registrada", "Sèrie ajustant subnotificació"),
       col=c("blue4", "red4"), lty=c(2,3), cex=0.8)
# MODEL BAYESIÀ
ar3_model_std <- function()
{
  x[1] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[2] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[3] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  for (t in 4:N) {
    x[t] ~ dnorm((delta_b0) + phi1*x[t-1] + phi2*x[t-2] + phi3*x[t-3],
                 (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  }
  delta_b0 ~ dnorm(0, 1/1000)
  phi1 ~ dunif(-1, 1)
  phi2 ~ dunif(-1, 1)
  phi3 ~ dunif(-1, 1)
  tau ~ dgamma(0.001, 0.001)
  sigma2_w <- 1/tau
}
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+
                                       df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
data_CMR <- list(x=df2$Incidence_CMR, N=length(df2$Incidence_CMR))

params <- c("delta_b0", "phi1", "phi2", "phi3", "sigma2_w")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_std_ar3 <- jags.parallel(data = data_CMR, inits = NULL,
                                                parameters.to.save = params, model.file = ar3_model_std,
                                                n.chains = 5, n.iter = 5000, n.burnin = 500,
                                                n.thin = 10, DIC = T, jags.module = "mix"))
df2$Incidence_EMR <- df2$deaths_EMR/(df2$deaths_EMR+df2$survivors_EMR+
                                       df2$arrivals_EMR)*100
df2$Incidence_EMR[is.nan(df2$Incidence_EMR)] <- 0.01
data_EMR <- list(x=df2$Incidence_EMR, N=length(df2$Incidence_EMR))
params <- c("delta_b0", "phi1", "phi2", "phi3", "sigma2_w")
clusterExport(cl, list("data_EMR", "params"))
system.time(emed_route_std_ar3 <- jags.parallel(data = data_EMR, inits = NULL,
                                                parameters.to.save = params, model.file = ar3_model_std,
                                                n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin =
                                                  10, DIC = T, jags.module = "mix"))
df2$Incidence_WMR <- df2$deaths_WMR/(df2$deaths_WMR+df2$survivors_WMR+
                                       df2$arrivals_WMR)*100
df2$Incidence_WMR[is.nan(df2$Incidence_WMR)] <- 0.01
data_WMR <- list(x=df2$Incidence_WMR, N=length(df2$Incidence_WMR))
params <- c("delta_b0", "phi1", "phi2", "phi3", "sigma2_w")
clusterExport(cl, list("data_WMR", "params"))
system.time(wmed_route_std_ar3 <- jags.parallel(data = data_WMR, inits = NULL,
                                                parameters.to.save = params, model.file = ar3_model_std,
                                                n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin =
                                                  10, DIC = T, jags.module = "mix"))
# MODEL BAYESIÀ SUBNOTIFICAT
ar3_model <- function()
{
  x[1] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[2] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  x[3] ~ dnorm((delta_b0)/(1-phi1-phi2-phi3), (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
  for (t in 4:N) {
    x[t] ~ dnorm((delta_b0) + phi1*x[t-1] + phi2*x[t-2] + phi3*x[t-3],
                 (sigma2_w/(1-phi1^2-phi2^2-phi3^2))^-1)
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
ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  
data_CMR <- list(y=df2$Incidence_CMR, N=length(df2$Incidence_CMR), At=df2$At_CMR)
params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1",
            "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_CMR", "params"))
system.time(cmed_route_ar3 <- jags.parallel(data = data_CMR, inits = NULL,
                                            parameters.to.save = params, model.file = ar3_model,
                                            n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10,
                                            DIC = T, jags.module = "mix"))
stopCluster(cl)

ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  
data_EMR <- list(y=df2$Incidence_EMR, N=length(df2$Incidence_EMR), At=df2$At_EMR)
params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1",
            "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_EMR", "params"))
system.time(emed_route_ar3 <- jags.parallel(data = data_EMR, inits = NULL,
                                            parameters.to.save = params, model.file = ar3_model,
                                            n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10,
                                            DIC = T, jags.module = "mix"))
stopCluster(cl)

ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores)  
data_WMR <- list(y=df2$Incidence_WMR, N=length(df2$Incidence_WMR), At=df2$At_WMR)
params <- c("w_b0", "w_b1", "q_b0", "q_b1", "delta_b0", "phi1",
            "phi2", "phi3", "sigma2_w", "component_chosen")
clusterExport(cl, list("data_WMR", "params"))
system.time(wmed_route_ar3 <- jags.parallel(data = data_WMR, inits = NULL,
                                            parameters.to.save = params, model.file = ar3_model,
                                            n.chains = 5, n.iter = 5000, n.burnin = 500, n.thin = 10,
                                            DIC = T, jags.module = "mix"))
stopCluster(cl)

# RECONSTRUCCIÓ
x_rec <- function(par, y, At, z) {
  y <- ifelse(y == 0, y + 0.01, y)
  x <- vector()
  q <- exp(par[1] + par[2] *
             At) / (1 + exp(par[1] + par[2] * At))
  x <- ifelse(z == 1, y, y / q)
}
cmed_route_mcmc <- as.mcmc(cmed_route_ar3)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)
resultc <- list()
for (j in 1:dim(cmed_route_mcmc_df)[1])
{
  z <- cmed_route_mcmc_df[j, grepl("component_chosen", colnames(cmed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  resultc[[j]] <- x_rec(par=c(cmed_route_mcmc_df$q_b0[j], cmed_rote_mcmc_df$q_b1[j]),
                        y=df2$Incidence_CMR, At=df2$At_CMR, z=z)
}
result_df_CMR <- list.rbind(resultc)
rec_central_CMR <- colQuantiles(result_df_CMR, probs=0.5)
rec_central25_CMR <- colQuantiles(result_df_CMR, probs=0.25)
rec_central75_CMR <- colQuantiles(result_df_CMR, probs=0.75)
df1c <- df2[, c("date", "Incidence_CMR")]
colnames(df1c) <- c("Date", "Incidence")
df3c <- data.frame(Date=df1c$Date, Incidence=rec_central_CMR,
                   Sèries=rep("Sèrie ajustant subnotificació", length(rec_central_CMR)))
df1c$Sèries <- "Sèrie registrada"

dfc <- rbind(df1c, df3c)
dfc$Deaths25_CMR <- c(rep(NA, 95), rec_central25_CMR)
dfc$Deaths75_CMR <- c(rep(NA, 95), rec_central75_CMR)
dfc$Date <- as.Date(dfc$Date)
CMed_plot <- dfc %>%
  ggplot(aes(x = Date, y = Incidence, ymin = Deaths25_CMR, ymax = Deaths75_CMR,
             col = Sèries)) +
  geom_lineribbon(size = 0.6) +
  ylab("") +
  ylim(0, 150) +
  scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.18, 0.84)
  ) +
  xlab("") +
  ggtitle("Reconstrucció ruta Mediterrani Central") +
  scale_color_manual(values = c("gold1", "blue4"))
dfc$Deaths_CMR <- (dfc$Incidence*df2$survivors_CMR+
                     dfc$Incidence*df2$arrivals_CMR)/(100-df2$Incidence_CMR)
emed_route_mcmc <- as.mcmc(emed_route_ar3)
emed_route_mcmc_df <- codatools::coda_df(emed_route_mcmc)
resulte <- list()
for (j in 1:dim(emed_route_mcmc_df)[1])
{
  z <- emed_route_mcmc_df[j, grepl("component_chosen", colnames(emed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  resulte[[j]] <- x_rec(par=c(emed_route_mcmc_df$q_b0[j], emed_route_mcmc_df$q_b1[j]),
                        y=df2$Incidence_EMR, At=df2$At_EMR, z=z)
}
result_df_EMR <- list.rbind(resulte)
rec_central_EMR <- colQuantiles(result_df_EMR, probs=0.5)
rec_central25_EMR <- colQuantiles(result_df_EMR, probs=0.25)
rec_central75_EMR <- colQuantiles(result_df_EMR, probs=0.75)

df1e <- df2[, c("date", "Incidence_EMR")]
colnames(df1e) <- c("Date", "Incidence")
df3e <- data.frame(Date=df1e$Date, Incidence=rec_central_EMR,
                   Sèries=rep("Sèrie ajustant subnotificació", length(rec_central_EMR)))
df1e$Sèries <- "Sèrie registrada"
dfe <- rbind(df1e, df3e)
dfe$Deaths25_EMR <- c(rep(NA, 96), rec_central25_EMR)
dfe$Deaths75_EMR <- c(rep(NA, 96), rec_central75_EMR)
dfe$Date <- as.Date(dfe$Date)
EMed_plot <- dfe %>%
  ggplot(aes(x = Date, y = Incidence, ymin = Deaths25_EMR, ymax = Deaths75_EMR,
             col = Sèries)) +
  geom_lineribbon(size = 0.6) +
  ylab("") +
  ylim(0, 2) +
  scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.18, 0.84)
  ) +
  xlab("") +
  ggtitle("Reconstrucció ruta Mediterrani Oriental") +
  scale_color_manual(values = c("gold1", "blue4"))
dfe$Deaths_EMR <- (dfe$Incidence*df2$survivors_EMR+
                     dfe$Incidence*df2$arrivals_EMR)/(100-df2$Incidence_EMR)
wmed_route_mcmc <- as.mcmc(wmed_route_ar3)
wmed_route_mcmc_df <- codatools::coda_df(wmed_route_mcmc)
resultw <- list()
for (j in 1:dim(wmed_route_mcmc_df)[1])
{
  z <- wmed_route_mcmc_df[j, grepl("component_chosen", colnames(wmed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  resultw[[j]] <- x_rec(par=c(wmed_route_mcmc_df$q_b0[j], wmed_route_mcmc_df$q_b1[j]),
                        y=df2$Incidence_WMR, At=df2$At_WMR, z=z)
}

result_df_WMR <- list.rbind(resultw)
rec_central_WMR <- colQuantiles(result_df_WMR, probs=0.5)
rec_central25_WMR <- colQuantiles(result_df_WMR, probs=0.25)
rec_central75_WMR <- colQuantiles(result_df_WMR, probs=0.75)
df1w <- df2[, c("date", "Incidence_WMR")]
colnames(df1w) <- c("Date", "Incidence")
df3w <- data.frame(Date=df1w$Date, Incidence=rec_central_WMR,
                   Sèries=rep("Sèrie ajustant subnotificació", length(rec_central_WMR)))
df1w$Sèries <- "Sèrie registrada"
dfw <- rbind(df1w, df3w)
dfw$Deaths25_WMR <- c(rep(NA, 95), rec_central25_WMR)
dfw$Deaths75_WMR <- c(rep(NA, 95), rec_central75_WMR)
dfw$Date <- as.Date(dfw$Date)
WMed_plot <- dfw %>%
  ggplot(aes(x = Date, y = Incidence, ymin = Deaths25_WMR, ymax = Deaths25_WMR,
             col = Sèries)) +
  geom_lineribbon(size = 0.6) +
  ylab("") +
  ylim(0, 15) +
  scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.18, 0.84)
  ) +
  xlab("") +
  ggtitle("Reconstrucció ruta Mediterrani Occidental") +
  scale_color_manual(values = c("gold1", "blue4"))
dfw$Deaths_WMR <- (dfw$Incidence*df2$survivors_WMR+
                     dfw$Incidence*df2$arrivals_WMR)/(100-df2$Incidence_WMR)

