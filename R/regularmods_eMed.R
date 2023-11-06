library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(astsa)
library(tseries)

options(scipen = 999)

df <- read.csv(file = "Data/df_daily_deaths_wide.csv", sep = ",")
df$date <- as.Date(df$date, format = "%Y-%m-%d")
df <- df[!duplicated(df$date), ]
df2 <- df %>% group_by(month = floor_date(date, "month")) %>%
  summarize(deaths_CMR = sum(deaths_CMR, na.rm = T), arrivals_CMR = mean(arrivals_CMR, na.rm = T),
            deaths_WMR = sum(deaths_WMR, na.rm = T), arrivals_WMR = mean(arrivals_WMR, na.rm = T),
            deaths_EMR = sum(deaths_EMR, na.rm = T), arrivals_EMR = mean(arrivals_EMR, na.rm = T))
missdates <- data.frame(month=c(as.Date("2013-02-01", format="%Y-%m-%d"), as.Date("2013-06-01", format="%Y-%m-%d")), deaths_CMR=c(0, 7),
                        arrivals_CMR=c(232, 3563), deaths_WMR=c(3, 1), arrivals_WMR=c(420, 617), 
                        deaths_EMR=c(2, 19), arrivals_EMR=c(931, 1782))
df3 <- rbind(df2, missdates)
df3 <- df3[order(df3$month), ]
rm(df, df2, missdates)
df3 <- df3[df3$month <= "2021-06-01", ]

### Kalman filter estimation
y<-df3$deaths_EMR
A<-1
mu0<-0
Sigma0<-10000
Phi<-0.1
cQ<-0.1
cR<-0.1

Linn<-function(par){
  Phi<-par[1]
  cQ<-par[2]
  cR<-par[3]
  kf<-Kfilter(y, A, mu0, Sigma0, Phi, cQ, cR)
  kf$like
}
initpar<-(c(0.1,0.1,0.1))
sol<-optim(initpar,Linn,NULL,method="BFGS",hessian=TRUE)
sol$par

### Reconstruccion de Kalman
Phi<-sol$par[1]
cQ<-sol$par[2]
cR<-sol$par[3]

kf<-Ksmooth(y, A, mu0, Sigma0, Phi, t(cQ), t(cR))
dCMRsmo<-ts(c(unlist(kf$Xs)),start=c(2009,1),freq=12)

### AR(1)
ar1 <- arima(dEMRsmo, order=c(1,0,0))
garch  <- garch(ar1$residuals, order=c(0, 1))
AIC_AR1 <- AIC(garch)

### AR(2)
ar2 <- arima(dEMRsmo, order=c(2,0,0))
garch  <- garch(ar2$residuals, order=c(0, 1))
AIC_AR2 <- AIC(garch)

### ARMA(1, 1)
arma11 <- arima(dEMRsmo, order=c(1,0,1))
garch  <- garch(arma11$residuals, order=c(0, 1))
AIC_ARMA11 <- AIC(garch)

### SAR(1)
sar1 <- arima(dEMRsmo, order=c(1,0,0), seasonal=list(order=c(1,0,0), period=12))
garch  <- garch(sar1$residuals, order=c(0, 1))
AIC_SAR1 <- AIC(garch)

### SAR(2)
sar2 <- arima(dEMRsmo, order=c(2,0,0), seasonal=list(order=c(1,0,0), period=12))
garch  <- garch(sar2$residuals, order=c(0, 1))
AIC_SAR2 <- AIC(garch)

### SARMA(1, 1)
sarma11 <- arima(dEMRsmo, order=c(1,0,1), seasonal=list(order=c(1,0,0), period=12))
garch  <- garch(sarma11$residuals, order=c(0, 1))
AIC_SARMA11 <- AIC(garch)


