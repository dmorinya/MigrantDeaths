library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(astsa)
library(BSL)
library(tseries)
library(quantspec)
library(doParallel)
library(mixtools)

options(scipen = 999)

source("R/Estep.R")
source("R/Mstep.R")
source("R/EM.R")


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
dEMRsmo<-ts(c(unlist(kf$Xs)),start=c(2009,1),freq=12)
xlim1<-c(2009, 2021)

ncores <- detectCores()
cl <- makeCluster(ncores-2)
registerDoParallel(cl)

logPrior <- function(theta)
{
  log(theta[1] > 0 & theta[2] > 0 & theta[2] < 1 &
        abs(theta[4]) < 1 & theta[5] > 0 & theta[6] > 0 & theta[6] < 2 & theta[8] > 0 &
        theta[9] > 0 & theta[9] < 1 & theta[7] > 0 & theta[7] < 1)
}

sim <- function(theta, T)
{
  library(TSA)
  library(quantspec)
  error <- vector()
  error[1] <- 0
  error <- ARCH1(T, theta[5], theta[6])
  x <- theta[1]+sarima.sim(ar=theta[2], ma=theta[3], sar=theta[4], S=12,
                           rand.gen=function(n, ...) rnorm(n, mean=0, sd=theta[8]), n=T)+error
  z  <- rbinom(T, 1, theta[7])
  q <- theta[9]
  y  <- x*(1-z)+q*z*x
  y[y<0] <- 0
  return(y)
}

st <- function(z){
  s1=mean(z); s2=sd(z); s3=acf(z,plot=F)$acf[2]
  s4=acf(z,plot=F)$acf[3]; s5=acf(z,plot=F)$acf[4]
  c(s1,s2,s3,s4,s5)}

addTaskCallback(function(...) {set.seed(123);TRUE})
init <- normalmixEM(dEMRsmo, maxit=10000)
q_init <- init$mu[init$mu==min(init$mu)]/init$mu[init$mu==max(init$mu)]
init_sigma <- init$sigma[init$sigma==max(init$sigma)]
ind <- ifelse(init$posterior[, init$mu==min(init$mu)]<0.5, 0, 1)
x <- ifelse(ind==1, dEMRsmo/q_init, dEMRsmo)
pr <- tryCatch(arima(x, order=c(1,0,1), seasonal=list(order=c(1,0,0), period=12)))
init_ar1  <- tryCatch((garch(!is.na(pr$resid), order=c(0, 1))$coef))

init_vals <- c(mean(x), pr$coef[1], pr$coef[2], pr$coef[3], init_ar1, init$lambda[init$mu==min(init$mu)],
               init_sigma, q_init)

library(quantspec)
model <- newModel(fnSim = sim, fnSum = st,
                  simArgs = list(T = length(dEMRsmo)), theta0 = init_vals,
                  fnLogPrior = logPrior, thetaNames=c(expression(phi[0]), expression(alpha[1]), expression(alpha[2]), "intercept",
                                                      expression(beta[1]), expression(gamma[1]), expression(omega),  expression(sigma),
                                                      "q"))

### Glasso penalty selection
ssy <- st(dEMRsmo)
lambda_all <- list(exp(seq(-7, 0.5, length.out = 20)))
set.seed(1234)
sp_bsl_glasso_eMed <- selectPenalty(ssy = ssy, n = 200,
                                    lambda_all, theta = init_vals, M = 100, sigma = 0.5, model = model,
                                    method = 'BSL', shrinkage = 'glasso')
sp_bsl_glasso_eMed ### Best penalty parameter: lambda = 0.0214

### Estimating the SL
set.seed(1234)
resultMigrant <- bsl(y = as.numeric(dEMRsmo), n = 200, M = 50000, model = model, 
                     diag(c(.0005^2,.0005^2,.0005^2,.0005^2,.0005^2,.0005^2,.0005^2,.0005^2,.0005^2)),
                     method = 'BSL', shrinkage="glasso", penalty=0.0214, parallel=FALSE, verbose=TRUE)

### Parameter estimates and 95% credible intervals
print(show(resultMigrant))

### Pseudo Akaike's Information Criterion
2*dim(resultMigrant@theta)[2]-2*median(resultMigrant@loglike)

### Parameter estimates distribution
est_distr <- "../Results/estad_distr_EMR_sarma101_100.pdf"
pdf(est_distr, width=8.5, height=6.5)
plot(resultMigrant, which = 2, thin = 30,
     options.density = list(color = 'blue4', fill = 'blue', alpha = 0.5),
     options.theme = list(panel.background = element_rect(fill = 'lightgrey'),
                          plot.margin = grid::unit(rep(0.05, 4), "npc")))
dev.off()

posterior_probs <- lapply(resultMigrant@theta[, 9], FUN=function(x){EM(dEMRsmo, 2, 1e-16, x)})
ur_indicator    <- lapply(posterior_probs, FUN=function(x){ifelse(x$gamma[, 1] > x$gamma[, 2], 1, 0)})
q <- resultMigrant@theta[, 9]
rec_values <- list()
for (j in 1:length(ur_indicator))
{
  den <- q[j]
  rec_values[[j]] <- ifelse(ur_indicator[[j]]==0, dEMRsmo, 
                            dEMRsmo/den)
}

pr <- as.data.frame(do.call(rbind, rec_values))

### Export the results
name_data <- "../Results/EMR_incidence_sarma101_100.csv"
write.csv(pr, name_data, row.names=FALSE)

plot(as.numeric(df3$deaths_EMR), type="l", ylim=c(0, 1200))
lines(as.numeric(colMeans(pr)), col="red")

sum(dEMRsmo)/sum(colMeans(pr))*100
sum(df3$deaths_EMR)/sum(colMeans(pr))*100
save(list=ls(), file="../Results/eMed_results_sarma101_100.RData")
