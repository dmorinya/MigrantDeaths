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

df <- read_excel("Data/df_monthly_med.xlsx", guess_max = 10000)
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")

df <- df %>% distinct()

df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals), names_sep = "_")

### Select time period with complete information (2014-2021)
df2 <- df2[df2$date>="2014-01-01" & df2$date < "2021-12-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)

x_rec <- function(par, y, At, z) { #par: q_b0, q_b1
  y <- ifelse(y==0, y+0.01, y)
  x <- vector()
  q <- exp(par[1]+par[2]*At)/(1+exp(par[1]+par[2]*At))
  x <- ifelse(z==1, y, y/q)
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
load("Results/cmed_route_covsAR2.RData")
cmed_route_mcmc <- as.mcmc(cmed_route)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)

result <- list()
for (j in 1:dim(cmed_route_mcmc_df)[1])
{
  z <- cmed_route_mcmc_df[j, grepl("component_chosen", colnames(cmed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  result[[j]] <- x_rec(par=c(cmed_route_mcmc_df$q_b0[j], cmed_route_mcmc_df$q_b1[j]),
                       y=df2$Incidence_CMR, At=df2$At_CMR, z=z) 
}

result_df <- list.rbind(result)
rec_central <- colQuantiles(result_df, probs=0.5)
rec_central25 <- colQuantiles(result_df, probs=0.25)
rec_central75 <- colQuantiles(result_df, probs=0.75)

df1 <- df2[, c("date", "Incidence_CMR")]
colnames(df1) <- c("Date", "Incidence")
df3 <- data.frame(Date=df1$Date, Incidence=rec_central, 
                  Series=rep("Estimated", length(rec_central)))
df1$Series <- "Registered"
df <- rbind(df1, df3)
df$Deaths25 <- c(rep(NA, 95), rec_central25)
df$Deaths75 <- c(rep(NA, 95), rec_central75)
df$Date <- as.Date(df$Date)

CMed_plot <- df %>% 
  ggplot(aes(x=Date, y=Incidence, ymin = Deaths25, ymax = Deaths75, col=Series)) + geom_lineribbon() + ylab("Mortality rate") +
  ylim(0, 125)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), 
        legend.position = "right") + xlab("") + ggtitle("CENTRAL MEDITERRANEAN") 

### Estimated evolution of the number of deaths:
df$Deaths <- (df$Incidence*df2$survivors_CMR+
                        df$Incidence*df2$arrivals_CMR)/(100-df2$Incidence_CMR)

### Cumulated total of deaths:
sum(df$Deaths[df$Series=="Registered"]) ### Registered
round(sum(df$Deaths[df$Series=="Estimated"])) ### Estimated

sum(df2$deaths_CMR)/round(sum(df$Deaths[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total

### EASTERN MEDITERRANEAN
df2$Incidence_EMR <- df2$deaths_EMR/(df2$deaths_EMR+df2$survivors_EMR+df2$arrivals_EMR)*100
df2$Incidence_EMR[is.nan(df2$Incidence_EMR)] <- 0.01
load("Results/emed_route_covsAR1.RData")
emed_route_mcmc <- as.mcmc(emed_route)
emed_route_mcmc_df <- codatools::coda_df(emed_route_mcmc)

result <- list()
for (j in 1:dim(emed_route_mcmc_df)[1])
{
  z <- emed_route_mcmc_df[j, grepl("component_chosen", colnames(emed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  result[[j]] <- x_rec(par=c(emed_route_mcmc_df$q_b0[j], emed_route_mcmc_df$q_b1[j]),
                       y=df2$Incidence_EMR, At=df2$At_EMR, z=z) 
}

result_df <- list.rbind(result)
rec_east <- colQuantiles(result_df, probs=0.5)
rec_east25 <- colQuantiles(result_df, probs=0.25)
rec_east75 <- colQuantiles(result_df, probs=0.75)

df1 <- df2[, c("date", "Incidence_EMR")]
colnames(df1) <- c("Date", "Incidence")
df3 <- data.frame(Date=df1$Date, Incidence=rec_east, 
                  Series=rep("Estimated", length(rec_east)))
df1$Series <- "Registered"
df <- rbind(df1, df3)
df$Deaths25 <- c(rep(NA, 95), rec_east25)
df$Deaths75 <- c(rep(NA, 95), rec_east75)
df$Date <- as.Date(df$Date)

EMed_plot <- df %>% 
  ggplot(aes(x=Date, y=Incidence, ymin = Deaths25, ymax = Deaths75, col=Series)) + geom_lineribbon() + ylab("Mortality rate") +
  ylim(0, 3)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), 
        legend.position = "right") + xlab("") + ggtitle("EASTERN MEDITERRANEAN") 

### Estimated evolution of the number of deaths:
df$Deaths <- (df$Incidence*df2$survivors_EMR+
                df$Incidence*df2$arrivals_EMR)/(100-df2$Incidence_EMR)

### Cumulated total of deaths:
sum(df$Deaths[df$Series=="Registered"]) ### Registered
round(sum(df$Deaths[df$Series=="Estimated"])) ### Estimated

sum(df2$deaths_EMR)/round(sum(df$Deaths[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total

### WESTERN MEDITERRANEAN
df2$Incidence_WMR <- df2$deaths_WMR/(df2$deaths_WMR+df2$survivors_WMR+df2$arrivals_WMR)*100
df2$Incidence_WMR[is.nan(df2$Incidence_WMR)] <- 0.01
load("Results/wmed_route_covsAR2.RData")
wmed_route_mcmc <- as.mcmc(wmed_route)
wmed_route_mcmc_df <- codatools::coda_df(wmed_route_mcmc)

result <- list()
for (j in 1:dim(wmed_route_mcmc_df)[1])
{
  z <- wmed_route_mcmc_df[j, grepl("component_chosen", colnames(wmed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  result[[j]] <- x_rec(par=c(wmed_route_mcmc_df$q_b0[j], wmed_route_mcmc_df$q_b1[j]),
                       y=df2$Incidence_WMR, At=df2$At_WMR, z=z) 
}

result_df <- list.rbind(result)
rec_west <- colQuantiles(result_df, probs=0.5)
rec_west25 <- colQuantiles(result_df, probs=0.25)
rec_west75 <- colQuantiles(result_df, probs=0.75)

df1 <- df2[, c("date", "Incidence_WMR")]
colnames(df1) <- c("Date", "Incidence")
df3 <- data.frame(Date=df1$Date, Incidence=rec_west, 
                  Series=rep("Estimated", length(rec_west)))
df1$Series <- "Registered"
df <- rbind(df1, df3)
df$Deaths25 <- c(rep(NA, 95), rec_west25)
df$Deaths75 <- c(rep(NA, 95), rec_west75)
df$Date <- as.Date(df$Date)

WMed_plot <- df %>% 
  ggplot(aes(x=Date, y=Incidence, ymin = Deaths25, ymax = Deaths75, col=Series)) + geom_lineribbon() + ylab("Mortality rate") +
  ylim(0, 35)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), 
        legend.position = "right") + xlab("") + ggtitle("WESTERN MEDITERRANEAN") 

### Estimated evolution of the number of deaths:
df$Deaths <- (df$Incidence*df2$survivors_WMR+
                df$Incidence*df2$arrivals_WMR)/(100-df2$Incidence_WMR)

### Cumulated total of deaths:
sum(df$Deaths[df$Series=="Registered"]) ### Registered
round(sum(df$Deaths[df$Series=="Estimated"])) ### Estimated

sum(df2$deaths_WMR)/round(sum(df$Deaths[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
