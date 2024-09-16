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

source("R/get_data.R")
df <- df_monthly_med %>% distinct()
names(df) <- c("date", "route", "deaths", "survivors", "At", "arrivals")

df2 <- df %>%
  pivot_wider(names_from = route, values_from = c(deaths, survivors, At, arrivals), names_sep = "_")

### Select time period with complete information (2014-2023)
df2 <- df2[df2$date>="2014-01-01" & df2$date <= "2023-12-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)

x_rec <- function(par, y, z) { #par: q_b0
  y <- ifelse(y==0, y+0.01, y)
  x <- vector()
  q <- exp(par[1])/(1+exp(par[1]))
  x <- ifelse(z==1, y, y/q)
}

### CENTRAL MEDITERRANEAN
df2$Incidence_CMR <- df2$deaths_CMR/(df2$deaths_CMR+df2$survivors_CMR+df2$arrivals_CMR)*100
df2$Incidence_CMR[is.nan(df2$Incidence_CMR)] <- 0.01
load("Results/cmed_route_covsAR3_3.RData")
cmed_route_mcmc <- as.mcmc(cmed_route)
cmed_route_mcmc_df <- codatools::coda_df(cmed_route_mcmc)

result <- list()
for (j in 1:dim(cmed_route_mcmc_df)[1])
{
  z <- cmed_route_mcmc_df[j, grepl("component_chosen", colnames(cmed_route_mcmc_df))==TRUE]
  ord <- str_sort(colnames(z), numeric = TRUE)
  z <- z[, ord]
  z <- as.numeric(z)
  result[[j]] <- x_rec(par=c(cmed_route_mcmc_df$q_b0[j]),
                       y=df2$Incidence_CMR, z=z) 
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
df$inc25 <- c(rep(NA, 120), rec_central25)
df$inc75 <- c(rep(NA, 120), rec_central75)
df$Date <- as.Date(df$Date)

CMed_plot <- df %>% 
  ggplot(aes(x=Date, y=Incidence, ymin = inc25, ymax = inc75, col=Series)) + geom_lineribbon() + ylab("Mortality rate") +
  ylim(0, 170)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), 
        legend.position = "right") + xlab("") + ggtitle("") 

### Estimated evolution of the number of deaths:
df$Deaths <- (df$Incidence*df2$survivors_CMR+
                        df$Incidence*df2$arrivals_CMR)/(100-df2$Incidence_CMR)
df$Deaths25 <- (df$inc25*df2$survivors_CMR+
                  df$inc25*df2$arrivals_CMR)/(100-df2$Incidence_CMR)
df$Deaths75 <- (df$inc75*df2$survivors_CMR+
                  df$inc75*df2$arrivals_CMR)/(100-df2$Incidence_CMR)

### Cumulated total of deaths:
sum(df$Deaths[df$Series=="Registered"]) ### Registered
round(sum(df$Deaths[df$Series=="Estimated"])) ### Estimated
round(sum(df$Deaths25[df$Series=="Estimated"])) ### Estimated (P25)
round(sum(df$Deaths75[df$Series=="Estimated"])) ### Estimated (P75)

sum(df2$deaths_CMR)/round(sum(df$Deaths[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
sum(df2$deaths_CMR)/round(sum(df$Deaths25[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
sum(df2$deaths_CMR)/round(sum(df$Deaths75[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
