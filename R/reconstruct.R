library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(lubridate)
library(R2jags)
library(doParallel)
library(stringr)
library(rlist)
library(matrixStats)
library(ggdist)

load("Data/df_update_12_2024.Rdata")

###
### Select time period with complete information (depending on the mortality rate definition)
df2 <- df[df$date_month >="2016-01-01" & df$date_month  <= "2023-11-01", ]

### Replace NAs by zeros
df2 <- df2 %>% replace(is.na(.), 0)

x_rec <- function(par, y, z, Ft) { #par: q_b0, q_b1
  y <- ifelse(y==0, y+0.01, y)
  x <- vector()
  q <- exp(par[1]+par[2]*Ft)/(1+exp(par[1]+par[2]*Ft))
  x <- ifelse(z==1, y, y/q)
}

### CENTRAL MEDITERRANEAN
load("Resultsx100/cmed_route_AR1_2.RData")
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
                       y=df2$mortality_rate*100, z=z, Ft=df2$hhi_deaths_month) 
}

result_df <- list.rbind(result)
rec_central <- colQuantiles(result_df, probs=0.5)
rec_central25 <- colQuantiles(result_df, probs=0.25)
rec_central75 <- colQuantiles(result_df, probs=0.75)

df1 <- df2[, c("date_month", "mortality_rate")]
colnames(df1) <- c("Date", "Incidence")
df1$Incidence <- df1$Incidence*100
df3 <- data.frame(Date=df1$Date, Incidence=rec_central, 
                  Series=rep("Estimated", length(rec_central)))
df1$Series <- "Registered"
df <- rbind(df1, df3)
df$inc25 <- c(rep(NA, 95), rec_central25)
df$inc75 <- c(rep(NA, 95), rec_central75)
df$Date <- as.Date(df$Date)

CMed_plot <- df %>% 
  ggplot(aes(x=Date, y=Incidence, ymin = inc25, ymax = inc75, col=Series)) + geom_lineribbon() + ylab("Mortality rate") +
  ylim(0, 100)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), 
        legend.position = "right") + xlab("") + ggtitle("") 

### Estimated evolution of the number of deaths:
df$Deaths <- (df$Incidence*df2$pushbacks+
                        df$Incidence*df2$total_arrivals_frontex)/(100-df$Incidence)
df$Deaths25 <- (df$inc25*df2$pushbacks+
                  df$inc25*df2$total_arrivals_frontex)/(100-df$Incidence)
df$Deaths75 <- (df$inc75*df2$pushbacks+
                  df$inc75*df2$total_arrivals_frontex)/(100-df$Incidence)

### Cumulated total of deaths:
sum(df$Deaths[df$Series=="Registered"]) ### Registered
round(sum(df$Deaths[df$Series=="Estimated"])) ### Estimated
round(sum(df$Deaths25[df$Series=="Estimated"])) ### Estimated (P25)
round(sum(df$Deaths75[df$Series=="Estimated"])) ### Estimated (P75)

sum(df2$total_dead_and_missing)/round(sum(df$Deaths[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
sum(df2$total_dead_and_missing)/round(sum(df$Deaths25[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
sum(df2$total_dead_and_missing)/round(sum(df$Deaths75[df$Series=="Estimated"]))*100 ### Percentage of registered deaths over the estimated total
