library(matrixStats)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(MMWRweek)

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

### Best fitting model is AR(1)
EMR <- read.table("../Results/EMR_incidence.csv", header=T, sep=",")
q2.5 <- colQuantiles(as.matrix(EMR), probs=c(0.025))
med <- colMedians(as.matrix(EMR))
q97.5 <- colQuantiles(as.matrix(EMR), probs=c(0.975))
resum <- data.frame(q2.5=c(q2.5, rep(NA, length(q2.5))), med=c(med, rep(NA, length(med))),
                    q97.5=c(q97.5, rep(NA, length(q97.5))))
resum$Date <- df3$month
resum$Value <- "Estimated"
resum$med[151:300] <- df3$deaths_EMR
resum$Value[151:300] <- "Registered"
resum$q2.5[151:300] <- resum$med[151:300]
resum$q97.5[151:300] <- resum$med[151:300]
resum$med2 <- c(resum$med[1:150], rep(NA, 150))
round(sum(resum$med[resum$Value=="Registered"])/sum(resum$med[resum$Value=="Estimated"])*100)

png("../Results/EMR_reconstructed.png", width=1200, height=800)
ggplot(data=resum, aes(x=Date, y=med, col=Value)) +
  geom_line()+xlab("Date")+ylab("Deaths")+geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "lightgrey", show.legend=F)+
  labs(col = "Value")+ggtitle("Eastern Mediterranean Route")+
  geom_line(data=resum[!is.na(resum$med2),], aes(x=Date, y=med2), color = "red", linetype = "dotted")+
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"), axis.text.x = element_text(hjust = 1), axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  labs(col = "Value")
dev.off()