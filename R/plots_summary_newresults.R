library(tidyverse)

# reading data
cmr <- read.csv("CMR_incidence_arma11.csv")
wmr <- read.csv("WMR_incidence_arma200.csv")
emr <- read.csv("EMR_incidence.csv")
df <- read.csv("df_daily_deaths_wide.csv")

# summarize to the month
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

df3_long <- df3 %>% 
  select(month, starts_with("death")) %>% 
  pivot_longer(2:4,
               names_prefix  = "deaths_",
               names_to = "route",
               values_to = "observed")

# computing the percentiles
cmr_median <- apply(cmr, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.25,0.5,0.975)))
wmr_median <- apply(wmr, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.25,0.5,0.975)))
emr_median <- apply(emr, MARGIN=2, FUN=function(x) quantile(x, probs=c(0.25,0.5,0.975)))

# transform to data frame
cmr_median_t <- as_tibble(t(cmr_median))
colnames(cmr_median_t) <- c("y_min", "y_median", "y_max")
cmr_median_t$route <- "CMR"
wmr_median_t <- as_tibble(t(wmr_median))
colnames(wmr_median_t) <- c("y_min", "y_median", "y_max")
wmr_median_t$route <- "WMR"
emr_median_t <- as_tibble(t(emr_median))
colnames(emr_median_t) <- c("y_min", "y_median", "y_max")
emr_median_t$route <- "EMR"

# bind these data frames
recstr <- bind_rows(cmr_median_t,wmr_median_t,emr_median_t)
recstr$date <- rep(df3$month, 3) # add the date variable

# join reconstructed with observed values
recstr_n <- left_join(recstr,df3_long %>% rename(date=month), by=c("date","route"))

# compute some stats
recstr_n <- recstr_n %>% 
  mutate(diff = y_median-observed,
         diff_min = y_min-observed,
         diff_max = y_max-observed)

# estimate of total deaths by route during the period
recstr_n %>% 
  group_by(route) %>% 
  summarise(sum(observed),
            sum(y_median))

# Plots
p1 <- ggplot(recstr_n, aes(x=date, y=y_median)) + 
  geom_line() + 
  geom_line(aes(x=date, y=observed, col="red"), show.legend = F) +
  ylab("Number of migrant deaths") +
  xlab("Date") +
  ggtitle("A.") +
  geom_ribbon(aes(ymin=y_min, ymax=y_max), linetype=2, alpha=0.1) +
  theme_bw() + 
  facet_wrap(~ route, scales="free", nrow=3)

p2 <- ggplot(recstr_n, aes(x=date, y=diff)) + 
  geom_line() + 
  geom_hline(yintercept = 0, lty=2, col="red") + 
  geom_ribbon(aes(ymin=diff_min, ymax=diff_max), linetype=2, alpha=0.1) +
  ylab("Difference between reconstructed and observed") +
  xlab("Date")+
  ggtitle("B.") +
  theme_bw() + 
  facet_wrap(~ route, scales="free", nrow=3)

gridExtra::grid.arrange(p1, p2, 
                        nrow=1, 
                        top="Reconstructed and observed number of migrant deaths in the Mediterranean")

