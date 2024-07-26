# Reading libraries

# we use the following libraries
'
library(gsheet)
library(readxl)
library(lubridate)
library(tidyverse)
'

# List of packages to install
packages <- c("gsheet", "readxl", "lubridate","tidyverse")

# Function to install packages if they're not already installed
install_if_not_installed <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  # Load the package
  library(package, character.only = TRUE)
}

# Install and load the packages
sapply(packages, install_if_not_installed)

# Missing migrants
url_mmp <- "https://missingmigrants.iom.int/sites/g/files/tmzbdl601/files/report-migrant-incident/Missing_Migrants_Global_Figures_allData.xlsx"
temp1 <- tempfile(fileext = ".xlsx")
download.file(url_mmp, temp1, mode = "wb")

# Read the Excel file
df_mmp <- read_excel(temp1)

# Transform data and aggregate by day
df_mmp_daily <- df_mmp %>% 
  mutate(date = ymd(`Incident Date`),
         route_name = `Migration Route`) %>% 
  group_by(route_name, date) %>% 
  summarise(total_dead_and_missing = sum(`Total Number of Dead and Missing`, na.rm=T),
            total_survivors = sum(`Number of Survivors`, na.rm=T)) %>% 
  ungroup()

df_mmp_daily_med <- df_mmp_daily %>% 
  filter(route_name %in% c("Central Mediterranean",
                           "Eastern Mediterranean",
                           "Western Mediterranean",
                           "Western Africa / Atlantic route to the Canary Islands")) %>% 
  mutate(route = case_when(
    route_name %in% c("Central Mediterranean") ~ "CMR",
    route_name %in% c("Eastern Mediterranean") ~ "EMR",
    route_name %in% c("Western Mediterranean") ~ "WMR",
    route_name %in% c("Western Africa / Atlantic route to the Canary Islands") ~ "WAR")) %>% 
  select(date, route, total_dead_and_missing, total_survivors)
  

# Migrant files
# Specify the URL of the Google Sheets file
url_migfiles <- "https://docs.google.com/spreadsheets/d/1YNqIzyQfEn4i_be2GGWESnG2Q80E_fLASffsXdCOftI/edit#gid=1085726718"

df_migfiles <- gsheet2tbl(url_migfiles, sheetid = "Events")

# Transform data and aggregate by day
df_migfiles_daily <- df_migfiles %>% 
  filter(!is.na(`route (Frontex)`)) %>% 
  mutate(date_short = substr(date, start=1, stop=10),
         date_recode = ymd(`date_short`),
         route_name = `route (Frontex)`) %>% 
  group_by(route_name, date_recode) %>% 
  summarise(total_dead_and_missing = sum(dead_and_missing, na.rm=T)) %>% 
  ungroup()

df_migfiles_daily_med <- df_migfiles_daily %>% 
  rename(date = date_recode) %>% 
  filter(route_name %in% c("Central Mediterranean route",
                           "Central Mediterranean Route",
                           "Eastern Mediterranean route",
                           "Western African route",
                           "Western Mediterranean route",
                           "Western Mediterranean Route")) %>% 
  mutate(route = case_when(
    route_name %in% c("Central Mediterranean route","Central Mediterranean Route") ~ "CMR",
    route_name %in% c("Eastern Mediterranean route") ~ "EMR",
    route_name %in% c("Western Mediterranean route","Western Mediterranean Route") ~ "WMR",
    route_name %in% c("Western African route") ~ "WAR")
  ) %>% 
  select(date, route, total_dead_and_missing)

# Frontex arrivals

url_frontex <- "https://www.frontex.europa.eu/assets/Migratory_routes/2024/Monthly_detections_of_IBC_2024_01_05.xlsx"
temp3 <- tempfile(fileext = ".xlsx")
download.file(url_frontex, temp3, mode = "wb")

# Read the Excel file
df_frontex <- read_excel(temp3)

# data transformations
df_frontex_monthly_med <- df_frontex %>% 
  pivot_longer(cols = 4:ncol(df_frontex)) %>% 
  mutate(date = ymd(paste0(substr(name, start=4, stop=8),"-",substr(name, start=1, stop=3),"-",01)),
         route = case_when(
           Route %in% c("Black Sea Route", "Circular Route from Albania to Greece","Eastern Borders Route", "Other","Western Balkan Route") ~ "OR",
           Route %in% c("Central Mediterranean Route") ~"CMR",
           Route %in% c("Eastern Mediterranean Route") ~ "EMR",
           Route %in% c("Western African Route") ~ "WAR",
           Route %in% c("Western Mediterranean Route") ~ "WMR",
         )) %>% 
  filter(route %in% c("CMR","WMR","EMR","WAR")) %>% 
  group_by(route, date) %>% 
  summarize(total_arrivals_frontex = sum(value, na.rm=T)) %>% 
  select(date, route, total_arrivals_frontex) %>% 
  ungroup()

ggplot(df_frontex_long, aes(x=date, y=total_arrivals_frontex)) + geom_line() + facet_grid(~route)

# Join data and aggregate to the month

df_daily_med <- bind_rows(df_migfiles_daily_med %>% filter(date < ymd(min(df_mmp_daily_med$date))),
                    df_mmp_daily_med)

df_daily_med <- df_daily_med %>% 
  mutate(date_month = ymd(paste0(year(date),"-",month(date),"-",01))) %>% 
  group_by(date_month, route) %>% 
  mutate(share_deaths_day = total_dead_and_missing/sum(total_dead_and_missing)) # used to compute HHI (fractionalization)

# Here we create the monthly file

df_monthly_med <- df_daily_med %>% 
  mutate(date_month = ymd(paste0(year(date),"-",month(date),"-",01))) %>% 
  group_by(date_month, route) %>% 
  summarise(total_dead_and_missing = sum(total_dead_and_missing, na.rm=T),
            total_survivors = ifelse(date_month < min(df_mmp_daily_med$date),NA,sum(total_survivors, na.rm = T)),
            hhi_deaths_month = sum(share_deaths_day^2)) %>% # HHI (fractionalization index)
  ungroup() %>% 
  rename(date = date_month) %>% 
  select(date,route,total_dead_and_missing,total_survivors,hhi_deaths_month)

# Join the frontex data on arrivals

df_monthly_med <- left_join(df_monthly_med, df_frontex_monthly_med, by=c("date","route"))

g1 <- ggplot(df_monthly_med, aes(x=date, y=total_dead_and_missing)) + geom_col() + facet_wrap(~route, scales="free") + ggtitle("Total death and missing by route")
g2 <- ggplot(df_monthly_med, aes(x=date, y=total_survivors)) + geom_col() + facet_wrap(~route, scales="free") + ggtitle("Total survivors by route")
g3 <- ggplot(df_monthly_med, aes(x=date, y=hhi_deaths_month)) + geom_line() + facet_wrap(~route, scales="free")  + ggtitle("Fractionalization index at the month by route")
g4 <- ggplot(df_monthly_med, aes(x=date, y=total_arrivals_frontex)) + geom_line() + facet_wrap(~route, scales="free")  + ggtitle("Total arrivals by route")

gridExtra::grid.arrange(g1,g2,g3,g4, ncol=4)
