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

df_mmp_monthly_med <- df_mmp_daily_med %>% 
  mutate(date_month = ymd(paste0(year(date),"-",month(date),"-",01))) %>% 
  group_by(date_month, route) %>% 
  summarise(total_dead_and_missing = sum(total_dead_and_missing, na.rm=T),
            total_survivors = sum(total_survivors, na.rm=T)) %>% 
  ungroup() %>% 
  rename(date = date_month) %>% 
  select(date,route,total_dead_and_missing,total_survivors)


plot(filter(df_mmp_daily_med, route=="CMR")$date,
     filter(df_mmp_daily_med, route=="CMR")$total_dead_and_missing, type="l", main="Daily deaths")

plot(filter(df_mmp_monthly_med, route=="CMR")$date,
     filter(df_mmp_monthly_med, route=="CMR")$total_dead_and_missing, type="l", main="Monthly deaths")

plot(log(filter(df_mmp_monthly_med, route=="CMR")$total_survivors),
     log(filter(df_mmp_monthly_med, route=="CMR")$total_dead_and_missing), main="Monthly deaths vs survivors log scale")
