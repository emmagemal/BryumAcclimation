# Data Analysis and Calculations for Seasonal Change 
# Emma Gemal, emmagemal@outlook.com
# University of Edinburgh, Stockholm University 

### Library ----
library(tidyverse)

### Seasonal Change Models ----
## NP, DR and GP across the season
season <- read.csv("Data/season_full.csv")

str(season)
season <- season %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y"))
str(season)

season_dr <- lm(DR ~ date, data = season)
season_np <- lm(NP ~ date, data = season)
season_gp <- lm(GP ~ date, data = season)

plot(season_dr)   
hist(resid(season_dr))  # not super normally distributed residuals, but it's fine 

plot(season_np)   
hist(resid(season_np))

plot(season_gp)   
hist(resid(season_gp))

# DR
summary(season_dr)   # slope = 0.122 ± 0.021
                     # p = 1.05e-5, t = 5.68 (significant)
anova(season_dr)

# NP
summary(season_np)   # slope = -0.013 ± 0.0065
                     # p = 0.061, t = -1.97 (NOT significant)
anova(season_np)

# GP
summary(season_gp)   # slope = 0.109 ± 0.0025
                     # p = 2.19e-4, t = 4.42 (significant)
anova(season_gp)


### Models for Climate (Temperature) ----
## Chlorophyll across the season
chl_season <- read.csv("Data/chl_season.csv")

str(chl_season)
chl_season <- chl_season %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
  dplyr::select(date, chl_mg, chl_mmol)
str(chl_season)

chl_lm <- lm(chl_mg ~ date, data = chl_season)

plot(chl_lm)   
hist(resid(chl_lm))

summary(chl_lm)
anova(chl_lm)


climate <- read.csv("Data/climate_combo.csv")

str(climate)
climate$date_time <- as.POSIXct(climate$date_time)
str(climate)

# only 2004/2005 season
climate <- climate %>% 
  filter(season == "2004/5")
summary(climate)


## Simple linear model 
temp_lm <- lm(temp ~ date_time, data = climate)
summary(temp_lm)  # 2.06e-6, p = <2e-16, very significant
                  # adjusted R2 = 0.574

# checking model assumptions 
plot(temp_lm)  # seems normally distributed 
hist(resid(temp_lm))  
bptest(temp_lm)  # there is heteroskedasticity in the model though
