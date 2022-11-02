# Data Analysis and Calculations for Seasonal Change 
# Emma Gemal, emmagemal@outlook.com
# University of Edinburgh, Stockholm University 

### Library ----
library(tidyverse)

### NP, DR and GP Seasonal Changes ----
season <- read.csv("Data/Seasons/season_full.csv")

str(season)
season <- season %>% 
            mutate(date = as.Date(date, format = "%d/%m/%Y"))
str(season)

# model creation 
season_dr <- lm(DR ~ date, data = season)
season_np <- lm(NP ~ date, data = season)
season_gp <- lm(GP ~ date, data = season)

plot(season_dr)   
hist(resid(season_dr))  # not super normally distributed residuals, but it's fine 

plot(season_np)   
hist(resid(season_np))

plot(season_gp)   
hist(resid(season_gp))

## Results 
# DR
summary(season_dr)   # slope = 0.122 ± 0.021
                     # p = 1.05e-5, t = 5.68 (significant)
                     # adjusted R2 = 0.576
anova(season_dr)     # F = 32.21, p = 1.045e-5 

# NP
summary(season_np)   # slope = -0.013 ± 0.0065
                     # p = 0.061, t = -1.97 (NOT significant)
                     # adjusted R2 = 0.112
anova(season_np)     # F = 3.89, p = 0.0614

# GP
summary(season_gp)   # slope = 0.109 ± 0.0025
                     # p = 2.19e-4, t = 4.42 (significant)
                     # adjusted R2 = 0.446
anova(season_gp)     # F = 19.498, p = 0.00022 


### Climate and Chlorophyll Changes ----
## Chlorophyll across the season
chl_season <- read.csv("Data/Seasons/chl_season.csv")

str(chl_season)

chl_season <- chl_season %>% 
                mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
                dplyr::select(date, chl_mg, chl_mmol)
str(chl_season)

# model 
chl_lm <- lm(chl_mg ~ date, data = chl_season)

plot(chl_lm)   
hist(resid(chl_lm))

# results 
summary(chl_lm)   # slope = -8.846
                  # p = 0.000112 (significant)
                  # adjusted R2 = 0.3577
anova(chl_lm)  # F = 19.374, p = 0.000112 


## Temperature across the season 
climate <- read.csv("Data/Seasons/climate_combo.csv")

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
