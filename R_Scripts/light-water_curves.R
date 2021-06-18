# Light curves & water curves for thermal acclimation of B. argenteum var. muticum
# by Emma Gemal, emmagemal@gmail.com
# The University of Edinburgh
# Last edited: 16/06/2021


## Library
library(tidyverse)
library(retistruct)

# loading the data
light <- read.csv("Data/lightcurves.csv", header = TRUE)
water_control <- read.csv("Data/watercurves.csv", header = TRUE)
water_treatment <- read.csv("Data/watercurves_programmed.csv")

weights_c <- read.csv("Data/WC_control_weights2.csv")


### Data Manipulation ----
## Light curves 
str(light)
summary(light)

light <- light %>% 
            dplyr::select(Date, Time, Object, dCO2MP, Area, PARtop, A) %>% 
            mutate(Object = as.factor(Object))

str(light)

light <- light %>% 
            mutate(name = case_when(Object == "3" ~ "T3",
                                    Object == "4" ~ "T4",
                                    Object == "5" ~ "C1",
                                    Object == "6" ~ "C2")) %>% 
            mutate(treatment_type = case_when(Object == "3" ~ "treatment",
                                              Object == "4" ~ "treatment",
                                              Object == "5" ~ "control",
                                              Object == "6" ~ "control"))

## Water curves (control)
str(water_control)
str(weights_c)

water_control <- water_control %>% 
                    filter(Object > 4) %>%    # getting rid of ZP readings, and treatments
                    dplyr::select(Date, Time, Object, Area, dCO2MP, PARtop, A) %>% 
                    mutate(Object = as.factor(Object))
str(water_control)

weights_c <- weights_c %>% 
                mutate(Object = as.factor(Object))
str(weights_c)


water_joined <- left_join(water_control, weights_c)
water_joined$Weight_g[1] <- 6.0548

water_joined <- water_joined %>% 
                  mutate(type = case_when(PARtop > 400 ~ "NP",
                                          PARtop < 400 ~ "DR"))

## Water curves (treatment)
str(water_treatment)

water_treatment <- water_treatment %>% 
                      filter(Object > 0) %>%    # getting rid of ZP readings, and treatments
                      dplyr::select(Date, Time, Object, Area, dCO2MP, PARtop, A) %>% 
                      mutate(Object = as.factor(Object)) %>% 
                      mutate(Date = as.Date(Date, format = "%d/%m/%Y"))
str(water_treatment)


### Light Curve Plot ----
(facet_light <- ggplot(light, aes(x = PARtop, y = A)) +
                  geom_point(aes(color = Object)) +
                  geom_line(aes(color = Object)) +
                  facet_wrap(~Object))

### LSP Determination ----
# Light saturation points (maximum assimilation)
summary(light$A[light$Object == "3"])  # max T3 (object 3) = 0.9754
summary(light$A[light$Object == "4"])  # max T4 (object 4) = 4.7267
summary(light$A[light$Object == "5"])  # max C1 (object 5) = 6.4634
summary(light$A[light$Object == "6"])  # max C2 (object 6) = 4.6967

# 90% LSP
0.9754*0.9    # T3 (object 3) = 0.87786
4.7267*0.9    # T4 (object 4) = 4.25403
6.4634*0.9    # C1 (object 5) = 5.81706
4.6967*0.9    # C2 (object 6) = 4.22703

## Determining the intersection points for 90% LSP
# visualizing the 90% light saturation points 
hline_light <- data.frame(z = c(0.87786, 4.25403, 5.81706, 4.22703), 
                          Object = factor(c("3", "4", "5", "6")))

(light_lsp <- facet_light +
                geom_hline(data = hline_light, aes(yintercept = z)))

# T3 intersection point
t3_a <- c(999.7, 0.8722355)
t3_b <- c(1499.8, 0.9753883)
t3_y_a <- c(900, 0.87786)
t3_y_b <- c(1500, 0.87786)

line.line.intersection(t3_a, t3_b, t3_y_a, t3_y_b, 
                       interior.only = FALSE)               # T3 = 1026.97 µmol m^-2 s^-1

# T4 intersection point
t4_a <- c(1499.0, 4.1941620)
t4_b <- c(2000.5, 4.7267050)
t4_y_a <- c(1200, 4.25403)
t4_y_b <- c(2000, 4.25403)

line.line.intersection(t4_a, t4_b, t4_y_a, t4_y_b, 
                       interior.only = FALSE)               # T4 = 1555.38 µmol m^-2 s^-1

# C1 intersection point
c1_a <- c(999.2, 5.3655250)
c1_b <- c(1498.8, 6.4634280)
c1_y_a <- c(900, 5.81706)
c1_y_b <- c(1500, 5.81706)

line.line.intersection(c1_a, c1_b, c1_y_a, c1_y_b, 
                       interior.only = FALSE)               # C1 = 1204.67 µmol m^-2 s^-1

# C2 intersection point
c2_a <- c(999.4, 3.8976100)
c2_b <- c(1499.3, 4.6967450)
c2_y_a <- c(900, 4.22703)
c2_y_b <- c(1500, 4.22703)

line.line.intersection(c2_a, c2_b, c2_y_a, c2_y_b, 
                       interior.only = FALSE)               # C2 = 1205.47 µmol m^-2 s^-1


### LCP Determination ----
# visualizing the light compensation points 
(light_lcp <- facet_light +
                geom_hline(yintercept = 0))

## Determining the intersection points for LCP
# T3 intersection point
t3_c <- c(101.0, -0.4919699)
t3_d <- c(200.2, 0.0779900)
t3_y_c <- c(0, 0)
t3_y_d <- c(500, 0)

line.line.intersection(t3_c, t3_d, t3_y_c, t3_y_d, 
                       interior.only = FALSE)               # x = 186.63 µmol m^-2 s^-1

# T4 intersection point
t4_c <- c(100.2, -1.1950990)
t4_d <- c(200.4, 0.1431190)
t4_y_c <- c(0, 0)
t4_y_d <- c(500, 0)

line.line.intersection(t4_c, t4_d, t4_y_c, t4_y_d, 
                       interior.only = FALSE)               # x = 189.68 µmol m^-2 s^-1

# C1 intersection point
c1_c <- c(50.2, -0.4950366)
c1_d <- c(100.4, 0.0397000)
c1_y_c <- c(0, 0)
c1_y_d <- c(500, 0)

line.line.intersection(c1_c, c1_d, c1_y_c, c1_y_d, 
                       interior.only = FALSE)               # x = 96.67 µmol m^-2 s^-1

# C2 intersection point
c2_c <- c(50.3, -0.7033754)
c2_d <- c(100.5, 0.1151478)
c2_y_c <- c(0, 0)
c2_y_d <- c(500, 0)

line.line.intersection(c2_c, c2_d, c2_y_c, c2_y_d, 
                       interior.only = FALSE)               # x = 93.44 µmol m^-2 s^-1

### Water Curve Plot ----
(water_curves <- ggplot(water_joined, aes(x = Weight_g, y = A, group = type)) +
                    geom_point(aes(color = type)) +
                    geom_line(aes(group = type, color = type)) +
                    facet_wrap(~Object, scales = "free"))

### Optimum Water Content Determination ----
# Maximum assimilation 
summary(water_joined$A[water_joined$Object == "5"])  # max C1 (object 5) = 4.348 
summary(water_joined$A[water_joined$Object == "6"])  # max C2 (object 6) = 3.9006

# 90% maximum
4.348*0.9    # C1 = 3.9132
3.9006*0.9    # C2 = 3.51054

## Determining the intersection points for 90% maximum
# visualizing the 90% maximum points  
hline_water <- data.frame(z = c(3.9132, 3.51054), 
                          Object = factor(c("5", "6")))

(water_max <- water_curves +
                  geom_hline(data = hline_water, aes(yintercept = z)))

# C1 intersection point (lower)
c1_e <- c(3.6169, 3.8132810)	
c1_f <- c(3.7025, 4.1573930)
c1_y_e <- c(3, 3.9132)
c1_y_f <- c(4, 3.9132)

line.line.intersection(c1_e, c1_f, c1_y_e, c1_y_f, 
                       interior.only = FALSE)               # C1 (lower threshold) = 3.6418 g

# C1 intersection point (upper)
c1_g <- c(3.8682, 4.2020030)	
c1_h <- c(4.6635, 3.2562640)
c1_y_g <- c(3, 3.9132)
c1_y_h <- c(5, 3.9132)

line.line.intersection(c1_g, c1_h, c1_y_g, c1_y_h, 
                       interior.only = FALSE)               # C1 (upper threshold) = 4.1111 g

# C2 intersection point (lower)
c2_e <- c(4.5305, 3.2474550)	
c2_f <- c(4.6645, 3.7807100)
c2_y_e <- c(3, 3.51054)
c2_y_f <- c(5, 3.51054)

line.line.intersection(c2_e, c2_f, c2_y_e, c2_y_f, 
                       interior.only = FALSE)               # C2 (lower threshold) = 4.5966 g

# C2 intersection point (upper)
c2_g <- c(4.7084, 3.9006400)	
c2_h <- c(4.7504, 3.3694930)
c2_y_g <- c(3, 3.51054)
c2_y_h <- c(5, 3.51054)

line.line.intersection(c2_g, c2_h, c2_y_g, c2_y_h, 
                       interior.only = FALSE)               # C2 (upper threshold) = 4.7392 g




