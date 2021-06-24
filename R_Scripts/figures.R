# Plots for Hallett publication 
# Emma Gemal, s1758915@sms.ed.ac.uk
# University of Edinburgh 

### Library ----
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggpubr)

### Temperature Response Curves ----
avgdata <- read.csv("Data/np_dr_averages.csv", header = TRUE)
str(avgdata)

avgdata <- avgdata %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type))
str(avgdata)

# plotting temperature response curves using dry weight  
(dw_plot <- ggplot(avgdata, aes(x = temp, y = avgDW, color = treatment_type)) +
              geom_hline(yintercept = 0, color = "grey", size = 0.8) +  # optional to keep              
              geom_point(aes(shape = type), size = 2) +
              geom_line(aes(linetype = type)) +
              ylab(label = expression(Assimilation~per~dry~weight~(nmol~g^-1~s^-1))) +
              xlab(label = "Temperature (˚C)") +              
              geom_errorbar(aes(ymin = avgDW-se_DW, ymax = avgDW+se_DW), width = 0.5) +
              theme_bw() +
              theme(axis.title.x = 
                      element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                    axis.title.y = 
                      element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                    panel.grid.minor = element_blank()) +
              theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
              scale_color_manual(values = c("#12A7B8", "#004452"),
                                 name = c("Treatment Type", "Process"),
                                 labels = c("Control", "Treatment")) +
              scale_linetype_discrete(name = c("Process", "Treatment Type")) +
              scale_shape_discrete(name = c("Process", "Treatment Type")))  

ggsave("Figures/t_response_DW.png", plot = dw_plot, 
       width = 6.8, height = 5.2, units = "in")

# plotting curves using chlorophyll content 
(chl_plot <- ggplot(avgdata, aes(x = temp, y = avgChl, color = treatment_type)) +
                geom_hline(yintercept = 0, color = "grey", size = 0.8) +  # optional to keep              
                geom_point(aes(shape = type), size = 2) +
                geom_line(aes(linetype = type)) +
                geom_errorbar(aes(ymin = avgChl-se_Chl, ymax = avgChl+se_Chl), width = 0.5) +
                ylab(label = 
                       expression(Assimilation~per~chlorophyll~content~(nmol~mg^-1~s^-1))) +
                xlab(label = "Temperature (˚C)") +              
                theme_bw() +
                theme(axis.title.x = 
                        element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                      axis.title.y = 
                        element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                      panel.grid.minor = element_blank()) +
                theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
                scale_color_manual(values = c("#12A7B8", "#004452"),
                                   name = c("Treatment Type", "Process"),
                                   labels = c("Control", "Treatment")) +
                scale_linetype_discrete(name = c("Process", "Treatment Type")) +
                scale_shape_discrete(name = c("Process", "Treatment Type")) +
                scale_y_continuous(breaks = seq(-20, 10, 5)))

ggsave("Figures/t_response_chl.png", plot = chl_plot, 
       width = 6.8, height = 5.2, units = "in")


### Carbon Gain Efficiency ----
cgain <- read.csv("Data/c_gain_long.csv") 

str(cgain)

# making stacked plots of DR:NP 
(stacked_both <- ggplot(cgain, aes(x = temp, y = percent, fill = type)) +
                    geom_bar(position = "fill", stat = "identity") +
                    ylab(label = "Carbon Use Efficiency (%)") +
                    xlab(label = "Temperature (˚C)") +              
                    facet_wrap(~treatment_type, nrow = 1) +
                    theme_bw() +
                    theme(axis.title.x = 
                            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = 
                            element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
                    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
                          panel.spacing = unit(1, "cm")) + 
                    scale_fill_manual(values = c("#FF6D33", "#7A292A"),
                                      name = "Process") +  # could change the name 
                    scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                                       labels = scales::percent_format(suffix = "")))  # OR
                #    scale_y_continuous(expand = c(0,0)))

ggsave("Figures/c_gain_stacked.png", plot = stacked_both, 
       width = 7.5, height = 5, units = "in")


### Acclimation Ratios ----
ratios <- read.csv("Data/acclimation_ratios.csv")

(ratio_plot <- ggplot(ratios, aes(x = temp, y = DWt.c)) +
                  geom_point(aes(color = type, shape = type), size = 2.5, alpha = 0.9) +
                  geom_hline(yintercept = 1, linetype = "dotted") +
                  ylab(label = "Acclimation ratio") +
                  xlab(label = "Temperature (˚C)") +
                  theme_bw() +
                  theme(axis.title.x = 
                          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                        axis.title.y = 
                          element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                        panel.grid.minor = element_blank()) +
                  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
                  scale_color_manual(values = c("#FF6D33", "#7A292A"),
                                     name = "Process") +
                  scale_shape_discrete(name = "Process") + # can change the name
                  scale_y_continuous(limits = c(0, 6)))

ggsave("Figures/acclim_ratio_plot.png", plot = ratio_plot, 
       width = 6.5, height = 5, units = "in")


### Light Response Curves ----
light <- read.csv("Data/full_lightresponses_revised.csv") 
str(light)

# calculating averages, standard deviation and standard error
light_sum <- light %>% 
                group_by(Lcuv, treatment_type) %>% 
                mutate(seCO2 = std.error(CO2)) %>% 
                summarise(avgCO2 = mean(CO2),
                          seCO2 = mean(seCO2)) %>% 
                na.omit()

light_sum2 <- light %>% 
  group_by(Lcuv, treatment_type) %>% 
  summarise(avgCO2 = mean(CO2),
            sdCO2 = sd(CO2)) %>%    # if I want to SD instead of SE
  na.omit()


# plotting the light response curves
(light_plots <- ggplot(light_sum, aes(x = Lcuv, y = avgCO2)) +
                  geom_hline(yintercept = 0, size = 0.5, linetype = "dotted") +               
                  geom_point(aes(color = treatment_type), size = 2.2) +
                  geom_line(aes(color = treatment_type)) +
                  geom_errorbar(aes(ymin = avgCO2-seCO2, ymax = avgCO2+seCO2, 
                                    color = treatment_type, width = 20), alpha = 0.8) + 
                  ylab(label = expression(paste(
                    "\u0394", "CO"[2], " (rel. ppm)"))) +  # check units 
                  xlab(label = expression(paste(
                    "PPFD ", "(µmol ", "m"^-2, " s"^-1, ")"))) +
                  theme_bw() +
                  theme(panel.grid.minor = element_blank(),
                        axis.title.x = 
                          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
                  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
                  scale_color_manual(values = c("#12A7B8", "#004452"),
                                     name = "Treatment Type",
                                     labels = c("Control", "Treatment"))) 

ggsave("Figures/light_response.png", plot = light_plots, 
       width = 7, height = 5.5, units = "in")


### Climate Data ----
climate <- read.csv("Data/climate_combo.csv")

# creating summary data 
# for year
sum_year <- climate %>% 
  group_by(year) %>% 
  summarise(avg_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp))

sum_year_long <- sum_year %>% 
  pivot_longer(cols = c(2:4),
               names_to = "type",
               values_to = "temp")

# for season
sum_season <- climate %>% 
  group_by(season) %>% 
  summarise(avg_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp))

sum_season_long <- sum_season %>% 
  pivot_longer(cols = c(2:4),
               names_to = "type",
               values_to = "temp")


# plotting the average, minimum and maximum temperature over time (by year) 
minor <- seq(2004, 2016, by = 4)   # making minor gridlines for the plot 

(temp_year <- ggplot(sum_year_long, aes(x = year, y = temp, 
                                        color = type, shape = type)) +
    geom_vline(xintercept = minor, color = "grey92") +                 
    geom_point(size = 2.5) +  
    geom_line(aes(group = type)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Year") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_color_manual(values = c("#004452", "#B5BA4F", "#FF6D33"),
                       labels = c("Mean", "Maximum", "Minimum")) + 
    scale_shape_discrete(labels = c("Mean", "Maximum", "Minimum")) +
    scale_x_continuous(n.breaks = 8))

ggsave("Figures/climate_plot_year.png", plot = temp_year, 
       width = 6.5, height = 5.5, units = "in")

# adding the mixed effects model to the plot 
temp_month_season <- lmer(temp ~ year + (1|month) + (1|season), data = climate, REML = F)

# creating model predictions 
pred.mm <- ggpredict(temp_month_season, terms = c("year"))  

# plotting the predictions 
(temp_year_model <- ggplot(pred.mm) + 
    geom_vline(xintercept = minor, color = "grey92") +                   
    geom_line(aes(x = x, y = predicted), color = "#004452") +   # slope
    geom_ribbon(aes(x = x, ymin = conf.low, 
                    ymax = conf.high), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(size = 2.5, data = sum_year_long, aes(x = year, y = temp, 
                                                     color = type, shape = type)) +
    geom_line(data = sum_year_long, aes(x = year, y = temp, 
                                        color = type, group = type)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Year") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_color_manual(values = c("#004452", "#B5BA4F", "#FF6D33"),
                       labels = c("Mean", "Maximum", "Minimum")) + 
    scale_shape_discrete(labels = c("Mean", "Maximum", "Minimum")) +
    scale_x_continuous(n.breaks = 8))

ggsave("Figures/climate_plot_model.png", plot = temp_year_model, 
       width = 6.5, height = 5.5, units = "in")

## For the appendix
# plotting the average, minimum and maximum temperature over time (by season) 
(temp_season <- ggplot(sum_season_long, aes(x = season, y = temp, 
                                            color = type, shape = type)) +
    geom_point(size = 2.5) +  
    geom_line(aes(group = type)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Season") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_color_manual(values = c("#004452", "#B5BA4F", "#FF6D33"),
                       labels = c("Mean", "Maximum", "Minimum")) + 
    scale_shape_discrete(labels = c("Mean", "Maximum", "Minimum")))

ggsave("Figures/climate_plot_season.png", plot = temp_season, 
       width = 6.5, height = 5.5, units = "in")

# creating a faceted plot of temperature over time
str(climate)
climate$date_time <- as.POSIXct(climate$date_time)

(temp_facet <- ggplot(climate, aes(x = date_time, y = temp, fill = season)) +   
    geom_line(aes(color = season)) +
    facet_wrap(vars(season), ncol = 3, scales = "free_x") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Date") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_x_datetime(date_labels = "%b", date_breaks = "1 month") +
    scale_color_manual(values = c("#DADDA7", "#D1D491", "#C8CC7B", "#BEC365", 
                                  "#B5BA4F", "#A2A741", "#8B8F38", "#74772E", 
                                  "#5D5F25", "#45481C", "#2E3013", "#171809", 
                                  "#000000")))

ggsave("Figures/climate_facet.png", plot = temp_facet, 
       width = 7, height = 7, units = "in")

### Microclimate Data ----
microlog <- read.csv("Data/microclimate.csv")

str(microlog)
microlog <- microlog %>% 
  dplyr::select("Date", "Time", "Inside.temp", "Outside.temp", "ground.temp",
                "inside.humidity", "Outside.humidity") %>% 
  rename(Ground.temp = ground.temp,
         Inside.humidity = inside.humidity) %>% 
  unite("date_time", c("Date", "Time"), remove = FALSE, sep = " ") %>% 
  mutate(date_time = as.POSIXct(date_time, format = "%d/%m/%Y %H:%M")) %>% 
  mutate(Inside.temp = as.numeric(Inside.temp)) %>% 
  mutate(Outside.temp = as.numeric(Outside.temp)) %>% 
  mutate(Ground.temp = as.numeric(Ground.temp)) %>% 
  mutate(Inside.humidity = as.numeric(Inside.humidity)) %>% 
  mutate(Outside.humidity = as.numeric(Outside.humidity))

# plotting outside temperatures
(outside_temp <- ggplot(microlog, aes(x = date_time, y = Outside.temp)) +
    geom_line(size = 0.7, color = "#FF6D33") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Date") +
    theme_bw() +
    theme(axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")))

# plotting ground temperatures
(ground_temp <- ggplot(microlog, aes(x = date_time, y = Ground.temp)) +
    geom_line(size = 0.7, color = "#E64100") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Date") +
    theme_bw() +
    theme(axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")))

# plotting inside temperatures
(inside_temp <- ggplot(microlog, aes(x = date_time, y = Inside.temp)) +
    geom_line(size = 0.7, color = "#FF9166") +
    ylab(label = "Temperature (˚C)") +
    xlab(label = "Date") +
    theme_bw() +
    theme(axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")))

# plotting outside relative humidity (DON'T INCLUDE)
(outside_rh <- ggplot(microlog, aes(x = date_time, y = Outside.humidity)) +
    geom_line())

# plotting inside relative humidity 
(inside_rh <- ggplot(microlog, aes(x = date_time, y = Inside.humidity)) +
    geom_line())

panel <- ggarrange(inside_temp, outside_temp, ground_temp, labels = c("A", "B", "C"),
                   nrow = 1)
panel2 <- ggarrange(inside_temp, outside_temp, ground_temp, labels = c("A", "B", "C"),
                    nrow = 2, ncol = 2)

ggsave("Figures/microclimate_panel_tall.png", plot = panel2, width = 8, height = 8, units = "in")


### Water Content Curves ----
water <- read.csv("Data/water_content_full.csv")
str(water)
water <- water %>% 
  pivot_longer(cols = c(1:2),
               names_to = "type",
               names_prefix = "CO2_",
               values_to = "CO2")

(water_curves <- ggplot(water, aes(x = weight, y = CO2)) +
    geom_point(aes(color = type, shape = type), size = 2.5) +
    geom_line(aes(color = type), size = 0.5) +
    facet_wrap(~treatment_type, scales = "free") +
    ylab(label = expression(paste("\u0394", "CO"[2], " (rel. ppm)"))) +
    xlab(label = "Weight (g)") +
    theme_bw() +
    theme(axis.title.x = 
            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = 
            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    scale_color_manual(values = c("#FF6D33", "#7A292A"),
                       name = "Process") +
    scale_shape_discrete(name = "Process"))

ggsave("Figures/water_content.png", plot = water_curves, width = 9, height = 6, units = "in")
