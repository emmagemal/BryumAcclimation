# Plots for Hallett publication 
# Emma Gemal, s1758915@sms.ed.ac.uk
# University of Edinburgh 

### Library ----
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggpubr)
library(plotrix)
library(gridExtra)
library(patchwork)

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

# plotting curves using surface area  
(sa_plot <- ggplot(avgdata, aes(x = temp, y = avgSA, color = treatment_type)) +
                geom_hline(yintercept = 0, color = "grey", size = 0.8) +  # optional to keep              
                geom_point(aes(shape = type), size = 2) +
                geom_line(aes(linetype = type)) +
                geom_errorbar(aes(ymin = avgSA-se_SA, ymax = avgSA+se_SA), width = 0.5) +
                ylab(label = 
                       expression(Assimilation~per~area~(µmol~m^-2~s^-1))) +
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

ggsave("Figures/t_response_sa.png", plot = sa_plot, 
       width = 6.8, height = 5.2, units = "in")


### Carbon Gain Efficiency ----
cgain <- read.csv("Data/c_gainSA.csv") 
str(cgain)

# making stacked plots of DR:NP 
(stacked_SA <- ggplot(cgain, aes(x = temp, y = percent, fill = type)) +
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

ggsave("Figures/c_gain_SA.png", plot = stacked_SA, 
       width = 7.5, height = 5, units = "in")

# using dry weight instead 
cgain_dw <- read.csv("Data/c_gain_long.csv") 
str(cgain_dw)

# making stacked plots of DR:NP 
(stacked_DW <- ggplot(cgain_dw, aes(x = temp, y = percent, fill = type)) +
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

ggsave("Figures/c_gain_DW.png", plot = stacked_DW, 
       width = 7.5, height = 5, units = "in")


grid.arrange(stacked_DW, stacked_SA, nrow = 2)


### Acclimation Ratios ----
ratios <- read.csv("Data/acclimation_ratios.csv")

(DWratio_plot <- ggplot(ratios, aes(x = temp, y = DWt.c)) +
                    geom_point(aes(color = type, shape = type), size = 2.5, alpha = 0.9) +
                    geom_hline(yintercept = 1, linetype = "dotted") +
                    ylab(label = "DW Acclimation ratio") +
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

ggsave("Figures/DWacclim_ratio_plot.png", plot = DWratio_plot, 
       width = 6.5, height = 5, units = "in")


(SAratio_plot <- ggplot(ratios, aes(x = temp, y = SAt.c)) +
                    geom_point(aes(color = type, shape = type), size = 2.5, alpha = 0.9) +
                    geom_hline(yintercept = 1, linetype = "dotted") +
                    ylab(label = "SA Acclimation ratio") +
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

ggsave("Figures/SAacclim_ratio_plot.png", plot = SAratio_plot, 
       width = 6.5, height = 5, units = "in")

grid.arrange(SAratio_plot, DWratio_plot, ncol = 2)


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


### Seasonal Changes ----
## NP, DR and GP
season_long <- read.csv("Data/season_long.csv")

str(season_long)
season_long <- season_long %>% 
                  mutate(date = as.Date(date, format = "%d/%m/%Y"))
str(season_long)

(season_plot <- ggplot(season_long, aes(x = date, y = rate, group = type)) +
                  geom_point(aes(color = type), size = 2.2) +
                  geom_line(aes(color = type)) +
                  geom_hline(aes(yintercept = 0), linetype = "dotted") +
                  geom_errorbar(aes(ymin = rate-se, ymax = rate+se, color = type), 
                                width = 1) +
                  ylab(label = expression(paste(
                    "CO"[2], " exchange (µmol ", "m"^-2, " s"^-1, ")"))) +  
                  xlab(label = "Date") +
                  theme_bw() +
                  theme(axis.title.x = 
                          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                        axis.title.y = 
                          element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                        panel.grid.minor = element_blank()) +
                  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
                  scale_color_manual(values = c("#FF6D33", "#B5BA4F", "#7A292A"),
                                     name = "Type"))

ggsave("Figures/season_rates.png", plot = season_plot, 
       width = 6.8, height = 5.2, units = "in")

## Chlorophyll content
chl_season <- read.csv("Data/chl_season.csv")

str(chl_season)
chl_season <- chl_season %>% 
                mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
                dplyr::select(date, chl_mg, chl_mmol) 
                
str(chl_season)

chl_sum <- chl_season %>% 
              group_by(date) %>% 
              mutate(chl_mg_se = std.error(chl_mg),
                     chl_mmol_se = std.error(chl_mmol)) %>% 
              summarise(chl_mg = mean(chl_mg),
                        chl_mmol = mean(chl_mmol),
                        chl_mg_se = mean(chl_mg_se),
                        chl_mmol_se = mean(chl_mmol_se)) %>% 
              mutate(chl_µg = chl_mg*1000) %>% 
              na.omit()

(chl_season_plot <- ggplot(chl_sum, aes(x = date, y = chl_mg)) +
                      geom_point(fill = "#5D5F25", color = "#5D5F25", size = 2.2) +
                      geom_line(color = "#5D5F25") +
                      geom_errorbar(aes(ymin = chl_mg-chl_mg_se, ymax = chl_mg+chl_mg_se), 
                                    width = 1, color = "#5D5F25") +
                      ylab(label = expression(Chlorophyll~content~(mg~m^-2))) +  
                      xlab(label = "Date") +
                      theme_bw() +
                      theme(axis.title.x = 
                              element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                            axis.title.y = 
                              element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                            panel.grid.minor = element_blank()) +
                      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")))

ggsave("Figures/season_chl.png", plot = chl_season_plot, 
       width = 6.2, height = 5.2, units = "in")

# creating a panel of the 2 plots 
panel <- grid.arrange(season_plot, chl_season_plot, ncol = 2, padding = 0, widths = c(1.3, 1))
ggsave("Figures/season_panel.png", plot = panel, width = 12, height = 5.5, units = "in")


# 2 y-axis plot instead of panel 
combo <- read.csv("Data/season_chl_combo.csv")

str(combo)
combo <- combo %>% 
            mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
            na.omit()
str(combo)

coeff <- 100
(mixed_plot <- ggplot(combo, aes(x = date)) +
                  geom_point(aes(y = rate, color = type)) +
                  geom_line(aes(y = rate, color = type)) +
                  geom_point(aes(y = chl_mg / coeff), color = "#5D5F25") +
                  geom_line(aes(y = chl_mg / coeff), color = "#5D5F25") +
                  xlab(label = "Date") +
                  theme_bw() +
                  scale_y_continuous(name = expression(paste(
                                            "CO"[2], " exchange (µmol ", "m"^-2, " s"^-1, ")")),
                                     sec.axis = sec_axis(~.*coeff, 
                                                         name = 
                                                         expression(paste(
                                                "Chlorophyll content (mg ", "m"^-2, ")")))) +
                  theme(axis.title.x = 
                          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                        axis.title.y = 
                          element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
                        axis.title.y.right = element_text(vjust = 3), 
                        panel.grid.minor = element_blank()) +
                  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
                  scale_color_manual(values = c("#5D5F25", "#7A292A", "#B5BA4F", "#FF6D33"),
                                     name = "Type",
                                     labels = c("Chlorophyll", "DR", "GP", "NP")))

ggsave("Figures/combo_season.png", plot = mixed_plot, 
       width = 6.8, height = 5.2, units = "in")


### Climate Data ----
climate <- read.csv("Data/climate_combo.csv")

str(climate)
climate$date_time <- as.POSIXct(climate$date_time)
str(climate)

# plotting only 2004/2005 season
climate <- climate %>% 
              filter(season == "2004/5")
summary(climate)

(temp2004 <- ggplot(climate, aes(x = date_time, y = temp)) +   
                geom_line(color = "#FDB321") +
                geom_smooth(method = "lm", se = TRUE, color = "#CA8702") +
                ylab(label = "Temperature (˚C)") +
                xlab(label = "Date") +
                theme_bw() +
                theme(legend.position = "none",
                      axis.title.x = 
                        element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                      axis.title.y = 
                        element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
                theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
                scale_x_datetime(date_labels = "%b", date_breaks = "1 month"))

ggsave("Figures/temp_plot.png", plot = temp2004, 
       width = 7, height = 5.2, units = "in")

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
