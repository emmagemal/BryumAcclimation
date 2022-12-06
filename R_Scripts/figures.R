# Plots for Hallett publication 
# Emma Gemal, emmagemal@outlook.com
# University of Edinburgh & Stockholm University 

### Library ----
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggpubr)
library(plotrix)  # for std.error()
library(gridExtra)
library(patchwork)
library(DescTools)  # ColToGrey
library(ggh4x)  # for changing facet_grid sizes in ggplot 
library(ggnewscale)  # for editing ggplot2 legend scales (creating groups) 
library(grid)  # for textGrob


### Temperature Response Curves ----
avgdata <- read.csv("Data/Pulse_Experiment/np_dr_averages.csv", header = TRUE)
fulldata <- read.csv("Data/Pulse_Experiment/raw_np_dr_data.csv", header = TRUE)

avgdata <- avgdata %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type))
str(avgdata)

fulldata <- fulldata %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type),
                     sample = as.factor(sample))

# calculating GP
fulldata_wide <- fulldata %>% 
                  dplyr::select(-c(np_DW, np_Chl)) %>% 
                  pivot_wider(names_from = type, values_from = np_SA) %>% 
                  mutate(GP = NP - DR) %>% 
                  dplyr::select(-c(NP, DR)) %>% 
                  group_by(temp, treatment_type) %>% 
                  summarize(avgSA = mean(GP),
                            se_SA = std.error(GP)) %>% 
                  mutate(type = "GP")

avgdata <- full_join(avgdata, fulldata_wide)

avgdata <- avgdata %>% 
              dplyr::select(-c(avgDW, avgChl))


## SA plot 
(sa_plot <- ggplot(avgdata, aes(x = temp, y = avgSA, 
                                group = interaction(type, treatment_type))) +
                geom_hline(yintercept = 0, size = 1, color = "grey80") +
                geom_point(data = subset(avgdata, treatment_type == "control"), 
                                         aes(shape = type, color = type, fill = type), size = 2) +
                geom_smooth(data = subset(avgdata, treatment_type == "control"), 
                            method = "loess", se = F, span = 1.5,
                            aes(linetype = type, color = type)) +
            #    geom_line(aes(linetype = treatment_type)) +
                geom_errorbar(aes(ymin = avgSA-se_SA, ymax = avgSA+se_SA,
                                  color = type), 
                              width = 0.5, linetype = "solid", alpha = 0.7) +
                ylab(label = expression(paste(
                  "DR, NP or GP µmol ", "CO"[2], " m"^"−2", " s"^"−1", ")"))) + 
                xlab(label = "Temperature (˚C)") +  
                annotate("text", label = "    Respiration                Assimilation",
                         x = -4, y = -1, angle = 90, size = 3) +
                coord_cartesian(clip = "off", xlim = c(2, 30), ylim = c(-12, 11)) +
                theme_bw() +
                theme(axis.title.x = 
                        element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                      axis.title.y = 
                        element_text(margin = margin(t = 0, r = 17, b = 0, l = 0)), 
                      panel.grid = element_blank(),
                      legend.key.width = unit(1,"cm"),
                      panel.border = element_rect(size = 0.3)) +
                theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")) + 
                scale_y_continuous(breaks = seq(-20, 10, 5)) +
                scale_x_continuous(n.breaks = 6) +
                scale_linetype_manual(name = "Control",
                                      values = c("solid", "solid", "solid"),
                                      labels = c("DR", "GP", "NP"),
                                      guide = guide_legend(order = 1)) +
                scale_shape_manual(name = "Control", values = c(1, 0, 2), 
                                   labels = c("DR", "GP", "NP"), 
                                   guide = guide_legend(order = 1)) +
                scale_color_manual(name = "Control", values = c("#EA7A0B", "#395493", "#3EACDC"),
                                   labels = c("DR", "GP", "NP"),
                                   guide = guide_legend(order = 1)) +
                scale_fill_manual(name = "Control", values = c("#EA7A0B", "#395493", "#3EACDC"),
                                  labels = c("DR", "GP", "NP"),
                                  guide = guide_legend(order = 1)) +
                new_scale_color() +
                new_scale_fill() +
                new_scale("shape") +
                new_scale("linetype") +
                geom_point(data = subset(avgdata, treatment_type == "treatment"), 
                           aes(shape = type, color = type, fill = type), size = 2) +
                geom_smooth(data = subset(avgdata, treatment_type == "treatment"), 
                            method = "loess", se = F, span = 1.5,
                            aes(linetype = type, color = type)) +
                scale_shape_manual(name = "Treatment", values = c(21, 22, 24), 
                                   labels = c("DR", "GP", "NP"), 
                                   guide = guide_legend(order = 2))+
                scale_color_manual(name = "Treatment", 
                                   values = c("#EA7A0B", "#395493", "#3EACDC"), 
                                   labels = c("DR", "GP", "NP"),
                                   guide = guide_legend(order = 2)) +
                scale_fill_manual(name = "Treatment", values = c("#EA7A0B", "#395493", "#3EACDC"),
                                   labels = c("DR", "GP", "NP"),
                                   guide = guide_legend(order = 2)) +
                scale_linetype_manual(name = "Treatment",
                                      values = c("22", "22", "22"),
                                      labels = c("DR", "GP", "NP"),
                                      guide = guide_legend(order = 2)))
             
#ggsave("Figures/t_response_sa.png", plot = sa_plot, 
#       width = 6, height = 5, units = "in")


### Acclimation Ratios ----
# calculation of AR
ARratio <- as.data.frame(matrix(ncol = 3, nrow = 12))
colnames(ARratio) <- c("temp", "type", "ratio")

# treatment/control
ARratio$ratio <- c((-0.6804083/-0.4568380), (-2.5739912/-0.7737948), (-4.2830330/-1.7684478), # DR
                   (-6.4027841/-3.0684388), (-7.9797452/-4.9741951), (-10.5753313/-6.1260173), # DR
                   (2.4157606/1.7340523), (4.6715594/1.3318238), (4.5930837/1.7676349), # NP
                   (2.7201547/0.8418550), (0.9378886/abs(-0.7544705)), (-2.1095351/-1.9636927)) # NP
ARratio$temp <- c(3.5, 7.5, 12.5, 17.5, 22.5, 27.5, 3.5, 7.5, 12.5, 17.5, 22.5, 27.5)
ARratio$type <- c(rep("DR", times = 6), rep("NP", times = 6))
ARratio <- ARratio[-1,]  # removing outliers


(SAratio_plot <- ggplot(ARratio, aes(x = temp, y = ratio)) +
                    geom_point(aes(color = type, shape = type), 
                               size = 2.5, alpha = 0.9) +
                    geom_hline(yintercept = 1, linetype = "dotted", color = "grey30") +
                    ylab(label = "Acclimation ratio") +
                    xlab(label = "Temperature (˚C)") +
                    theme_bw() +
                    theme(axis.title.x = 
                            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = 
                            element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          panel.grid = element_blank()) +
                    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")) +
                    scale_x_continuous(limits = c(2, 30)) +
                    scale_y_continuous(limits = c(0.9, 4)) +
                    scale_color_manual(values = c("#EA7A0B", "#3EACDC"),
                                       name = "Process") +
                    scale_shape_manual(name = "Process",
                                       values = c(16, 17)))


### Light Response Curves ----
light <- read.csv("Data/Pulse_Experiment/lightresponses_revised.csv") 
str(light)
light <- light %>% 
            mutate(sample_nr = case_when(sample == "C1" ~ "1",
                                         sample == "C2" ~ "2",
                                         sample == "C3" ~ "3",
                                         sample == "T1" ~ "1",
                                         sample == "T2" ~ "2",
                                         sample == "T3" ~ "3")) %>% 
            mutate(sample = as.factor(sample)) %>% 
            mutate(treatment_type = case_when(treatment_type == "control" ~ "Control",
                                              treatment_type == "treatment" ~ "Treatment"))

# creating multiple photosynthetic light response (A-Q) curves 
fit_AQ_curve <- function(df, group_id, Photo, PARi, fit_type = "onls"){
  AQ_curve_fits <- data.frame(ID = character(),
                              Asat = numeric(),
                              Phi = numeric(),
                              Rd = numeric(),
                              theta = numeric(),
                              resid_SSs = numeric(),
                              LCP = numeric(),
                              Q_sat_75 = numeric(),
                              Q_sat_85 = numeric(),  
                              stringsAsFactors = FALSE
  )
  if(fit_type == "onls"){
    if(require("onls")){
      print("onls is loaded correctly")
    } else {
      print("trying to install onls")
      install.packages("onls")
      if(require("onls")){
        print("onls installed and loaded")
      } else {
        stop("could not install onls")
      }
    }
    library("onls")      
    for(i in seq_along(unique(df[[group_id]]))){
      tryCatch({
        AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
        # Subset by group_ID iteratively:
        single_curve1 <- df[df[[group_id]] == unique(df[[group_id]])[i],]
        single_curve1$assim <- single_curve1[[Photo]]
        single_curve1$PAR <- single_curve1[[PARi]]
        single_curve = single_curve1[order(single_curve1$PAR),]
        phi.as.slope <- with(single_curve,
                             as.numeric(coef(lm(
                               assim[1:5] ~ PAR[1:5]))[2]))
        # Fit the curve:
        temp.fit <- with(single_curve, # use the subset of a single curve
                         onls(assim ~ ((Phi * PAR + Asat - 
                                          sqrt((Phi * PAR + Asat)^2 - 
                                                 4 * Phi * theta * 
                                                 Asat * PAR ))
                         )/(2*theta) - Rd,
                         start=list(
                           Asat = (max(assim)),
                           Phi = phi.as.slope,
                           Rd = -min(assim),
                           theta = 0.5),
                         control = list(maxiter = 50)#,
                         #algorithm = "port"
                         )
        )
        AQ_curve_fits[i, 2] <- as.numeric(coef(temp.fit)[1]) # asat 
        AQ_curve_fits[i, 3] <- as.numeric(coef(temp.fit)[2]) # Phi
        AQ_curve_fits[i, 4] <- as.numeric(coef(temp.fit)[3]) # Rd
        AQ_curve_fits[i, 5] <- as.numeric(coef(temp.fit)[4]) # theta
        AQ_curve_fits[i, 6] <- sum(resid(temp.fit)^2)
        AQ_curve_fits[i, 7] <- (as.numeric(coef(temp.fit)[3]) *(
          as.numeric(coef(temp.fit)[3]) * as.numeric(coef(temp.fit)[4]) - 
            as.numeric(coef(temp.fit)[1]))
        ) / (as.numeric(coef(temp.fit)[2]) * (
          as.numeric(coef(temp.fit)[3]) - as.numeric(coef(temp.fit)[1])
        ))
        AQ_curve_fits[i, 8] <- (
          (as.numeric(coef(temp.fit)[1]) * 0.75 + 
             (as.numeric(coef(temp.fit)[3]))) * (
               as.numeric(coef(temp.fit)[1]) * 0.75 *
                 as.numeric(coef(temp.fit)[4]) +
                 as.numeric(coef(temp.fit)[3]) *
                 as.numeric(coef(temp.fit)[4]) -
                 as.numeric(coef(temp.fit)[1])
             )) / (
               as.numeric(coef(temp.fit)[2])* (
                 as.numeric(coef(temp.fit)[1]) * 0.75 +
                   as.numeric(coef(temp.fit)[3]) -
                   as.numeric(coef(temp.fit)[1])
               ))
        
        AQ_curve_fits[i, 9] <- (
          (as.numeric(coef(temp.fit)[1]) * 0.85 + 
             (as.numeric(coef(temp.fit)[3]))) * (
               as.numeric(coef(temp.fit)[1]) * 0.85 *
                 as.numeric(coef(temp.fit)[4]) +
                 as.numeric(coef(temp.fit)[3]) *
                 as.numeric(coef(temp.fit)[4]) -
                 as.numeric(coef(temp.fit)[1])
             )) / (
               as.numeric(coef(temp.fit)[2])* (
                 as.numeric(coef(temp.fit)[1]) * 0.85 +
                   as.numeric(coef(temp.fit)[3]) -
                   as.numeric(coef(temp.fit)[1])
               ))
      }, error = function(E){cat("Error: ", conditionMessage(E), "\n")})
    }
    return(AQ_curve_fits)
  } else{
    if(fit_type == "nls"){
      for(i in seq_along(unique(df[[group_id]]))){
        tryCatch({
          AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
          # Subset by group_ID iteratively:
          single_curve1 <- df[df[[group_id]] == unique(df[[group_id]])[i],]
          single_curve1$assim <- single_curve1[[Photo]]
          single_curve1$PAR <- single_curve1[[PARi]]
          single_curve = single_curve1[order(single_curve1$PAR),]
          phi.as.slope <- with(single_curve,
                               as.numeric(coef(lm(
                                 assim[1:5] ~ PAR[1:5]))[2]))
          # Fit the curve:
          temp.fit <- with(single_curve, 
                           nls(assim ~ ((Phi * PAR + Asat - 
                                           sqrt((Phi * PAR + Asat)^2 - 
                                                  4 * Phi * theta * 
                                                  Asat * PAR ))
                           )/(2*theta) - Rd,
                           start=list(
                             Asat = (max(assim)),
                             Phi = phi.as.slope,
                             Rd = -min(assim),
                             theta = 0.5),
                           control = list(maxiter = 50),
                           algorithm = "port")
          )
          AQ_curve_fits[i, 2] <- as.numeric(coef(temp.fit)[1]) # asat 
          AQ_curve_fits[i, 3] <- as.numeric(coef(temp.fit)[2]) # Phi
          AQ_curve_fits[i, 4] <- as.numeric(coef(temp.fit)[3]) # Rd
          AQ_curve_fits[i, 5] <- as.numeric(coef(temp.fit)[4]) # theta
          AQ_curve_fits[i, 6] <- sum(resid(temp.fit)^2)
          AQ_curve_fits[i, 7] <- (as.numeric(coef(temp.fit)[3]) *(
            as.numeric(coef(temp.fit)[3]) * 
              as.numeric(coef(temp.fit)[4]) - 
              as.numeric(coef(temp.fit)[1]))
          ) / (as.numeric(coef(temp.fit)[2]) * (
            as.numeric(coef(temp.fit)[3]) - 
              as.numeric(coef(temp.fit)[1])
          ))
          AQ_curve_fits[i, 8] <- (
            (as.numeric(coef(temp.fit)[1]) * 0.75 + 
               (as.numeric(coef(temp.fit)[3]))) * (
                 as.numeric(coef(temp.fit)[1]) * 0.75 *
                   as.numeric(coef(temp.fit)[4]) +
                   as.numeric(coef(temp.fit)[3]) *
                   as.numeric(coef(temp.fit)[4]) -
                   as.numeric(coef(temp.fit)[1])
               )) / (
                 as.numeric(coef(temp.fit)[2])* (
                   as.numeric(coef(temp.fit)[1]) * 0.75 +
                     as.numeric(coef(temp.fit)[3]) -
                     as.numeric(coef(temp.fit)[1])
                 ))
          AQ_curve_fits[i, 9] <- (
            (as.numeric(coef(temp.fit)[1]) * 0.85 + 
               (as.numeric(coef(temp.fit)[3]))) * (
                 as.numeric(coef(temp.fit)[1]) * 0.85 *
                   as.numeric(coef(temp.fit)[4]) +
                   as.numeric(coef(temp.fit)[3]) *
                   as.numeric(coef(temp.fit)[4]) -
                   as.numeric(coef(temp.fit)[1])
               )) / (
                 as.numeric(coef(temp.fit)[2])* (
                   as.numeric(coef(temp.fit)[1]) * 0.85 +
                     as.numeric(coef(temp.fit)[3]) -
                     as.numeric(coef(temp.fit)[1])
                 ))
        }, error = function(E){
          cat("Error: ", conditionMessage(E), "\n")})
      }
      return(AQ_curve_fits)      
    } else{print("ERROR: 'fit_type' specified incorrectly.")}
  }
}

my.fits <- fit_AQ_curve(df = light,
                        Photo = "NP_SA", PARi = "Lcuv", group_id = "sample", fit_type = "onls")

str(my.fits)
my.fits <- my.fits %>% 
              mutate(sample = case_when(ID == "1" ~ "C1",
                                        ID == "2" ~ "C2",
                                        ID == "3" ~ "C3",
                                        ID == "4" ~ "T1",
                                        ID == "5" ~ "T2",
                                        ID == "6" ~ "T3"))

# doing it manually using output from 'my.fits'
# C1 
curve.c1 <- function(PARi){
              (0.07287260 * PARi + 8.675147 - 
                sqrt((0.07287260 * PARi + 8.675147)^2 - 4 *
                       0.07287260 * -1.2748040 * PARi *
                       8.675147)
              ) / (2*-1.2748040) - 2.873080
}

par.c1 <- data.frame(Lcuv = 0:1550,
                     curve = curve.c1(PARi = 0:1550),
                     sample = "C1")
# C2
curve.c2 <- function(PARi){
              (0.08371293 * PARi + 4.422790 - 
                sqrt((0.08371293 * PARi + 4.422790)^2 - 4 *
                       0.08371293 * -1.3723880 * PARi *
                       4.422790)
              ) / (2*-1.3723880) - 2.562400
}

par.c2 <- data.frame(Lcuv = 0:1550,
                     curve = curve.c2(PARi = 0:1550),
                     sample = "C2")
# C3
curve.c3 <- function(PARi){
              (0.07674525 * PARi + 13.465161 - 
                 sqrt((0.07674525* PARi + 13.465161)^2 - 4 *
                        0.07674525 * -1.1614485 * PARi *
                        13.465161)
              ) / (2*-1.1614485) - 3.212828
}

par.c3 <- data.frame(Lcuv = 0:1550,
                     curve = curve.c3(PARi = 0:1550),
                     sample = "C3")
# T1
curve.t1 <- function(PARi){
              (0.11458331 * PARi + 13.982050 - 
                 sqrt((0.11458331 * PARi + 13.982050)^2 - 4 *
                        0.11458331 * 0.7072400 * PARi *
                        13.982050)
              ) / (2*0.7072400) - 2.823915
}

par.t1 <- data.frame(Lcuv = 0:1550,
                     curve = curve.t1(PARi = 0:1550),
                     sample = "T1")
# T2
curve.t2 <- function(PARi){
              (0.08046531 * PARi + 10.091202 - 
                 sqrt((0.08046531 * PARi + 10.091202)^2 - 4 *
                        0.08046531 * 0.9372664 * PARi *
                        10.091202)
              ) / (2*0.9372664) - 2.030062
}

par.t2 <- data.frame(Lcuv = 0:1550,
                     curve = curve.t2(PARi = 0:1550),
                     sample = "T2")
# T3
curve.t3 <- function(PARi){
              (0.16754053 * PARi + 18.954037 - 
                 sqrt((0.16754053 * PARi + 18.954037)^2 - 4 *
                        0.16754053 * 0.1385859 * PARi *
                        18.954037)
              ) / (2*0.1385859) - 3.666333
}

par.t3 <- data.frame(Lcuv = 0:1550,
                     curve = curve.t3(PARi = 0:1550), 
                     sample = "T3")

par.all <- rbind(par.c1, par.c2, par.c3, par.t1, par.t2, par.t3)
par.all <- par.all %>% 
              mutate(treatment_type = ifelse(grepl("C", sample), "Control", "Treatment"))

# vertical plot 
(light_curves <- ggplot(par.all, aes(x = Lcuv, y = curve, group = sample)) +
                    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
                    geom_point(data = light, aes(x = Lcuv, y = NP_SA, color = treatment_type,
                                                 shape = treatment_type, fill = treatment_type), 
                               alpha = 0.5, size = 2) +
                    geom_line(aes(color = treatment_type, linetype = treatment_type),
                              size = 0.8) +
                    facet_wrap(~treatment_type, dir = "v") +
                    ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
                    xlab(label = expression(paste(
                      "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
                    theme_bw() +
                    theme(axis.title.x = 
                            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = 
                            element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                          panel.grid = element_blank(),
                          legend.position = "bottom",
                          legend.key.width = unit(1, "cm"),
                          legend.key.height = unit(0.1, "cm"),
                          legend.title = element_text(margin = margin(r = 10))) +
                    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
                    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
                    scale_color_manual(values = c("black", "black"),
                                       name = "Treatment Type",
                                       labels = c("Control", "Treatment")) +
                    scale_fill_manual(values = c("black", "black"),
                                      name = "Treatment Type",
                                      labels = c("Control", "Treatment")) +
                    scale_linetype_manual(name = "Treatment Type",
                                          labels = c("Control", "Treatment"),
                                          values = c("solid", "22")) +
                    scale_shape_manual(name = "Treatment Type",
                                       labels = c("Control", "Treatment"),
                                       values = c(1, 21)) +
                    scale_y_continuous(limits = c(-5, 16)))

#ggsave("Figures/light_response_curves.png", plot = light_curves, 
#       height = 7, width = 4.5, units = "in")

# horizontal plot 
(light_curves2 <- ggplot(par.all, aes(x = Lcuv, y = curve, group = sample)) +
                    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
                    geom_point(data = light, aes(x = Lcuv, y = NP_SA, color = treatment_type,
                                                 shape = treatment_type, fill = treatment_type), 
                               alpha = 0.5, size = 2) +
                    geom_line(aes(color = treatment_type, linetype = treatment_type),
                              size = 0.8) +
                    facet_wrap(~treatment_type) +
                    ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
                    xlab(label = expression(paste(
                      "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
                    theme_bw() +
                    theme(axis.title.x = 
                            element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y = 
                            element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                          panel.grid = element_blank(),
                          legend.key.width = unit(1,"cm"),
                          legend.position = "right") +
                    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
                    scale_color_manual(values = c("black", "black"),
                                       name = "Treatment Type",
                                       labels = c("Control", "Treatment")) +
                    scale_fill_manual(values = c("black", "black"),
                                      name = "Treatment Type",
                                      labels = c("Control", "Treatment")) +
                    scale_linetype_manual(name = "Treatment Type",
                                          labels = c("Control", "Treatment"),
                                          values = c("solid", "22")) +
                    scale_shape_manual(name = "Treatment Type",
                                       labels = c("Control", "Treatment"),
                                       values = c(1, 21)) +
                    scale_y_continuous(limits = c(-5, 16)))


# calculating averages, standard deviation and standard error
light_sum <- light %>% 
                group_by(Lcuv, treatment_type) %>% 
                mutate(seNP = std.error(NP_SA)) %>% 
                summarise(avgNP = mean(NP_SA),
                          seNP = mean(seNP)) %>% 
                na.omit()


### Seasonal Changes ----
## NP, DR and GP
season_long <- read.csv("Data/Seasons/season_long.csv")
str(season_long)
season_long <- season_long %>% 
                  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
                  mutate(type = as.factor(type)) %>% 
                  group_by(date, type) %>% 
                  summarize(weight.avg = mean(weight),
                            se.weight = (sd(weight)/sqrt(5)),
                            rate.avg = mean(rate),
                            se.rate = (sd(rate)/sqrt(5)))

season_long <- season_long %>% 
                  mutate(rate.avg = ifelse(type == "DR", -rate.avg, rate.avg))

# write.csv(season_long,"Data/Seasons/season_long.csv", row.names = FALSE)


(season_plot <- ggplot(season_long, aes(x = date, y = rate.avg, group = type)) +
                  geom_point(aes(color = type, fill = type, shape = type), 
                             size = 3) +
                 # geom_smooth(method = "loess", span = 1, se = FALSE, aes(color = type)) +
                  geom_line(aes(color = type), size = 1) +
                  geom_hline(aes(yintercept = 0), linetype = "dotted") +
                  geom_errorbar(aes(ymin = rate.avg-se.rate, ymax = rate.avg+se.rate, 
                                    color = type), 
                                width = 1, alpha = 0.7) +
                  ylab(label = expression(paste(
                    "NP, DR or GP (µmol ", "CO"[2], " m"^"−2", " s"^"−1", ")"))) +  
                  xlab(label = "Date") +
                  theme_bw() +
                  theme(axis.title.x = 
                          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                        axis.title.y = 
                          element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                        panel.grid = element_blank(),
                        axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
                  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")) + 
                  scale_color_manual(values = c("#EA7A0B", "#395493", "#3EACDC"),
                                     name = "Process") +
                  scale_fill_manual(values = c("#EA7A0B", "#395493", "#3EACDC"),
                                    name = "Process") +
                  scale_shape_manual(name = "Process",
                                     values = c(21, 22, 24)) +
                  scale_x_date(date_labels = "%d/%m/%y"))

#ggsave("Figures/season_rates.png", plot = season_plot, 
#       width = 5.5, height = 4.5, units = "in")


## Chlorophyll content
chl_season <- read.csv("Data/Seasons/chl_season.csv")

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
                      geom_point(fill = "#5D5F25", color = "#5D5F25", size = 3) +
                   #   geom_smooth(method = "loess", se = F, span = 1, color = "#5D5F25") +
                      geom_line(color = "#5D5F25", size = 1, ) +
                      geom_errorbar(aes(ymin = chl_mg-chl_mg_se, ymax = chl_mg+chl_mg_se), 
                                    width = 1, color = "#5D5F25", alpha = 0.7) +
                      ylab(label = expression(paste("Chlorophyll content "("mg m"^"−2")))) +  
                      xlab(label = "Date") +
                      theme_bw() +
                      theme(axis.title.x = 
                              element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                            axis.title.y = 
                              element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
                            panel.grid = element_blank(),
                            axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
                      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.75), "cm")) +
                      scale_x_date(date_labels = "%d/%m/%y"))

#ggsave("Figures/season_chl.png", plot = chl_season_plot, 
#       width = 4.5, height = 4.5, units = "in")


### Climate Data ----
climate <- read.csv("Data/Seasons/climate_combo.csv")

str(climate)
climate$date_time <- as.POSIXct(climate$date_time)
str(climate)

# plotting only 2004/2005 season
climate <- climate %>% 
              filter(season == "2004/5")
summary(climate)

(temp2004 <- ggplot(climate, aes(x = date_time, y = temp)) +   
                geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                geom_line(color = "grey20") +
                ylab(label = "Temperature (˚C)") +
                xlab(label = "Date") +
                theme_bw() +
                theme(legend.position = "none",
                      axis.title.x = 
                        element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                      axis.title.y = 
                        element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                      panel.grid = element_blank()) +
                theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
                scale_x_datetime(date_labels = "%d/%m/%y"))

#ggsave("Figures/temp_plot.png", plot = temp2004, 
#       width = 6, height = 4.5, units = "in")


### Microclimate Data ----
microlog <- read.csv("Data/Pulse_Experiment/microclimate.csv")
str(microlog)
summary(microlog)

microlog <- microlog %>% 
              dplyr::select("Date", "Time", "Inside.temp", "Outside.temp",
                            "inside.humidity", "Outside.humidity") %>% 
              rename(Inside.humidity = inside.humidity) %>% 
              unite("date_time", c("Date", "Time"), remove = FALSE, sep = " ") %>% 
              mutate(date_time = as.POSIXct(date_time, format = "%d/%m/%Y %H:%M")) %>% 
              mutate(Inside.temp = as.numeric(Inside.temp)) %>% 
              mutate(Outside.temp = as.numeric(Outside.temp)) %>% 
              mutate(Inside.humidity = as.numeric(Inside.humidity)) %>% 
              mutate(Outside.humidity = as.numeric(Outside.humidity))
str(microlog)

# temperature only
microlog.t <- microlog %>% 
                dplyr::select(-Inside.humidity, -Outside.humidity) %>% 
                mutate(Temp.diff = Inside.temp-Outside.temp) %>% 
                pivot_longer(cols = c("Inside.temp", "Outside.temp", "Temp.diff"), 
                             names_to = "location", values_to = "temp")

# relative humidity only 
microlog.rh <- microlog %>% 
                dplyr::select(-Inside.temp, -Outside.temp) %>% 
                pivot_longer(cols = c("Inside.humidity", "Outside.humidity"),
                             names_to = "location", values_to = "humidity") %>% 
                filter(!(Date < "27/11/2018" & location == "Outside.humidity"))

# calculating water vapor pressure
vp.fun <- function(temp){ 
            6.1078^((17.269*temp)/(temp + 237.3))  # Tetens formula in millibars
}

vpd.fun <- function(vp, rh){ 
              vp-(rh*vp/100)   # VPD equation in millibars 
}

microlog.vp <- microlog %>% 
                mutate(Inside.VP = vp.fun(temp = Inside.temp)*0.1) %>% 
                mutate(Outside.VP = vp.fun(temp = Outside.temp)*0.1) %>% 
                mutate(Inside.VPD = vpd.fun(vp = Inside.VP, rh = Inside.humidity)*0.1) %>% 
                mutate(Outside.VPD = vpd.fun(vp = Outside.VP, rh = Outside.humidity)*0.1)

microlog.vp <- microlog.vp %>% 
                  dplyr::select(-Inside.temp, -Outside.temp, 
                                -Inside.humidity, -Outside.humidity) %>% 
                  pivot_longer(cols = c("Inside.VP", "Outside.VP", "Inside.VPD", "Outside.VPD"),
                               names_to = "location", values_to = "vapor") %>% 
                  mutate(vapor_type = case_when(location == "Inside.VP" ~ "VP",
                                                location == "Inside.VPD" ~ "VPD",
                                                location == "Outside.VP" ~ "VP",
                                                location == "Outside.VPD" ~ "VPD")) %>% 
                  mutate(location_type = case_when(location == "Inside.VP" ~ "Inside",
                                                   location == "Inside.VPD" ~ "Inside",
                                                   location == "Outside.VP" ~ "Outside",
                                                   location == "Outside.VPD" ~ "Outside"))

## Plots 
# plotting temperatures
(temp <- ggplot(microlog.t, aes(x = date_time, y = temp)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            geom_hline(yintercept = 5, linetype = "dashed", color = "grey80") +
            geom_line(size = 0.7, aes(color = location)) +
            ylab(label = "Temperature (˚C)") +
            xlab(label = "Date") +
            theme_bw() +
            theme(axis.title.x = 
                    element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y = 
                    element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
            theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
            scale_x_datetime(date_labels = "%d/%m/%y") +
            scale_color_manual(name = "Logger Location",
                               label = c("Inside", "Outside", "Difference"),
                               values = c("#DA4D10", "#395493", "#F4B076")))

# plotting relative humidity 
(rh <- ggplot(microlog.rh, aes(x = date_time, y = humidity)) +
          geom_line(aes(color = location), size = 0.7) +
          ylab(label = "Relative humidity (%)") +
          xlab(label = "Date") +
          theme_bw() +
          theme(axis.title.x = 
                  element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = 
                  element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                panel.grid = element_blank(),
                axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
          scale_x_datetime(date_labels = "%d/%m/%y") +
          scale_color_manual(name = "Logger Location",
                             label = c("Inside", "Outside"),
                             values = c("#DA4D10", "#395493")))

# plotting vapor pressure 
(vp <- ggplot(microlog.vp, aes(x = date_time, y = vapor)) +
          geom_line(aes(color = location_type, linetype = vapor_type), size = 0.7) +
          ylab(label = "VP and VPD (kPa)") +
          xlab(label = "Date") +
          facet_wrap(~location_type, scales = "free_y", dir = "v") +
          force_panelsizes(rows = c(1, 0.4)) +  # change relative sizes of the plots 
          theme_bw() +
          theme(axis.title.x = 
                  element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = 
                  element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                panel.grid = element_blank(),
                legend.position = "top",
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
          theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
          scale_y_continuous(n.breaks = 3) +
          scale_x_datetime(date_labels = "%d/%m/%y") +
          scale_linetype_manual(name = "",
                                values = c("solid", "11")) +
          scale_color_manual(name = "",
                             label = c("Inside", "Outside"),
                             values = c("#DA4D10", "#395493")) +
          guides(color = "none"))

# another version of VP (legend on right)
(vp2 <- ggplot(microlog.vp, aes(x = date_time, y = vapor)) +
          geom_line(aes(color = location_type, linetype = vapor_type), size = 0.7) +
          ylab(label = "VP and VPD (kPa)") +
          xlab(label = "Date") +
          facet_wrap(~location_type, scales = "free_y", dir = "v") +
          force_panelsizes(rows = c(1, 0.4)) +  # change relative sizes of the plots 
          theme_bw() +
          theme(axis.title.x = 
                  element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = 
                  element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                panel.grid = element_blank(),
                legend.position = "right",
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
          scale_y_continuous(n.breaks = 3) +
          scale_x_datetime(date_labels = "%d/%m/%y") +
          scale_linetype_manual(name = "",
                                values = c("solid", "11")) +
          scale_color_manual(name = "",
                             label = c("Inside", "Outside"),
                             values = c("#DA4D10", "#395493")) +
          guides(color = "none"))

# another version of VP (facet horizontal)
(vp3 <- ggplot(microlog.vp, aes(x = date_time, y = vapor)) +
          geom_line(aes(color = location_type, linetype = vapor_type), size = 0.7) +
          ylab(label = "VP and VPD (kPa)") +
          xlab(label = "Date") +
          facet_wrap(~location_type) +
          force_panelsizes(rows = c(1, 0.4)) +  # change relative sizes of the plots 
          theme_bw() +
          theme(axis.title.x = 
                  element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = 
                  element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                panel.grid = element_blank(),
                legend.position = "top",
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7)) +
          theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
          scale_y_continuous(n.breaks = 3) +
          scale_x_datetime(date_labels = "%d/%m/%y") +
          scale_linetype_manual(name = "",
                                values = c("solid", "11")) +
          scale_color_manual(name = "",
                             label = c("Inside", "Outside"),
                             values = c("#DA4D10", "#395493")) +
          guides(color = "none"))

### OTC Absorbance Spectra ----
absorbance <- read.csv("Data/Pulse_Experiment/OTC_absorbance_spectra.csv")

absorbance <- absorbance %>% 
                mutate(otc_arbunits2 = otc_arbunits/1000)

(otc <- ggplot(absorbance, aes(x = wavelength_nm, y = otc_arbunits2)) +
          geom_line(color = "grey20", size = 0.7) +
          ylab(label = "Transmittance") +
          xlab(label = "Wavlength (nm)") +
          theme_bw() +
          theme(axis.title.x = 
                  element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                axis.title.y = 
                  element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                axis.text.x = element_text(margin = margin(t = 5)),
                axis.ticks.length.x = unit(.25, "cm"),
                panel.grid = element_blank()) +
          theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +
          scale_x_continuous(limits = c(190, 650)))

(otc2 <- ggplot(absorbance, aes(x = wavelength_nm, y = otc_arbunits2)) +
            geom_line(color = "grey20", size = 0.7) +
            ylab(label = "Arbitrary units") +
            xlab(label = "Wavlength (nm)") +
            theme_bw() +
            theme(axis.title.x = 
                    element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y = 
                    element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  panel.grid = element_blank()) +
            theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +  # different margin
            scale_x_continuous(limits = c(190, 650)))


### Water Content Curves ----
water <- read.csv("Data/Pulse_Experiment/water_content_full.csv")
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


### Main Panel (control vs. treatment, acclimation ratios, light curves) ----
panel.all <- ggarrange(light_curves,
                       ggarrange(sa_plot, 
                                 ggarrange(SAratio_plot, NULL, ncol = 2, widths = c(1, 0.03)),
                                 nrow = 2, labels = c("b", "c"),
                                 heights = c(1, 0.78),
                                 font.label = list(size = 20, face = "bold")),
                       ncol = 2, labels = "a", widths = c(0.75, 1),
                       font.label = list(size = 20, face = "bold")) 
(panel.all <- panel.all + bgcolor("#ffffff"))
ggsave("Figures/pulse_results_panel.png", panel.all, width = 8.2, height = 6, units = "in")

# with light curves at top instead
panel.all.top <- ggarrange(light_curves2,
                           ggarrange(sa_plot2, SAratio_plot, ncol = 2, labels = c("b", "c"),
                                     font.label = list(size = 20, face = "bold"), 
                                     widths = c(1, 0.9)),
                           nrow = 2, heights = c(0.9, 1), labels = "a",
                           font.label = list(size = 20, face = "bold")) 
panel.all.top
#ggsave("Figures/pulse_results_panel3.png", panel.all.top, width = 9, height = 7, units = "in")


### Season Panel (season, chlorophyll, light curves) ----
# panel with season and chlorophyll only
panel.season <- ggarrange(season_plot, chl_season_plot, 
                   ncol = 2, widths = c(1.2, 1), labels = c("a", "b"), 
                   font.label = list(size = 20, face = "bold"))
panel.season
#ggsave("Figures/season_panel.png", plot = panel.season, width = 9, height = 3.8, units = "in")

# panel with temperature 
panel.season2 <- ggarrange(nrow = 2, heights = c(1, 0.8),
                  ggarrange(season_plot, chl_season_plot, 
                            ncol = 2, widths = c(1, 0.75), labels = c("a", "b"), 
                            font.label = list(size = 20, face = "bold")),
                  ggarrange(temp2004, labels = "c", ncol = 1,
                            font.label = list(size = 20, face = "bold")))
panel.season2
ggsave("Figures/season_panel2.png", plot = panel.season2, width = 8, height = 6.5, units = "in")



### Microclimate Panel (temperature, cloche spectra, rH, VP/ VPD) ----
panel_otc <- ggarrange(nrow = 2,
                       ggarrange(temp, rh, common.legend = T, legend = "top",
                       labels = c("a", "b"), ncol = 2, 
                       widths = c(1, 1), 
                       font.label = list(size = 20, face = "bold")),
                       ggarrange(vp, otc,
                                 labels = c("c", "d"), 
                                 widths = c(1, 1), ncol = 2,
                                 font.label = list(size = 20, face = "bold")))
panel_otc
#ggsave("Figures/otc_panel.png", plot = panel_otc, width = 6.75, height = 6.75, units = "in")

# alternative VP plot
panel_otc2 <- ggarrange(nrow = 2,
                       ggarrange(temp, rh, common.legend = T, legend = "top",
                                 labels = c("a", "b"), ncol = 2, 
                                 widths = c(1, 1), 
                                 font.label = list(size = 20, face = "bold")),
                       ggarrange(vp2, otc2,
                                 labels = c("c", "d"), 
                                 widths = c(1, 0.8), ncol = 2,
                                 font.label = list(size = 20, face = "bold")))
panel_otc2
#ggsave("Figures/otc_panel2.png", plot = panel_otc2, width = 6.5, height = 6.5, units = "in")


# with horizontal vp 
panel_otc3 <- ggarrange(nrow = 2,
                       ggarrange(NULL, temp, NULL, rh, NULL,
                                 common.legend = T, legend = "top",
                                 labels = c("", "a", "b", ""), ncol = 5, 
                                 widths = c(0.1, 1, 0.1, 1, 0.1), 
                                 font.label = list(size = 20, face = "bold")),
                       ggarrange(vp3, otc,
                                 labels = c("c", "d"), 
                                 widths = c(1, 0.75), ncol = 2,
                                 font.label = list(size = 20, face = "bold")))
(panel_otc3 <- panel_otc3 + bgcolor("#ffffff"))
ggsave("Figures/otc_panel3.png", plot = panel_otc3, width = 6.5, height = 6.5, units = "in")
