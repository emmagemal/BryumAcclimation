# Data Analysis and Calculations for Pulse Warming Experiment 
# Emma Gemal, emmagemal@outlook.com
# University of Edinburgh, Stockholm University 

### Library ----
library(tidyverse)
library(retistruct)
library(plotrix)  # for std.error 
library(lme4)
library(lmerTest)
library(lmtest)
library(car)
library(sjPlot)
library(MuMIn)

### NP/DR/GP Data Manipulation ----
fulldata <- read.csv("Data/Pulse_Experiment/raw_np_dr_data.csv", header = TRUE)
# avgdata <- read.csv("Data/Pulse_Experiment/np_dr_averages.csv", header = TRUE)

str(fulldata)
# str(avgdata)

fulldata <- fulldata %>% 
              dplyr::select(-c(np_DW, np_Chl)) %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type),
                     sample = as.factor(sample)) %>% 
              mutate(np_SA = abs(np_SA))

# calculating GP and % max rate data 
fulldata <- fulldata %>% 
                pivot_wider(names_from = type, values_from = np_SA) %>% 
                mutate(GP = NP + DR) %>% 
                pivot_longer(cols = c("DR", "NP", "GP"), names_to = "type", 
                             values_to = "np_SA")

avgdata <- fulldata %>% 
              group_by(temp, treatment_type, type) %>% 
              summarize(avgSA = mean(np_SA),
                        se_SA = std.error(np_SA))

# subsetting average NP and DR for calculations
np_only <- avgdata %>% 
              filter(type == "NP")
dr_only <- avgdata %>% 
              filter(type == "DR")
gp_only <- avgdata %>% 
              filter(type == "GP")

np_control <- np_only %>% 
                filter(treatment_type == "control")
np_treatment <- np_only %>% 
                  filter(treatment_type == "treatment")
dr_control <- dr_only %>% 
                filter(treatment_type == "control")
dr_treatment <- dr_only %>% 
                  filter(treatment_type == "treatment")

# subsetting full NP and DR for models and calculations  
np_full <- fulldata %>% 
              filter(type == "NP")
dr_full <- fulldata %>% 
              filter(type == "DR")

np_full_control <- np_full %>% 
                      filter(treatment_type == "control")
np_full_treatment <- np_full %>% 
                        filter(treatment_type == "treatment")
### Carbon Use Efficiency ----
# calculating CUE
cuedata <- fulldata %>% 
              pivot_wider(names_from = type, values_from = np_SA) %>% 
              group_by(sample, temp, treatment_type) %>% 
              mutate(CUE = NP/(GP)*100) %>% 
              dplyr::select(!c(GP, DR, NP)) 
            
cue_sum <- cuedata %>% 
              group_by(temp, treatment_type) %>% 
              summarize(avgCUE = mean(CUE),
                        seCUE = std.error(CUE))

ggplot(cue_sum, aes(x = temp, y = avgCUE)) +
  geom_point(aes(color = treatment_type))

cue_sum2 <- cuedata %>% 
              group_by(treatment_type) %>% 
              summarize(avgCUE = mean(CUE),
                        seCUE = std.error(CUE))

# t-test
t.test(CUE ~ treatment_type, data = cuedata)  # t = -1.2055, df = 67.17, p-value = 0.2322

# t-test (averages)
t.test(avgCUE ~ treatment_type, data = cue_sum)  # t = -1.4134, df = 10.088, p-value = 0.1876

## T-test for each temperature step separately 
# 2˚C
cue2 <- cuedata %>% filter(temp == 2)
t.test(CUE ~ treatment_type, data = cue2)
  # t = -1.6484, df = 3.4419, p-value = 0.186

# 5˚C
cue5 <- cuedata %>% filter(temp == 5)
t.test(CUE ~ treatment_type, data = cue5)
  # t = -1.0324, df = 9.905, p-value = 0.3265

# 10˚C
cue10 <- cuedata %>% filter(temp == 10)
t.test(CUE ~ treatment_type, data = cue10)
  # -1.4864, df = 9.3348, p-value = 0.1701

# 15˚C
cue15 <- cuedata %>% filter(temp == 15)
t.test(CUE ~ treatment_type, data = cue15)
  # t = -1.8075, df = 9.5971, p-value = 0.1021

# 20˚C
cue20 <- cuedata %>% filter(temp == 20)
t.test(CUE ~ treatment_type, data = cue20)
  # t = -0.36116, df = 5.9937, p-value = 0.7304

# 25˚C
cue25 <- cuedata %>% filter(temp == 25)
t.test(CUE ~ treatment_type, data = cue25)
  # t = 0.038116, df = 9.8632, p-value = 0.9704

# 30˚C
cue30 <- cuedata %>% filter(temp == 30)
t.test(CUE ~ treatment_type, data = cue30)
  # t = 1.0786, df = 3.7646, p-value = 0.345


### Acclimation Ratios & Q10 ----
## Calculation of AR
ARratio <- as.data.frame(matrix(ncol = 3, nrow = 12))
colnames(ARratio) <- c("temp", "type", "ratio")

# equation = treatment at T + 5˚C / control at T
ARratio$ratio <- c((-0.6804083/-0.4568380), (-2.5739912/-0.7737948), (-4.2830330/-1.7684478),   # DR
                   (-6.4027841/-3.0684388), (-7.9797452/-4.9741951), (-10.5753313/-6.1260173),  # DR
                   (2.4157606/1.7340523), (4.6715594/1.3318238), (4.5930837/1.7676349),         # NP
                   (2.7201547/0.8418550), (0.9378886/abs(-0.7544705)), (-2.1095351/-1.9636927)) # NP
ARratio$temp <- c(3.5, 7.5, 12.5, 17.5, 22.5, 27.5, 3.5, 7.5, 12.5, 17.5, 22.5, 27.5)
                    # want to represent that it's a range/ calculated across temperatures 
ARratio$type <- c(rep("DR", times = 6), rep("NP", times = 6))
ARratio <- ARratio[-1,]  # removing outliers

# AR 10-20˚C (DR)
-6.4027841/-1.7684478   # AR treat20/cont10
    # 3.620567


# AR statistics 
summary(ARratio$ratio[ARratio$type == "DR"])
#  Min.    Median    Mean     Max. 
# 1.604    2.087    2.233    3.326 

summary(ARratio$ratio[ARratio$type == "NP"])
#  Min.    Median    Mean     Max. 
# 1.074    1.996    2.175     3.508 

std.error(ARratio$ratio[ARratio$type == "DR"])   # 0.3086133
std.error(ARratio$ratio[ARratio$type == "NP"])   # 0.4382548




### Optimum Temperature (Max NP) Statistics ----
## Control samples optimum temperature calculations 
# maximum net photosynthesis for each sample 
summary(np_full_control$np_SA[np_full_control$sample == "C1"])  # max C1 = 3.1626
summary(np_full_control$np_SA[np_full_control$sample == "C2"])  # max C2 = 0.7603
summary(np_full_control$np_SA[np_full_control$sample == "C3"])  # max C3 = 2.63050 
summary(np_full_control$np_SA[np_full_control$sample == "C4"])  # max C4 = 4.576 
summary(np_full_control$np_SA[np_full_control$sample == "C5"])  # max C5 = 1.21684 
summary(np_full_control$np_SA[np_full_control$sample == "C6"])  # max C6 = 1.3988 

# visualizing the intersection points 
hline_c <- data.frame(z = c(3.1626, 0.7603, 2.63050, 4.576, 1.21684, 1.3988), 
                      sample = factor(c("C1", "C2", "C3", "C4", "C5", "C6")))

(np_control_plot <- ggplot(np_full_control, aes(x = temp, y = np_SA)) +
                      geom_point(aes(color = sample), 
                                 size = 2, alpha = 0.85) +
                      geom_line(aes(color = sample)) +
                      facet_wrap(~sample) +
                      geom_hline(data = hline_c, aes(yintercept = z)))

## Treatment samples optimum temperature calculations
# maximum net photosynthesis for each sample 
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T1"])  # max T1 = 3.3194
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T2"])  # max T2 = 4.9480
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T3"])  # max T3 = 6.2579
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T4"])  # max T4 = 4.794
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T5"])  # max T5 = 6.9659
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T6"])  # max T6 = 9.2087

# visualizing the intersection points 
hline_t <- data.frame(z = c(3.3194, 4.9480, 6.2579, 4.794, 6.9659, 9.2087), 
                      sample = factor(c("T1", "T2", "T3", "T4", "T5", "T6")))

(np_treatment_plot <- ggplot(np_full_treatment, aes(x = temp, y = np_SA)) +
                        geom_point(aes(color = sample), 
                                   size = 2, alpha = 0.85) +
                        geom_line(aes(color = sample)) +
                        facet_wrap(~sample) +
                        geom_hline(data = hline_t, aes(yintercept = z)))

## Creating a combined dataframe for analysis 
opt_temp <- c(15, 2, 2, 10, 10, 10, 2, 15, 15, 10, 10, 10)
sample <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3", "T4", "T5", "T6")

opt_temp <- data.frame(opt_temp, sample)

opt_temp <- opt_temp %>% 
              mutate(treatment_type = case_when(grepl("C", sample) ~ "control",
                                                grepl("T", sample) ~ "treatment"))

## Testing significance of optimum temperatures 
t.test(opt_temp ~ treatment_type, data = opt_temp)   # control = 8.2, treatment = 10.3
# p = 0.467, t = -0.75638, DF = 10 (9.9376)

opt_temp_sum <- opt_temp %>% 
                  group_by(treatment_type) %>%
                  summarise(avg = mean(opt_temp),
                            se = std.error(opt_temp))


### Negative NP Statistics ----
## Determining the control intersection points
# visualizing the intersection points
(negnp_control_plot <- ggplot(np_full_control, aes(x = temp, y = np_SA)) +
                           geom_point(aes(color = sample), 
                                      size = 2, alpha = 0.85) +
                           geom_line(aes(color = sample)) +
                           facet_wrap(~sample) +
                           geom_hline(yintercept = 0))

# C1 intersection points
c1aneg <- c(25, 0.05172152)
c1bneg <- c(30, -2.88872544	)
control_y_c1neg <- c(20, 0)
control_y2_c1neg <- c(30, 0)

line.line.intersection(c1aneg, c1bneg, control_y_c1neg, control_y2_c1neg, 
                       interior.only = FALSE)               # x = 25.088˚C

# C2 intersection points
c2aneg <- c(15, 0.49103820)
c2bneg <- c(20, -2.00524056)
control_y_c2neg <- c(10, 0)
control_y2_c2neg <- c(25, 0)

line.line.intersection(c2aneg, c2bneg, control_y_c2neg, control_y2_c2neg, 
                       interior.only = FALSE)               # x = 15.984˚C

# C3 intersection points
c3aneg <- c(15, 2.42153802)
c3bneg <- c(20, -1.95855397)
control_y_c3neg <- c(15, 0)
control_y2_c3neg <- c(20, 0)

line.line.intersection(c3aneg, c3bneg, control_y_c3neg, control_y2_c3neg, 
                       interior.only = FALSE)               # x = 17.764˚C

# C4 intersection points
c4aneg <- c(20, 2.68520787)
c4bneg <- c(25, -0.29280247)
control_y_c4neg <- c(20, 0)
control_y2_c4neg <- c(30, 0)

line.line.intersection(c4aneg, c4bneg, control_y_c4neg, control_y2_c4neg, 
                       interior.only = FALSE)               # x = 24.508˚C

# C5 intersection points
c5aneg <- c(10, 1.21684173)
c5bneg <- c(15, -2.28702710)
control_y_c5neg <- c(10, 0)
control_y2_c5neg <- c(20, 0)

line.line.intersection(c5aneg, c5bneg, control_y_c5neg, control_y2_c5neg, 
                       interior.only = FALSE)               # x = 11.736˚C
# C6 intersection points
c6aneg <- c(10, 1.39880761)
c6bneg <- c(15, -1.34316394)
control_y_c6neg <- c(10, 0)
control_y2_c6neg <- c(15, 0)

line.line.intersection(c6aneg, c6bneg, control_y_c6neg, control_y2_c6neg, 
                       interior.only = FALSE)               # x = 12.551˚C

## Determining the treatment intersection points
# visualizing the intersection points
(negnp_treatment_plot <- ggplot(np_full_treatment, aes(x = temp, y = np_SA)) +
                            geom_point(aes(color = sample), 
                                       size = 2, alpha = 0.85) +
                            geom_line(aes(color = sample)) +
                            facet_wrap(~sample) +
                            geom_hline(yintercept = 0))

# T1 intersection points
t1aneg <- c(15, 2.8796464)
t1bneg <- c(20, -0.6428889)
treatment_y_t1neg <- c(15, 0)
treatment_y2_t1neg <- c(25, 0)

line.line.intersection(t1aneg, t1bneg, treatment_y_t1neg, treatment_y2_t1neg, 
                       interior.only = FALSE)               # x = 19.087˚C

# T2 intersection points
t2aneg <- c(20, 4.3679585)
t2bneg <- c(25, -0.9695074)
treatment_y_t2neg <- c(20, 0)
treatment_y2_t2neg <- c(30, 0)

line.line.intersection(t2aneg, t2bneg, treatment_y_t2neg, treatment_y2_t2neg, 
                       interior.only = FALSE)               # x = 24.092˚C

# T3 intersection points
t3aneg <- c(25, 1.4059476)
t3bneg <- c(30, -0.7325323)
treatment_y_t3neg <- c(20, 0)
treatment_y2_t3neg <- c(30, 0)

line.line.intersection(t3aneg, t3bneg, treatment_y_t3neg, treatment_y2_t3neg, 
                       interior.only = FALSE)               # x = 28.287˚C

## Creating a combined dataframe for analysis 
negNP <- c(25.088, 15.984, 17.764, 24.508, 11.736, 12.551, 19.087, 24.092, 28.287)

sample <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3")

negNP_stats <- data.frame(negNP, sample)
str(negNP_stats)

negNP_stats <- negNP_stats %>% 
                  mutate(treatment_type = case_when(grepl("C", sample) ~ "control",
                                                    grepl("T", sample) ~ "treatment"))

# including the 3 treatment samples with no negative NP
negNP2 <- c(25.088, 15.984, 17.764, 24.508, 11.736, 12.551, 19.087, 24.092, 28.287, 30, 30, 30)

sample2 <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3", "T4", "T5", "T6")

negNP_stats2 <- data.frame(negNP2, sample2)
str(negNP_stats2)

negNP_stats2 <- negNP_stats2 %>% 
                  mutate(treatment_type = case_when(grepl("C", sample) ~ "control",
                                                    grepl("T", sample) ~ "treatment"))

## Testing significance of negative NP threshold 
t.test(negNP ~ treatment_type, data = negNP_stats)   # control = 17.939, treatment = 23.822
# t = -1.658, DF = 5.0987, p = 0.1571

# with 3 extra samples
t.test(negNP2 ~ treatment_type, data = negNP_stats2)   # control = 17.939, treatment = 26.911
# t = -3.0172, df = 9.4182, p-value = 0.01383

## Determining SE of the average
negNP_sum <- negNP_stats %>% 
                group_by(treatment_type) %>%
                mutate(se = std.error(negNP)) %>% 
                summarise(avg = mean(negNP),
                          se = mean(se),
                          sd = sd(negNP))

negNP_sum2 <- negNP_stats2 %>% 
                group_by(treatment_type) %>%
                mutate(se = std.error(negNP2)) %>% 
                summarise(avg = mean(negNP2),
                          se = std.error(negNP2))


### Maximum NP Statistics ----
# maximum net photosynthesis for each control sample 
summary(np_full_control$np_SA[np_full_control$sample == "C1"])  # max C1 = 3.1626
summary(np_full_control$np_SA[np_full_control$sample == "C2"])  # max C2 = 0.76027
summary(np_full_control$np_SA[np_full_control$sample == "C3"])  # max C3 = 2.63050
summary(np_full_control$np_SA[np_full_control$sample == "C4"])  # max C4 = 4.5760
summary(np_full_control$np_SA[np_full_control$sample == "C5"])  # max C5 = 1.2168
summary(np_full_control$np_SA[np_full_control$sample == "C6"])  # max C6 = 1.3988

# maximum net photosynthesis for each treatment sample 
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T1"])  # max T1 = 3.3194
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T2"])  # max T2 = 4.9480
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T3"])  # max T3 = 6.2579
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T4"])  # max T4 = 4.794
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T5"])  # max T5 = 6.9659
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T6"])  # max T6 = 9.2087

maxNP <- c(3.1626, 0.76027, 2.63050, 4.5760, 1.2168, 1.3988, 
           3.3194, 4.9480, 6.2579, 4.794, 6.9659, 9.2087)
sample <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3", "T4", "T5", "T6")

maxNP_stats <- data.frame(maxNP, sample)
str(maxNP_stats)

maxNP_stats <- maxNP_stats %>% 
                  mutate(treatment_type = case_when(grepl("C", sample) ~ "control",
                                                    grepl("T", sample) ~ "treatment"))

## Testing significance of maximum NP rates 
t.test(maxNP ~ treatment_type, data = maxNP_stats)

## Determining SE of the average
maxNP_sum <- maxNP_stats %>% 
                group_by(treatment_type) %>%
                mutate(se = std.error(maxNP)) %>% 
                summarise(avg = mean(maxNP),
                          se = mean(se),
                          sd = sd(maxNP))



### DR and GP Rate Comparisons ----
## DR
# difference in DR between control and treatment 
lmdr <- lm(avgSA ~ temp+treatment_type, data = dr_only)
summary(lmdr)   # treatment significantly greater DR than control (estimate = -1.025, p = 0.0116)

lmdr2 <- lm(avgSA ~ temp*treatment_type, data = dr_only)
summary(lmdr2)

ggplot(dr_only, aes(y = avgSA, x = temp)) +
  geom_point(aes(color = treatment_type)) +
  geom_smooth(method = "lm", aes(color = treatment_type))

# difference in % 
dr_perc <- dr_only %>% 
              dplyr::select(-se_SA) %>% 
              pivot_wider(names_from = treatment_type, values_from = avgSA) %>% 
              mutate(perc_diff = ((treatment-control)/control)*100)

mean(dr_perc$perc_diff)   # 28.62236
std.error(dr_perc$perc_diff)   # 7.776387

t.test(avgSA ~ treatment_type, paired = TRUE, data = dr_only)  # t = 3.6304, df = 6, 
                                                               # p-value = 0.01096 (significant)


## GP
# difference in DR between control and treatment 
lmgp <- lm(avgSA ~ temp + treatment_type, data = gp_only)
summary(lmgp)   # treatment significantly higher than control (estimate = 3.48, p = 0.00016)

t.test(avgSA ~ treatment_type, paired = TRUE, data = gp_only)  # t = -5.996, df = 6, 
                                                               # p-value = 0.00097 (significant)

# difference in % 
gp_perc <- gp_only %>% 
              dplyr::select(-se_SA) %>% 
              pivot_wider(names_from = treatment_type, values_from = avgSA) %>% 
              mutate(perc_diff = ((treatment-control)/control)*100)

mean(gp_perc$perc_diff)   # 94.39904
std.error(gp_perc$perc_diff)   # 12.41984


### NP Statistics ----
## Checking data assumptions  
# normally distributed response variable 
(hist <- ggplot(np_full, aes(x = np_SA)) +
           geom_histogram(color = "black") +
           theme_classic() +
           scale_y_continuous(expand = c(0,0)))
shapiro.test(np_full$np_DW)   # is normally distributed, p > 0.05

# homogeneity of regression slopes 
anova(lm(np_SA ~ temp*treatment_type, data = np_full))  # interaction not significant

## Models
lm_treat_temp <- lm(np_SA ~ temp + treatment_type, data = np_full)
mixed_null_np <- lmer(np_SA ~ 1 + (1|sample), data = np_full, REML = F)
mixed_notemp <- lmer(np_SA ~ treatment_type + (1|sample), data = np_full, REML = F)
mixed_sample_np <- lmer(np_SA ~ temp + treatment_type + (1|sample), data = np_full, REML = F)
# sample as a random effect because sample needs to be controlled for, but I am not
# interested in the direct relationship of it with np_DW 

# comparing the models 
AIC(mixed_null_np, mixed_notemp, mixed_sample_np)  # mixed: lowest AIC = mixed_sample_np

anova(mixed_null_np, mixed_sample_np)  # mixed_sample_np is better than the null model
anova(mixed_sample_np, lm_treat_temp)  # accounting for sample variation is best (mixed_sample_np)

# checking model assumptions 
plot(mixed_sample_np)

par(mfrow = c(1,2))
qqnorm(ranef(mixed_sample_np)$sample[, 1], main = "Random effects of sample")
qqnorm(resid(mixed_sample_np), main = "Residuals")

# model outputs
confint(mixed_sample_np)  # temp: -0.214 to -0.116
# treatment_type: 0.861 to 4.616

summary(mixed_sample_np)  # sample explains quite a bit of the residual variance 
# temp: -0.165 ± 0.0245 (within the confidence intervals)
# treatment_type: 2.728 ± 0.881 (within the confidence intervals)
# treatment_type and temp significantly effect NP
# they're significantly different from 0 according to the CI's
1.781/(1.781+3.220)    # sample explains 35.6% of the residual variation 

anova(mixed_sample_np)   # temp: F = 45.23, p = 7.465e-9, DF = 59.78
# treatment_type: F = 9.59, p = 0.0095, DF = 11.739

### DR Statistics ----
## Checking data assumptions  
# normally distributed response variable 
(hist <- ggplot(dr_full, aes(x = np_SA)) +
           geom_histogram(color = "black") +
           theme_classic() +
           scale_y_continuous(expand = c(0,0)))
shapiro.test(dr_full$np_DW)   # is NOT normally distributed, p << 0.05

# homogeneity of regression slopes 
anova(lm(np_SA ~ temp*treatment_type, data = dr_full))  # interaction not significant

## Models
lm_treat_temp_dr <- lm(np_SA ~ temp + treatment_type, data = dr_full)
mixed_null_dr <- lmer(np_SA ~ 1 + (1|sample), data = dr_full, REML = F)
mixed_notemp_dr <- lmer(np_SA ~ treatment_type + (1|sample), data = dr_full, REML = F)
mixed_sample_dr <- lmer(np_SA ~ temp + treatment_type + (1|sample), data = dr_full, REML = F)

# comparing the models 
AIC(mixed_null_dr, mixed_notemp_dr, mixed_sample_dr)  # mixed_sample_dr = lowest AIC

anova(mixed_sample_dr, mixed_null_dr)  # mixed_sample_dr is better than the null model 
anova(mixed_sample_dr, lm_treat_temp_dr)  # better to include the random effect 

# checking model assumptions 
plot(mixed_sample_dr)   # potential violation of linearity 

par(mfrow = c(1,2))
qqnorm(ranef(mixed_sample_dr)$sample[, 1], main = "Random effects of sample")
qqnorm(resid(mixed_sample_dr), main = "Residuals")

# model outputs 
confint(mixed_sample_dr)  # temp:  -0.360 to -0.293 
# treatment_type: -2.831 to 0.637 

summary(mixed_sample_dr)  # sample explains a lot of the residual variance 
# temp: -0.327 ± 0.0170
# treatment_type: -1.094 ± 0.813
# temp is significant, treatment_type is not 
1.720/(1.720+1.546)  # sample explains 52.7% of the residual variation

anova(mixed_sample_dr)   # temp: F = 369.69, p = <2e-16, DF = 59.69
# treatment_type: F = 1.811, p = 0.204, DF = 11.68
# treatment_type is not significant, temp is significant 


### NP Non-linear Curve Fits ----
# control
poly1c <- lm(avgSA ~ poly(temp, 2), data = np_control)

poly1c.pred <- data.frame(avgSA = rep(NA, 7), temp = rep(NA, 7))
poly1c.pred$avgSA <- predict(poly1c)
poly1c.pred$temp <- c(2,5,10,15,20,25,30)

ggplot(np_control, aes(x = temp, y = avgSA)) +
  geom_point() +
  geom_line(data = poly1c.pred, aes(x = temp, y = avgSA))

summary(poly1c)
# adjusted R2 = 0.958
# p = 0.0003 

# treatment
poly1t <- lm(avgSA ~ poly(temp, 2), data = np_treatment)

poly1t.pred <- data.frame(avgSA = rep(NA, 7), temp = rep(NA, 7))
poly1t.pred$avgSA <- predict(poly1t)
poly1t.pred$temp <- c(2,5,10,15,20,25,30)

ggplot(np_treatment, aes(x = temp, y = avgSA)) +
  geom_point() +
  geom_line(data = poly1t.pred, aes(x = temp, y = avgSA))

summary(poly1t)
# adjusted R2 = 0.876
# p = 0.0067

### DR Non-linear Curve Fits ----
# control
poly2c <- lm(avgSA ~ poly(temp, 2), data = dr_control)

poly2c.pred <- data.frame(avgSA = rep(NA, 7), temp = rep(NA, 7))
poly2c.pred$avgSA <- predict(poly2c)
poly2c.pred$temp <- c(2,5,10,15,20,25,30)

ggplot(dr_control, aes(x = temp, y = avgSA)) +
  geom_point() +
  geom_line(data = poly2c.pred, aes(x = temp, y = avgSA))

summary(poly2c)
# adjusted R2 = 0.992
# p = 1.01e-5

# treatment
poly2t <- lm(avgSA ~ poly(temp, 2), data = dr_treatment)

poly2t.pred <- data.frame(avgSA = rep(NA, 7), temp = rep(NA, 7))
poly2t.pred$avgSA <- predict(poly2t)
poly2t.pred$temp <- c(2,5,10,15,20,25,30)

ggplot(dr_treatment, aes(x = temp, y = avgSA)) +
  geom_point() +
  geom_line(data = poly2t.pred, aes(x = temp, y = avgSA))

summary(poly2t)
# adjusted R2 = 0.992
# p = 1.02e-5


### GP Non-linear Curve Fits ----

## ADD IF NECESSARY ##

## Light Response Curve LSP Statistics ----
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

# create Asat, Phi, Rd, theta, resid_SSs, LCP, Q_sat_75 and Q_sat_85 for each sample 
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
# equation used is:
  # (Phi * PAR + Asat - sqrt((Phi * PAR + Asat)^2 - 4 * Phi * theta * Asat * PAR ))/(2*theta) - Rd)

# C1 
curve.c1 <- function(PARi){
  (0.07287260 * PARi + 8.675147 - 
     sqrt((0.07287260 * PARi + 8.675147)^2 - 4 *
            0.07287260 * -1.2748040 * PARi *
            8.675147)
  ) / (2*-1.2748040) - 2.873080
}

par.c1 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.c1(PARi = seq(0,1550, by = 50)),
                     sample = "C1")
# C2
curve.c2 <- function(PARi){
  (0.08371293 * PARi + 4.422790 - 
     sqrt((0.08371293 * PARi + 4.422790)^2 - 4 *
            0.08371293 * -1.3723880 * PARi *
            4.422790)
  ) / (2*-1.3723880) - 2.562400
}

par.c2 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.c2(PARi = seq(0,1550, by = 50)),
                     sample = "C2")
# C3
curve.c3 <- function(PARi){
  (0.07674525 * PARi + 13.465161 - 
     sqrt((0.07674525* PARi + 13.465161)^2 - 4 *
            0.07674525 * -1.1614485 * PARi *
            13.465161)
  ) / (2*-1.1614485) - 3.212828
}

par.c3 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.c3(PARi = seq(0,1550, by = 50)),
                     sample = "C3")
# T1
curve.t1 <- function(PARi){
  (0.11458331 * PARi + 13.982050 - 
     sqrt((0.11458331 * PARi + 13.982050)^2 - 4 *
            0.11458331 * 0.7072400 * PARi *
            13.982050)
  ) / (2*0.7072400) - 2.823915
}

par.t1 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.t1(PARi = seq(0,1550, by = 50)),
                     sample = "T1")
# T2
curve.t2 <- function(PARi){
  (0.08046531 * PARi + 10.091202 - 
     sqrt((0.08046531 * PARi + 10.091202)^2 - 4 *
            0.08046531 * 0.9372664 * PARi *
            10.091202)
  ) / (2*0.9372664) - 2.030062
}

par.t2 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.t2(PARi = seq(0,1550, by = 50)),
                     sample = "T2")
# T3
curve.t3 <- function(PARi){
  (0.16754053 * PARi + 18.954037 - 
     sqrt((0.16754053 * PARi + 18.954037)^2 - 4 *
            0.16754053 * 0.1385859 * PARi *
            18.954037)
  ) / (2*0.1385859) - 3.666333
}

par.t3 <- data.frame(Lcuv = seq(0,1550, by = 50),
                     curve = curve.t3(PARi = seq(0,1550, by = 50)), 
                     sample = "T3")

par.all <- rbind(par.c1, par.c2, par.c3, par.t1, par.t2, par.t3)
par.all <- par.all %>% 
              mutate(treatment_type = ifelse(grepl("C", sample), "Control", "Treatment"))

light_sum<- light %>% 
            group_by(sample) %>% 
            summarize(maxNP = max(NP_SA)) %>% 
            mutate(LSPx90 = 0.9*maxNP)

par_sum<- par.all %>% 
            group_by(sample) %>% 
            summarize(maxNP = max(curve)) %>% 
            mutate(LSPx90 = 0.9*maxNP)


## Determining the intersection points for 90% LSP
# visualizing the 90% light saturation points for all samples
(light_plot <- ggplot(par.all, aes(x = Lcuv, y = curve)) +
                  geom_line(aes(color = sample, linetype = treatment_type),
                            size = 0.8) +
                  geom_hline(data = light_sum, aes(yintercept = LSPx90, color = sample)) +
                  facet_wrap(~sample))

# C1 intersection point
c1_1 <- c(1000, 4.115659)
c1_2 <- c(1050, 4.1752662)
control_y_c1 <- c(500, 3.960)
control_y2_c1 <- c(1200, 3.960)

line.line.intersection(c1_1, c1_2, control_y_c1, control_y2_c1, 
                       interior.only = FALSE)               # x = 1019.08 µE m^-2 s^-1

# C2 intersection point
c2_1 <- c(500, 1.15)
c2_2 <- c(1000, 1.42)
control_y_c2 <- c(500, 1.278)
control_y2_c2 <- c(1000, 1.278)

line.line.intersection(c2_1, c2_2, control_y_c2, control_y2_c2, 
                       interior.only = FALSE)               # x = 737.04 µE m^-2 s^-1

# C3 intersection point
c3_1 <- c(500, 5.11)
c3_2 <- c(1000, 7.25)
control_y_c3 <- c(500, 6.660)
control_y2_c3 <- c(1000, 6.660)

line.line.intersection(c3_1, c3_2, control_y_c3, control_y2_c3, 
                       interior.only = FALSE)               # x = 862.15 µE m^-2 s^-1

# T1 intersection point
t1_1 <- c(200, 7.98)
t1_2 <- c(500, 11.16)
treatment_y_t1 <- c(200, 10.044)
treatment_y2_t1 <- c(500, 10.044)

line.line.intersection(t1_1, t1_2, treatment_y_t1, treatment_y2_t1, 
                       interior.only = FALSE)               # x = 394.72 µE m^-2 s^-1

# T2 intersection point
t2_1 <- c(200, 7.14)
t2_2 <- c(500, 9.99)
treatment_y_t2 <- c(200, 8.991)
treatment_y2_t2 <- c(500, 8.991)

line.line.intersection(t2_1, t2_2, treatment_y_t2, treatment_y2_t2, 
                       interior.only = FALSE)               # x = 394.84 µE m^-2 s^-1

# T3 intersection point
t3_1 <- c(200, 8.82)
t3_2 <- c(500, 12.33)
treatment_y_t3 <- c(200, 12.303)
treatment_y2_t3 <- c(500, 12.303)

line.line.intersection(t3_1, t3_2, treatment_y_t3, treatment_y2_t3, 
                       interior.only = FALSE)               # x = 497.69 µE m^-2 s^-1

light_lsp <- light_sum
light_lsp$LSPcuv <- c(845.83, 737.04, 862.15, 394.72, 394.84, 497.69)
light_lsp$treatment_type <- c("control", "control", "control", 
                              "treatment", "treatment", "treatment")              

## T-Test for LSP
t.test(LSPcuv ~ treatment_type, data = light_lsp)  # p = 0.0019 (significant difference!)
# t = 7.402, df = 3.929
# mean control LSP = 815.0067
# mean treatment LSP = 429.0833

## Determining SE of the average
str(light_lsp)
light_sum_lsp <- light_lsp %>% 
                    group_by(treatment_type) %>%
                    mutate(se = std.error(LSPcuv)) %>% 
                    summarise(avg = mean(LSPcuv),
                              se = mean(se),
                              sd = sd(LSPcuv))


### Light Response Curve LCP Statistics ----
# visualizing the light compensation point for all samples
(lightcomp_plot_facet <- ggplot(light, aes(x = Lcuv, y = NP_SA)) +
                           geom_point(aes(color = treatment_type)) +
                           geom_line(aes(color = treatment_type)) +
                           facet_wrap(~sample) +
                           geom_hline(yintercept = 0))

## Determining the intersection points for LCP
# C1 intersection point
c1_1b <- c(50, -0.73)
c1_2b <- c(100, 0.25)
control_y_c1b <- c(50, 0)
control_y2_c1b <- c(200, 0)

line.line.intersection(c1_1b, c1_2b, control_y_c1b, control_y2_c1b, 
                       interior.only = FALSE)               # x = 87.25 µE m^-2 s^-1

# C2 intersection point
c2_1b <- c(100, -0.32)
c2_2b <- c(200, 0.37)
control_y_c2b <- c(100, 0)
control_y2_c2b <- c(200, 0)

line.line.intersection(c2_1b, c2_2b, control_y_c2b, control_y2_c2b, 
                       interior.only = FALSE)               # x = 146.38 µE m^-2 s^-1

# C3 intersection point
c3_1b <- c(50, -0.60)
c3_2b <- c(100, 0.82)
control_y_c3b <- c(50, 0)
control_y2_c3b <- c(100, 0)

line.line.intersection(c3_1b, c3_2b, control_y_c3b, control_y2_c3b, 
                       interior.only = FALSE)               # x = 71.13 µE m^-2 s^-1

# T1 intersection point
t1_1b <- c(25, -0.04)
t1_2b <- c(50, 2.26)
treatment_y_t1b <- c(10, 0)
treatment_y2_t1b <- c(50, 0)

line.line.intersection(t1_1b, t1_2b, treatment_y_t1b, treatment_y2_t1b, 
                       interior.only = FALSE)               # x = 25.43 µE m^-2 s^-1

# T2 intersection point
t2_1b <- c(12, -0.96)
t2_2b <- c(25, 0.23)
treatment_y_t2b <- c(10, 0)
treatment_y2_t2b <- c(50, 0)

line.line.intersection(t2_1b, t2_2b, treatment_y_t2b, treatment_y2_t2b, 
                       interior.only = FALSE)               # x = 22.49 µE m^-2 s^-1

# T3 intersection point
t3_1b <- c(25, -0.32)
t3_2b <- c(50, 2.45)
treatment_y_t3b <- c(20, 0)
treatment_y2_t3b <- c(60, 0)

line.line.intersection(t3_1b, t3_2b, treatment_y_t3b, treatment_y2_t3b, 
                       interior.only = FALSE)               # x = 27.89 µE m^-2 s^-1

light_lcp <- light_sum
light_lcp$LCPcuv <- c(87.25, 146.38, 71.13, 25.43, 22.49, 27.89)
light_lcp$treatment_type <- c("control", "control", "control", 
                              "treatment", "treatment", "treatment")   

## T-Test for LCP
t.test(LCPcuv ~ treatment_type, data = light_lcp)  # p = 0.079 (NOT significantly different)
# t = 3.329, DF = 2.0186 
# mean control LCP = 101.58
# mean treatment LCP = 25.27

## Determining SE of the average
light_sum_lcp <- light_lcp %>% 
                    group_by(treatment_type) %>%
                    mutate(se = std.error(LCPcuv)) %>% 
                    summarise(avg = mean(LCPcuv),
                              se = mean(se),
                              sd = sd(LCPcuv))


### Optimal Water Content Ranges ----
water <- read.csv("Data/Pulse_Experiment/water_content_full.csv")

water_long <- water %>% 
  pivot_longer(cols = c(1:2),
               names_to = "type",
               names_prefix = "CO2_",
               values_to = "CO2") %>% 
  filter(type == "NP")

summary(water_long$CO2[water_long$treatment_type == "control"])  # max control = 6.700
summary(water_long$CO2[water_long$treatment_type == "treatment"])  # max treatment = 2.400

6.7*0.9   # control = 6.03
2.4*0.9   # treatment = 2.16

(water_curves <- ggplot(water_long, aes(x = weight, y = CO2)) +
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
    scale_color_manual(values = "#7A292A",
                       name = "Process") +
    scale_shape_discrete(name = "Process"))

hline_water <- data.frame(z = c(6.7, 2.4), treatment_type = factor(c("control", "treatment")))

(water_curves <- water_curves + geom_hline(data = hline_water, aes(yintercept = z)))

## Determining the control intersection points
c_a <- c(12.8084, 1.6)
c_b <- c(10.8021, 6.7)
c_c <- c(10.5568, 4.7)
control_y <- c(10, 6.03)
control_y2 <- c(14, 6.03)

line.line.intersection(c_a, c_b, control_y, control_y2, interior.only = FALSE)    # x = 11.066 g
line.line.intersection(c_b, c_c, control_y, control_y2, interior.only = FALSE)    # x = 10.719 g

## Determining the treatment intersection points
t_a <- c(15.6690, -6.3)
t_b <- c(13.1344, 2.4)
t_c <- c(12.7396, 0.1)
treatment_y <- c(12, 2.16)
treatment_y2 <- c(16, 2.16)

line.line.intersection(t_a, t_b, treatment_y, treatment_y2, 
                       interior.only = FALSE)                 # x = 13.204 g
line.line.intersection(t_b, t_c, treatment_y, treatment_y2, 
                       interior.only = FALSE)                 # x = 13.093 g


### Chlorophyll Content ----
chlr <- read.csv("Data/Pulse_Experiment/chlorophyll.csv")

str(chlr)

t.test(chl_SA ~ type, data = chlr)  # p = 0.047 (is significantly different)
# t = -2.2996, DF = 8.7
# control = 736475.4, treatment = 1250650.3

