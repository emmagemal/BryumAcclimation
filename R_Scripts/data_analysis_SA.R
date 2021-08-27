# Data Analysis and Calculations USING SURFACE AREA (not DW)
# Emma Gemal, s1758915@sms.ed.ac.uk
# University of Edinburgh 

### Library ----
library(tidyverse)
library(retistruct)
library(plotrix)
library(lme4)
library(lmerTest)
library(lmtest)
library(car)
library(sjPlot)
library(MuMIn)

### NP/DR Data Manipulation ----
fulldata <- read.csv("Data/raw_np_dr_data.csv", header = TRUE)
avgdata <- read.csv("Data/np_dr_averages.csv", header = TRUE)

str(fulldata)
str(avgdata)

fulldata <- fulldata %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type),
                     sample = as.factor(sample))

avgdata <- avgdata %>% 
              mutate(treatment_type = as.factor(treatment_type),
                     type = as.factor(type))

# subsetting average NP and DR for calculations
np_only <- avgdata %>% 
              filter(type == "NP")
dr_only <- avgdata %>% 
              filter(type == "DR")

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

### Calculating Average Optimum Temperature Ranges ----
summary(np_control)  # max control NP = 1.7676
summary(np_treatment)  # max treatment NP = 4.672

summary(dr_control)  # max control DR (min. number) = -8.8291
summary(dr_treatment)  # max treatment DR = -10.5753

# calculating 90% of the maximum net photosynthesis 
0.9*1.7676  # control = 1.59084
0.9*4.672  # treatment = 4.2048

# visualizing the intersection points 
(np_plot <- ggplot(np_only, aes(x = temp, y = avgSA)) +
              geom_point(aes(color = treatment_type, shape = treatment_type), 
                         size = 2.5, alpha = 0.85) +
              geom_line(aes(color = treatment_type)) +
              facet_wrap(~treatment_type))

hline <- data.frame(z = c(1.59084, 4.2048), treatment_type = factor(c("control", "treatment")))

(np_plot <- np_plot + geom_hline(data = hline, aes(yintercept = z)))

## Determining the intersection points
# control intersection points
c1 <- c(5, 1.3318238)
c2 <- c(10, 1.7676349)
c3 <- c(10, 1.7676349)
c4 <- c(15, 0.8418550)
control_y <- c(0, 1.59084)
control_y2 <- c(20, 1.59084)

line.line.intersection(c1, c2, control_y, control_y2, 
                       interior.only = FALSE)               # x = 7.972˚C
line.line.intersection(c3, c4, control_y, control_y2, 
                       interior.only = FALSE)               # x = 10.955˚C

# treatment intersection points 
t1 <- c(5, 2.4157606)
t2 <- c(10, 4.6715594)
t3 <- c(15, 4.5930837)
t4 <- c(20, 2.7201547)
treatment_y <- c(0, 4.2048)
treatment_y2 <- c(30, 4.2048)

line.line.intersection(t1, t2, treatment_y, treatment_y2, 
                       interior.only = FALSE)               # x = 8.965˚C
line.line.intersection(t3, t4, treatment_y, treatment_y2, 
                       interior.only = FALSE)               # x = 16.037˚C

### Optimum Temperature Range Statistics ----
## Control samples optimum temperature calculations 
# maximum net photosynthesis for each sample 
summary(np_full_control$np_SA[np_full_control$sample == "C1"])  # max C1 = 3.1626
summary(np_full_control$np_SA[np_full_control$sample == "C2"])  # max C2 = 0.7603
summary(np_full_control$np_SA[np_full_control$sample == "C3"])  # max C3 = 2.63050 
summary(np_full_control$np_SA[np_full_control$sample == "C4"])  # max C4 = 4.576 
summary(np_full_control$np_SA[np_full_control$sample == "C5"])  # max C5 = 1.21684 
summary(np_full_control$np_SA[np_full_control$sample == "C6"])  # max C6 = 1.3988 

# calculating 90% of the maximum net photosynthesis 
0.9*3.1626  # C1 = 2.84634
0.9*0.7603  # C2 = 0.68427
0.9*2.63050  # C3 = 2.36745
0.9*4.576  # C4 = 4.1184
0.9*1.21684  # C5 = 1.095156
0.9*1.3988  # C6 = 1.25892

# visualizing the intersection points 
hline_c <- data.frame(z = c(2.84634, 0.68427, 2.36745, 4.1184, 1.095156, 1.25892), 
                      sample = factor(c("C1", "C2", "C3", "C4", "C5", "C6")))

(np_control_plot <- ggplot(np_full_control, aes(x = temp, y = np_SA)) +
                      geom_point(aes(color = sample), 
                                 size = 2, alpha = 0.85) +
                      geom_line(aes(color = sample)) +
                      facet_wrap(~sample) +
                      geom_hline(data = hline_c, aes(yintercept = z)))

## Determining the control intersection points
# C1 intersection points
c1a <- c(10, 2.06268241)
c1b <- c(15, 3.16261293)
c1c <- c(20, 1.28562409)
control_y_c1 <- c(5, 2.84634)
control_y2_c1 <- c(25, 2.84634)

line.line.intersection(c1a, c1b, control_y_c1, control_y2_c1, 
                       interior.only = FALSE)               # x = 13.562˚C
line.line.intersection(c1b, c1c, control_y_c1, control_y2_c1, 
                       interior.only = FALSE)               # x = 15.843˚C

# C2 intersection points
c2a <- c(2, 0.76027316)
c2b <- c(5, 0.51264099)
control_y_c2 <- c(0, 0.68427)
control_y2_c2 <- c(10, 0.68427)

line.line.intersection(c2a, c2b, control_y_c2, control_y2_c2, 
                       interior.only = FALSE)               # x = 2.921˚C

# C3 intersection points 

# (new intersection, barely intersects but it does)
c3c <- c(10, 1.28028326)
c3d <- c(15, 2.42153802)
c3e <- c(20, -1.95855397)
control_y_c3b <- c(10, 2.36745)
control_y2_c3b <- c(20, 2.36745)

line.line.intersection(c3c, c3d, control_y_c3b, control_y2_c3b, 
                       interior.only = FALSE)               # x = 14.763˚C
line.line.intersection(c3d, c3e, control_y_c3b, control_y2_c3b, 
                       interior.only = FALSE)               # x = 15.062˚C

# C4 intersection points
c4a <- c(5, 2.74946746)
c4b <- c(10, 4.57602446)
c4c <- c(15, 2.60613202)
control_y_c4 <- c(2, 4.1184)
control_y2_c4 <- c(20, 4.1184)

line.line.intersection(c4a, c4b, control_y_c4, control_y2_c4, 
                       interior.only = FALSE)               # x = 8.747˚C
line.line.intersection(c4b, c4c, control_y_c4, control_y2_c4, 
                       interior.only = FALSE)               # x = 11.162˚C

# C5 intersection points
c5a <- c(5, 0.20439546)
c5b <- c(10, 1.21684173)
c5c <- c(15, -2.28702710)
control_y_c5 <- c(2, 1.095156)
control_y2_c5 <- c(20, 1.095156)

line.line.intersection(c5a, c5b, control_y_c5, control_y2_c5, 
                       interior.only = FALSE)               # x = 9.399˚C
line.line.intersection(c5b, c5c, control_y_c5, control_y2_c5, 
                       interior.only = FALSE)               # x = 10.174˚C

# C6 intersection points
c6a <- c(5, 0.94727218)
c6b <- c(10, 1.39880761)
c6c <- c(15, -1.34316394)
control_y_c6 <- c(2, 1.25892)
control_y2_c6 <- c(20, 1.25892)

line.line.intersection(c6a, c6b, control_y_c6, control_y2_c6, 
                       interior.only = FALSE)               # x = 8.451˚C
line.line.intersection(c6b, c6c, control_y_c6, control_y2_c6, 
                       interior.only = FALSE)               # x = 10.255˚C

## Treatment samples optimum temperature calculations
# maximum net photosynthesis for each sample 
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T1"])  # max T1 = 3.3194
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T2"])  # max T2 = 4.9480
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T3"])  # max T3 = 6.2579
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T4"])  # max T4 = 4.794
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T5"])  # max T5 = 6.9659
summary(np_full_treatment$np_SA[np_full_treatment$sample == "T6"])  # max T6 = 9.2087

# calculating 90% of the maximum net photosynthesis 
0.9*3.3194  # T1 = 2.98746
0.9*4.9480  # T2 = 4.4532
0.9*6.2579  # T3 = 5.63211
0.9*4.794  # T4 = 4.3146
0.9*6.9659  # T5 = 6.26931
0.9*9.2087  # T6 = 8.28783

# visualizing the intersection points 
hline_t <- data.frame(z = c(2.98746, 4.4532, 5.63211, 4.3146, 6.26931, 8.28783), 
                      sample = factor(c("T1", "T2", "T3", "T4", "T5", "T6")))

(np_treatment_plot <- ggplot(np_full_treatment, aes(x = temp, y = np_SA)) +
                        geom_point(aes(color = sample), 
                                   size = 2, alpha = 0.85) +
                        geom_line(aes(color = sample)) +
                        facet_wrap(~sample) +
                        geom_hline(data = hline_t, aes(yintercept = z)))

## Determining the treatment intersection points
# T1 intersection points
t1a <- c(2, 3.3193799)
t1b <- c(5, 2.3092860)
treatment_y_t1 <- c(0, 2.98746)
treatment_y2_t1 <- c(10, 2.98746)

line.line.intersection(t1a, t1b, treatment_y_t1, treatment_y2_t1, 
                       interior.only = FALSE)               # x = 2.986˚C

# T2 intersection points
t2a <- c(10, 1.2574992)
t2b <- c(15, 4.9480402)
t2c <- c(20, 4.3679585)
treatment_y_t2 <- c(0, 4.4532)
treatment_y2_t2 <- c(25, 4.4532)

line.line.intersection(t2a, t2b, treatment_y_t2, treatment_y2_t2, 
                       interior.only = FALSE)               # x = 14.330˚C
line.line.intersection(t2b, t2c, treatment_y_t2, treatment_y2_t2, 
                       interior.only = FALSE)               # x = 19.265˚C

# T3 intersection points
t3a <- c(10, 4.0639881)
t3b <- c(15, 6.2579372)
t3c <- c(20, 1.8971447)
treatment_y_t3 <- c(5, 5.63211)
treatment_y2_t3 <- c(25, 5.63211)

line.line.intersection(t3a, t3b, treatment_y_t3, treatment_y2_t3, 
                       interior.only = FALSE)               # x = 13.574˚C
line.line.intersection(t3b, t3c, treatment_y_t3, treatment_y2_t3, 
                       interior.only = FALSE)               # x = 15.718˚C

# T4 intersection points
t4a <- c(5, 3.3360742)
t4b <- c(10, 4.7935230)
t4c <- c(15, 3.2997424)
treatment_y_t4 <- c(2, 4.3146)
treatment_y2_t4 <- c(20, 4.3146)

line.line.intersection(t4a, t4b, treatment_y_t4, treatment_y2_t4, 
                       interior.only = FALSE)               # x = 8.357˚C
line.line.intersection(t4b, t4c, treatment_y_t4, treatment_y2_t4, 
                       interior.only = FALSE)               # x = 11.603˚C

# T5 intersection points
t5a <- c(5, 0.4346008)
t5b <- c(10, 6.9658958)
t5c <- c(15, 3.1370500)
treatment_y_t5 <- c(2, 6.26931)
treatment_y2_t5 <- c(20, 6.26931)

line.line.intersection(t5a, t5b, treatment_y_t5, treatment_y2_t5, 
                       interior.only = FALSE)               # x = 9.467˚C
line.line.intersection(t5b, t5c, treatment_y_t5, treatment_y2_t5, 
                       interior.only = FALSE)               # x = 10.910˚C

# T6 intersection points
t6a <- c(5, 0.5804601)
t6b <- c(10, 9.2087301)
t6c <- c(15, 7.0360861)
treatment_y_t6 <- c(2, 8.28783)
treatment_y2_t6 <- c(20, 8.28783)

line.line.intersection(t6a, t6b, treatment_y_t6, treatment_y2_t6, 
                       interior.only = FALSE)               # x = 9.466˚C
line.line.intersection(t6b, t6c, treatment_y_t6, treatment_y2_t6, 
                       interior.only = FALSE)               # x = 12.119˚C


## Creating a combined dataframe for analysis 
opt_temp_max <- c(15.843, 2.921, 15.062, 11.162, 10.174, 10.255,
                  2.986, 19.265, 15.718, 11.603, 10.910, 12.119)
opt_temp_min <- c(13.562, 2.921, 14.763, 8.747, 9.399, 8.451, 
                  2.986, 14.330, 13.574, 8.357, 9.467, 9.466)
sample_max <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3", "T4", "T5", "T6")
sample_min <- c("C1", "C2", "C3", "C4", "C5", "C6", "T1", "T2", "T3", "T4", "T5", "T6")

opt_stats <- data.frame(opt_temp_max, opt_temp_min, sample_max, sample_min)
str(opt_stats)

opt_stats <- opt_stats %>% 
                mutate(treatment_type = case_when(grepl("C", sample_max) ~ "control",
                                                  grepl("T", sample_max) ~ "treatment"))

## Testing significance of optimum temperature ranges 
t.test(opt_temp_max ~ treatment_type, data = opt_stats)   # control = 10.9, treatment = 12.1
# p = 0.689, t = -0.412, DF = 10 (9.736)

t.test(opt_temp_min ~ treatment_type, data = opt_stats)   # control = 9.6, treatment = 9.7
# p = 0.98, t = -0.023, DF = 10 (9.989)

## Determining SE of the averages
opt_sum <- opt_stats %>% 
              group_by(treatment_type) %>%
              mutate(se_max = std.error(opt_temp_max),
                     se_min = std.error(opt_temp_min)) %>% 
              summarise(avg_max = mean(opt_temp_max),
                        avg_min = mean(opt_temp_min),
                        se_max = mean(se_max),
                        se_min = mean(se_min),
                        sd_max = sd(opt_temp_max),
                        sd_min = sd(opt_temp_min))


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

## Testing significance of optimum temperature ranges 
t.test(negNP ~ treatment_type, data = negNP_stats)   # control = 17.939, treatment = 23.822
# t = -1.658, DF = 5.0987, p = 0.1571

## Determining SE of the average
negNP_sum <- negNP_stats %>% 
                group_by(treatment_type) %>%
                mutate(se = std.error(negNP)) %>% 
                summarise(avg = mean(negNP),
                          se = mean(se),
                          sd = sd(negNP))


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


### Models for NP ----
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

### Models for DR ----
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

### Average Light Response Curve Calculations ----
light <- read.csv("Data/full_lightresponses_revised.csv")

light_avg <- light %>% 
  group_by(Lcuv, treatment_type) %>% 
  summarise(avgCO2 = mean(CO2)) %>% 
  filter(Lcuv <= 1000)

summary(light_avg$avgCO2[light_avg$treatment_type == "control"])  # max control = 9.593
summary(light_avg$avgCO2[light_avg$treatment_type == "treatment"])  # max treatment = 29.287 

0.9*9.593  # 90% max control = 8.6337
0.9*29.287  # 90% max treatment = 26.3583

## Determining the intersection points for 90% LSP
# visualizing the 90% light saturation point 
(light_plot <- ggplot(light_avg, aes(x = Lcuv, y = avgCO2)) +
    geom_point(aes(color = treatment_type)) +
    geom_line(aes(color = treatment_type)) +
    facet_wrap(~treatment_type))

hline_light <- data.frame(z = c(8.6337, 26.3583), treatment_type = factor(c("control", "treatment")))

(lightsat_plot <- light_plot + 
    geom_hline(data = hline_light, aes(yintercept = z)))

# control intersection point
c1_l <- c(500, 7.1400000)
c2_l <- c(1000, 9.5933333)
control_y_l <- c(500, 8.6337)
control_y2_l <- c(1000, 8.6337)

line.line.intersection(c1_l, c2_l, control_y_l, control_y2_l, 
                       interior.only = FALSE)               # x = 804.4 µE m^-2 s^-1

# treatment intersection point
t1_l <- c(200, 21.2966667)
t2_l <- c(500, 29.2866667)
treatment_y_l <- c(200, 26.3583)
treatment_y2_l <- c(500, 26.3583)

line.line.intersection(t1_l, t2_l, treatment_y_l, treatment_y2_l, 
                       interior.only = FALSE)               # x = 390.0 µE m^-2 s^-1

## Determining the intersection points for LCP
# visualizing the light compensation point  
(lightcomp_plot <- light_plot + 
    geom_hline(yintercept = 0))

# control intersection point
c1_lc <- c(50, -2.2333333)
c2_lc <- c(100, 0.1766667)
control_y_lc <- c(50, 0)
control_y2_lc <- c(200, 0)

line.line.intersection(c1_lc, c2_lc, control_y_lc, control_y2_lc, 
                       interior.only = FALSE)               # x = 96.3 µE m^-2 s^-1

# treatment intersection point
t1_lc <- c(25, -0.1700000)
t2_lc <- c(50, 6.0200000)
treatment_y_lc <- c(20, 0)
treatment_y2_lc <- c(50, 0)

line.line.intersection(t1_lc, t2_lc, treatment_y_lc, treatment_y2_lc, 
                       interior.only = FALSE)               # x = 25.7 µE m^-2 s^-1


### Light Response Curve LSP Statistics ----
light_sum <- light %>% 
  group_by(sample) %>% 
  summarize(maxCO2 = max(CO2)) %>% 
  mutate(LSPx90 = 0.9*maxCO2)

## Determining the intersection points for 90% LSP
# visualizing the 90% light saturation points for all samples
(light_plot_facet <- ggplot(light, aes(x = Lcuv, y = CO2)) +
    geom_point(aes(color = treatment_type)) +
    geom_line(aes(color = treatment_type)) +
    facet_wrap(~sample) +
    geom_hline(data = light_sum, aes(yintercept = LSPx90)))

# C1 intersection point
c1_1 <- c(500, 7.14)
c1_2 <- c(1000, 9.59)
control_y_c1 <- c(500, 8.865)
control_y2_c1 <- c(1100, 8.865)

line.line.intersection(c1_1, c1_2, control_y_c1, control_y2_c1, 
                       interior.only = FALSE)               # x = 852.0 µE m^-2 s^-1

# C2 intersection point
c2_1 <- c(500, 4.33)
c2_2 <- c(1000, 5.36)
control_y_c2 <- c(500, 4.824)
control_y2_c2 <- c(1000, 4.824)

line.line.intersection(c2_1, c2_2, control_y_c2, control_y2_c2, 
                       interior.only = FALSE)               # x = 739.8 µE m^-2 s^-1

# C3 intersection point
c3_1 <- c(500, 9.95)
c3_2 <- c(1000, 13.83)
control_y_c3 <- c(500, 12.960)
control_y2_c3 <- c(1000, 12.960)

line.line.intersection(c3_1, c3_2, control_y_c3, control_y2_c3, 
                       interior.only = FALSE)               # x = 887.9 µE m^-2 s^-1

# T1 intersection point
t1_1 <- c(200, 21.30)
t1_2 <- c(500, 29.29)
treatment_y_t1 <- c(200, 26.361)
treatment_y2_t1 <- c(500, 26.361)

line.line.intersection(t1_1, t1_2, treatment_y_t1, treatment_y2_t1, 
                       interior.only = FALSE)               # x = 390.0 µE m^-2 s^-1

# T2 intersection point
t2_1 <- c(200, 17.54)
t2_2 <- c(500, 23.57)
treatment_y_t2 <- c(200, 21.213)
treatment_y2_t2 <- c(500, 21.213)

line.line.intersection(t2_1, t2_2, treatment_y_t2, treatment_y2_t2, 
                       interior.only = FALSE)               # x = 382.7 µE m^-2 s^-1

# T3 intersection point
t3_1 <- c(200, 25.05)
t3_2 <- c(500, 35.00)
treatment_y_t3 <- c(200, 34.920)
treatment_y2_t3 <- c(500, 34.920)

line.line.intersection(t3_1, t3_2, treatment_y_t3, treatment_y2_t3, 
                       interior.only = FALSE)               # x = 497.6 µE m^-2 s^-1

light_lsp <- light_sum
light_lsp$LSPcuv <- c(852.0, 739.8, 887.9, 390.0, 382.7, 497.6)
light_lsp$treatment_type <- c("control", "control", "control", 
                              "treatment", "treatment", "treatment")              

## T-Test for LSP
t.test(LSPcuv ~ treatment_type, data = light_lsp)  # p = 0.0025 (significant difference!)
# t = 6.9453, df = 3.8731 
# using the means of the calculated LSPcuv, control = 826.6 and treatment = 423.4
# (different to pre-averaged data)

## Determining SE of the average
light_sum_lsp <- light_lsp %>% 
  group_by(treatment_type) %>%
  mutate(se = std.error(LSPcuv)) %>% 
  summarise(avg = mean(LSPcuv),
            se = mean(se),
            sd = sd(LSPcuv))


### Light Response Curve LCP Statistics ----
# visualizing the light compensation point for all samples
(lightcomp_plot_facet <- ggplot(light, aes(x = Lcuv, y = CO2)) +
   geom_point(aes(color = treatment_type)) +
   geom_line(aes(color = treatment_type)) +
   facet_wrap(~sample) +
   geom_hline(yintercept = 0))

## Determining the intersection points for LCP
# C1 intersection point
c1_1b <- c(50, -2.23)
c1_2b <- c(100, 0.18)
control_y_c1b <- c(50, 0)
control_y2_c1b <- c(200, 0)

line.line.intersection(c1_1b, c1_2b, control_y_c1b, control_y2_c1b, 
                       interior.only = FALSE)               # x = 96.3 µE m^-2 s^-1

# C2 intersection point
c2_1b <- c(100, -1.25)
c2_2b <- c(200, 1.39)
control_y_c2b <- c(100, 0)
control_y2_c2b <- c(200, 0)

line.line.intersection(c2_1b, c2_2b, control_y_c2b, control_y2_c2b, 
                       interior.only = FALSE)               # x = 147.3 µE m^-2 s^-1

# C3 intersection point
c3_1b <- c(50, -1.18)
c3_2b <- c(100, 1.60)
control_y_c3b <- c(500, 0)
control_y2_c3b <- c(100, 0)

line.line.intersection(c3_1b, c3_2b, control_y_c3b, control_y2_c3b, 
                       interior.only = FALSE)               # x = 71.2 µE m^-2 s^-1

# T1 intersection point
t1_1b <- c(25, -0.17)
t1_2b <- c(50, 6.02)
treatment_y_t1b <- c(10, 0)
treatment_y2_t1b <- c(50, 0)

line.line.intersection(t1_1b, t1_2b, treatment_y_t1b, treatment_y2_t1b, 
                       interior.only = FALSE)               # x = 25.7 µE m^-2 s^-1

# T2 intersection point
t2_1b <- c(12, -2.27)
t2_2b <- c(25, 0.56)
treatment_y_t2b <- c(10, 0)
treatment_y2_t2b <- c(50, 0)

line.line.intersection(t2_1b, t2_2b, treatment_y_t2b, treatment_y2_t2b, 
                       interior.only = FALSE)               # x = 22.4 µE m^-2 s^-1

# T3 intersection point
t3_1b <- c(25, -0.90)
t3_2b <- c(50, 6.95)
treatment_y_t3b <- c(20, 0)
treatment_y2_t3b <- c(60, 0)

line.line.intersection(t3_1b, t3_2b, treatment_y_t3b, treatment_y2_t3b, 
                       interior.only = FALSE)               # x = 27.9 µE m^-2 s^-1

light_lcp <- light_sum
light_lcp$LCPcuv <- c(96.3, 147.3, 71.2, 25.7, 22.4, 27.9)

## T-Test for LCP
t.test(LCPcuv ~ treatment_type, data = light_lcp)  # p = 0.07 (NOT significantly different)
# t = 3.5464, DF = 2.0204 
# using the means of the calculated LCPcuv, control = 104.93 and treatment = 25.33 
# (slightly different to pre-averaged data)

## Determining SE of the average
light_sum_lcp <- light_lcp %>% 
  group_by(treatment_type) %>%
  mutate(se = std.error(LCPcuv)) %>% 
  summarise(avg = mean(LCPcuv),
            se = mean(se),
            sd = sd(LCPcuv))


### Optimal Water Content Ranges ----
water <- read.csv("Data/water_content_full.csv")

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
chlr <- read.csv("Data/chlorophyll.csv")

str(chlr)

t.test(chl_SA ~ type, data = chlr)  # p = 0.047 (is significantly different)
# t = -2.2996, DF = 8.7
# control = 736475.4, treatment = 1250650.3

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

summary(season_dr)   # slope = 0.122 ± 0.021
                     # p = 1.05e-5, t = 5.68 (significant)
anova(season_dr)

summary(season_np)   # slope = -0.013 ± 0.0065
                     # p = 0.061, t = -1.97 (NOT significant)
anova(season_np)

summary(season_gp)   # slope = 0.109 ± 0.0025
                     # p = 2.19e-4, t = 4.42 (significant)
anova(season_gp)

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


### Models for Climate (Temperature) ----
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
