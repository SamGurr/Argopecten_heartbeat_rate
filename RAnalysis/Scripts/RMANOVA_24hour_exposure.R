# script written by Sam J Gurr
# Last edit: 20200522
# RM ANOVA for the 24-hour hypoxia exposure

rm(list=ls()) # remove all data - start fresh

# Install packages if not already in your library-----------------------------------------------------------------------------------------------
require(car) # get the right sums of squares calculations
require(dplyr) # for manipulating our data
require(ggplot2) # for plotting and for our dataset
require(sjstats) # save us time computing key ANOVA stats beyond car
require(broom) # nice for making our results in neat tibbles
require(emmeans) # for marginal means calculations
# Load packages and pacage version/date/import/depends info
library(dplyr) # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(lsmeans) # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(ggpubr) 
library(reshape2) 
library(broom)
library(mvnormtest)
library(rlang)
library(car)
library(sjstats)
library(CGPfunctions)
library(rstatix)
library(ez)

#set working directory -------------------------------------------------------------- #
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")


# load datasets
HB.first.TWOhours <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/initial_hb_response/Summary_initial_hypoxia_exposure.csv", header=T) # first two hours O2 decline and initial hour of hypoxia 
HB.24hr.hypoxia <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/24-hour_hb_response/Summary_24_hour_exposure.csv", header=T) # full 24 hour mean data (every 10 minutes)  of O2 <= 2.0 mg L-1
HB.last.TWOhours <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/recovery_hb_response/Summary_recovery_to_normoxia.csv", header=T) # last two hours of O2 increase and initial hour of recovery to normoxia 
------------------------------------------------------------- #
# call the first hour O2 decline
HB.hr.O2.decline <-  HB.first.TWOhours %>%  dplyr::select(Site, hb_sensor, treatment, block, min_10, min_20, min_30, min_40, min_50, min_60) # data from the initial hour of O2 DECLINE to 2.0 mg L-1 (every 10 minutes) 
HB.hr.O2.decline.RMANOVA.OM <-  na.omit(HB.hr.O2.decline)
HB.hr.O2.decline.RMANOVA.OM <-  melt(HB.hr.O2.decline.RMANOVA.OM, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.O2.decline.RMANOVA.OM  <- HB.hr.O2.decline.RMANOVA.OM  %>% dplyr::select(-block) # ommit sensor and block
HB.hr.O2.decline <- melt(HB.hr.O2.decline, id.vars=c("Site", "hb_sensor", "treatment", "block")) # melt into a vertical dataset 'variable' = time; 'value' = heartbeat rate
HB.hr.O2.decline.2 <- HB.hr.O2.decline %>% dplyr::select(-hb_sensor, -block) # ommit sensor and block
# call the first hour of hypoxia
HB.hr.initial.hypoxia <- HB.first.TWOhours %>% dplyr::select(Site, hb_sensor, treatment, block, min_70, min_80, min_90, min_100, min_110, min_120) # initial hour of O2 <= 2.0 mg L-1 (every 10 minutes) 
HB.hr.initial.hypoxia.RMANOVA.OM <-  na.omit(HB.hr.initial.hypoxia)
HB.hr.initial.hypoxia.RMANOVA.OM <-  melt(HB.hr.initial.hypoxia.RMANOVA.OM, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.initial.hypoxia.RMANOVA.OM  <- HB.hr.initial.hypoxia.RMANOVA.OM  %>% dplyr::select(-block) # ommit sensor and block
HB.hr.initial.hypoxia <- melt(HB.hr.initial.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block")) # melt into a vertical dataset 'variable' = time; 'value' = heartbeat rate
HB.hr.initial.hypoxia.2 <- HB.hr.initial.hypoxia %>% dplyr::select(-hb_sensor, -block) # ommit sensor and block

# 24 PROLONGED HYPOXIA -------------------------------------------------------------- #
HB.24hr.hypoxia # view data
HB.24hr.hypoxia.OM <-  na.omit(HB.24hr.hypoxia)
HB.24hr.hypoxia.OM <-  melt(HB.24hr.hypoxia.OM, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.24hr.hypoxia.OM  <- HB.24hr.hypoxia.OM  %>% dplyr::select(-block) # ommit sensor and block
HB.24hr.hypoxia.2 <- melt(HB.24hr.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block")) # melt into a vertical dataset 'variable' = time; 'value' = heartbeat rate
HB.24hr.hypoxia.2 <- HB.24hr.hypoxia.2 %>% dplyr::select(-hb_sensor, -block) # ommit sensor and block

# ANOVA tests on Monebogue and Quantuck under prolonged hypoxia 
HB.24hr.hypoxia.ommit<-na.omit(HB.24hr.hypoxia.2) # ommit NAs
hourly.average <- HB.24hr.hypoxia.ommit %>%  # get the average data for sites 
  dplyr::select(treatment, variable, Site, value) %>% # call column to summarize 
  dplyr::filter(variable %in% 'hr_9') # NOTE "hr_9" starts sig diffs for Quantuck and hr_17 start sig diff for MB
hourly.average_MB <- hourly.average %>% filter(Site %in% 'MB') # set above to hour 17 - last hour not sig diff
summary(aov(lm(value~treatment, data = hourly.average_MB)))
hourly.average_Q <- hourly.average %>% filter(Site %in% 'Quantuck') # set above to hour 9 - last hour not sig diff
summary(aov(lm(value~treatment, data = hourly.average_Q)))


################################################################################## #
# RM ANOVA tests ----------------------------------------------------------- #
################################################################################## 

# RUN  ON SITES INDIVIDUALLY
HB.24hr.hypoxia # view data

HB.24hr.hypoxia.TEST <-  melt(HB.24hr.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.24hr.hypoxia.TEST <-  na.omit(HB.24hr.hypoxia.TEST)
HB.24hr.hypoxia.TEST  <- HB.24hr.hypoxia.TEST  %>% dplyr::select(-block) # ommit sensor and block
HB.24hr.hypoxia.TEST$Subject_ID <- paste(HB.24hr.hypoxia.TEST$Site, 
                                         HB.24hr.hypoxia.TEST$hb_sensor, sep ="_")
HB.24hr.hypoxia.TEST$Subject_ID <- as.factor(HB.24hr.hypoxia.TEST$Subject_ID)
# Fire island
FI_24hr <- HB.24hr.hypoxia.TEST %>% dplyr::filter(Site %in% "FI") # used dplyr to filter for FI
FI_24hr$Subject_ID <- as.factor(FI_24hr$Subject_ID)
sapply(FI_24hr, class)
colnames(FI_24hr)[5] <- "heartbeat"
colnames(FI_24hr)[4] <- "time"
ggboxplot(FI_24hr, x = "time", y = "heartbeat", color = "treatment")
table(FI_24hr$treatment, FI_24hr$time)
ggqqplot(FI_24hr, "heartbeat", ggtheme = theme_bw()) + facet_grid(time ~ treatment, labeller = "label_both")
FI_24hr %>% group_by(treatment, time) %>% identify_outliers(heartbeat) # ID outliers
Shapiro_FI<-FI_24hr %>% rstatix::group_by(treatment,time) %>% rstatix::shapiro_test(heartbeat) # shapiro test
lm(heartbeat~Subject_ID+time:treatment,data=FI_24hr)
FIres.aov <- rstatix::anova_test(data = FI_24hr, 
                        dv = heartbeat, 
                        wid = Subject_ID, 
                        within = c(treatment, time),
                        type=1)
FIres.aov


# Nicoll Bay
NB_24hr <- HB.24hr.hypoxia.OM %>%  dplyr::filter(Site %in% "Nicoll Bay") # used dplyr to filter for Nicoll
sapply(NB_24hr, class)
colnames(NB_24hr)[5] <- "heartbeat"
colnames(NB_24hr)[4] <- "time"
ggboxplot(NB_24hr, x = "time", y = "heartbeat", color = "treatment")
table(NB_24hr$treatment, NB_24hr$time)
ggqqplot(NB_24hr, "heartbeat", ggtheme = theme_bw()) + facet_grid(time ~ treatment, labeller = "label_both")
table(NB_24hr$treatment, NB_24hr$time)
NB_24hr %>% group_by(treatment, time) %>% identify_outliers(heartbeat) # ID outliers
Shapiro_NB<-NB_24hr %>% rstatix::group_by(treatment,time) %>% rstatix::shapiro_test(heartbeat) # shapiro test
lm(heartbeat~Subject_ID+time:treatment,data=NB_24hr)
NB_24hr.aov <- rstatix::anova_test(data = NB_24hr, 
                                 dv = heartbeat, 
                                 wid = Subject_ID, 
                                 within = c(treatment, time),
                                 type=1)
NB_24hr.aov

# Quantuck
Q_24hr <- HB.24hr.hypoxia.OM %>% dplyr::filter(Site %in% "Quantuck") # used dplyr to filter for Quantuck
sapply(Q_24hr, class)
colnames(Q_24hr)[5] <- "heartbeat"
colnames(Q_24hr)[4] <- "time"
ggboxplot(Q_24hr, x = "time", y = "heartbeat", color = "treatment")
table(Q_24hr$treatment, Q_24hr$time)
ggqqplot(Q_24hr, "heartbeat", ggtheme = theme_bw()) + facet_grid(time ~ treatment, labeller = "label_both")
table(Q_24hr$treatment, Q_24hr$time)
Q_24hr %>% group_by(treatment, time) %>% identify_outliers(heartbeat) # ID outliers
Shapiro_Q<-Q_24hr %>% rstatix::group_by(treatment,time) %>% rstatix::shapiro_test(heartbeat) # shapiro test
lm(heartbeat~Subject_ID+time:treatment,data=Q_24hr)
Q_24hr.aov <- rstatix::anova_test(data = Q_24hr, 
                                   dv = heartbeat, 
                                   wid = Subject_ID, 
                                   within = c(treatment, time),
                                   type=1)
Q_24hr.aov


# Sag Harbor
Sag_24hr <- HB.24hr.hypoxia.OM %>% dplyr::filter(Site %in% "Sag Harbor") # used dplyr to filter for Sag
sapply(Sag_24hr, class)
colnames(Sag_24hr)[5] <- "heartbeat"
colnames(Sag_24hr)[4] <- "time"
ggboxplot(Sag_24hr, x = "time", y = "heartbeat", color = "treatment")
ggqqplot(Sag_24hr, "heartbeat", ggtheme = theme_bw()) + facet_grid(time ~ treatment, labeller = "label_both")
table(Sag_24hr$treatment, Sag_24hr$time)
Sag_24hr %>% group_by(treatment, time) %>% identify_outliers(heartbeat) # ID outliers
Shapiro_sag<-Sag_24hr %>% rstatix::group_by(treatment,time) %>% rstatix::shapiro_test(heartbeat) # shapiro test
lm(heartbeat~Subject_ID+time:treatment,data=Sag_24hr)
Sag_24hr.aov <- rstatix::anova_test(data = Sag_24hr, 
                                   dv = heartbeat, 
                                   wid = Subject_ID, 
                                   within = c(treatment, time),
                                   type=1)
Sag_24hr.aov

# Moneybogue
MB_24hr <- HB.24hr.hypoxia.OM %>% dplyr::filter(Site %in% "MB") # used dplyr to filter for Mpneybogue
sapply(MB_24hr, class)
colnames(MB_24hr)[5] <- "heartbeat"
colnames(MB_24hr)[4] <- "time"
ggboxplot(MB_24hr, x = "time", y = "heartbeat", color = "treatment")
table(MB_24hr$treatment, MB_24hr$time)
ggqqplot(MB_24hr, "heartbeat", ggtheme = theme_bw()) + facet_grid(time ~ treatment, labeller = "label_both")
table(MB_24hr$treatment, MB_24hr$time)
MB_24hr %>% group_by(treatment, time) %>% identify_outliers(heartbeat) # ID outliers
Shapiro_MB<-MB_24hr %>% rstatix::group_by(treatment,time) %>% rstatix::shapiro_test(heartbeat) # shapiro test
lm(heartbeat~Subject_ID+time:treatment,data=MB_24hr)
MB_24hr.aov <- rstatix::anova_test(data = MB_24hr, 
                                   dv = heartbeat, 
                                   wid = Subject_ID, 
                                   within = c(treatment, time),
                                   type=1)
MB_24hr.aov


# ONE HOUR OF RECPOVERY FROM HYPOXIA - 6 YIMES EVERY 10 MINUTES OVER 1 HOUR
HB.hr.recovery.hypoxia.OM$Subject_ID <- paste(HB.hr.recovery.hypoxia.OM$Site, 
                                              HB.hr.recovery.hypoxia.OM$hb_sensor, sep ="_")
HB.hr.recovery.hypoxia.OM$Site <- as.factor(HB.hr.recovery.hypoxia.OM$Site)
HB.hr.recovery.hypoxia.OM$Subject_ID <- as.factor(HB.hr.recovery.hypoxia.OM$Subject_ID)
HB.hr.recovery.hypoxia.OM$Subject_ID <- as.numeric(HB.hr.recovery.hypoxia.OM$Subject_ID)
HB.hr.recovery.hypoxia.OM$variable <- as.factor(HB.hr.recovery.hypoxia.OM$variable)
HB.hr.recovery.hypoxia.OM$value <- as.numeric(HB.hr.recovery.hypoxia.OM$value)
HB.hr.recovery.hypoxia.OM$treatment <- factor(HB.hr.recovery.hypoxia.OM$treatment)
HB.hr.recovery.hypoxia.OM$treatment <- as.numeric(HB.hr.recovery.hypoxia.OM$treatment)

mod.RECOVERY.RMANOVA <- ezANOVA(data=HB.hr.recovery.hypoxia.OM, 
                              dv=value, 
                              wid=Subject_ID, 
                              within=variable, 
                              between=.(Site,treatment)) 

print(mod.RECOVERY.RMANOVA)
