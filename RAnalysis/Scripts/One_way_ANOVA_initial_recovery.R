#Author: Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20190822
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("lmtest" %in% rownames(installed.packages()) == 'FALSE') install.packages('lmtest') 
if ("car" %in% rownames(installed.packages()) == 'FALSE') install.packages('car') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('car')
library('lmtest')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('car')
library('lmtest')
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(ggpubr)         # Version: 0.1.8 Date: 2018-08-30, Depends: R (>= 3.1.0), ggplot2, magrittrImports: ggrepel, grid, ggsci, stats, utils, tidyr, purrr, dplyr(>=0.7.1), cowplot, ggsignif, scales, gridExtra, glue, polynom
library(Rmisc)          # Version: 1.5 Packaged: 2013-10-21, Depends: lattice, plyr
library(plotrix)        # Version: 3.7-4, Date/Publication: 2018-10-03

# Set Working Directory:
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")

DAT.initial <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/initial_hb_response/Summary_initial_hypoxia_exposure_vertical.csv", header=T) #read Size.info DAT.initiala # choose 180202_initial_heartbeat_response_for_R
names(DAT.initial)
DAT_24_hr <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/24-hour_hb_response/Summary_24_hour_exposure.csv", header=T) #read Size.info DAT.initiala # choose 180202_initial_heartbeat_response_for_R
names(DAT_24_hr)
DAT.recovery <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/recovery_hb_response/Summary_recovery_to_normoxia_vertical.csv", header=T) #read Size.info DAT.initiala # choose 180202_initial_heartbeat_response_for_R
names(DAT.recovery)
names(DAT.recovery) <- c("site","mesocosm","site_trmt","Treatment","av_HB")

DAT.initial <- DAT.initial %>% select(site,Treatment,av_HB)
DAT.recovery <- DAT.recovery %>% select(site,Treatment,av_HB)
DAT.ALL <- rbind(DAT.initial,DAT.recovery)
################################ #
#subset DAT.initial
################################ #
DAT.ALL # view DAT.ALL
# both 1 hour in H and N (initial)
initial.hour.ALL <- DAT.ALL %>% 
  filter(Treatment %in% c('H_first_hour','N_first_hour'))
# just inital hour in H
initial.hour.HYP <-  DAT.ALL %>% 
  filter(Treatment %in% 'H_first_hour')
# just inital hour in N
initial.hour.NORM <-  DAT.ALL %>% 
  filter(Treatment %in% 'N_first_hour')
################################ #
#subseta
################################ #
DAT.ALL # view DAT.ALL
# both 1 hour in H and N (initial)
recovery.hour.ALL <- DAT.ALL %>% 
  filter(Treatment %in% c('Norm','HYP_post'))
# just inital hour in H
recovery.hour.HYP <-  DAT.ALL %>% 
  filter(Treatment %in% 'HYP_post')
# just inital hour in N
recovery.hour.NORM <-  DAT.ALL %>% 
  filter(Treatment %in% 'Norm')
# pre
pre.experiment <-  DAT.ALL %>% 
  filter(Treatment %in% 'H_pre_trmt_normoxia')
pre.experiment.ANOVA <- summary(aov(av_HB ~ site, data = pre.experiment))

#################################################### #
# One way ANOVA   INITIAL EXPOSURE HYP VS. NORMOXIA
#################################################### #

# FI
FI.initial <- dplyr::filter(initial.hour.ALL, site == "FI")
# X
X.initial <- dplyr::filter(initial.hour.ALL, site == "X")
# MBC
MBC.initial <- dplyr::filter(initial.hour.ALL, site == "MBC")
# NB
NB.initial <- dplyr::filter(initial.hour.ALL, site == "NB")
# Q
Q.initial <- dplyr::filter(initial.hour.ALL, site == "Q")

# ONE-WAY ANOVA    #################################################### #
par(mfrow=c(1,3)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots

# FI - ONE-WAY ANOVA Initial exoposure (1 hour mean)
FI.ONEWAY <-aov(av_HB ~ Treatment, data = FI.initial)
summary(FI.ONEWAY)
hist(residuals(FI.ONEWAY)) #plot histogram of residuals
boxplot(residuals(FI.ONEWAY)) #plot boxplot of residuals
plot(fitted(FI.ONEWAY),residuals(FI.ONEWAY))
qqnorm(residuals(FI.ONEWAY)) # qqplot
leveneTest(FI.ONEWAY) 
shapiro.test(residuals(FI.ONEWAY))
boxplot(av_HB ~ Treatment, data = FI.initial)

# Sag Harbor - ONE-WAY ANOVA Initial exoposure (1 hour mean)
X.ONEWAY <-aov(av_HB ~ Treatment, data = X.initial)
summary(X.ONEWAY)
hist(residuals(X.ONEWAY)) #plot histogram of residuals
boxplot(residuals(X.ONEWAY)) #plot boxplot of residuals
plot(fitted(X.ONEWAY),residuals(X.ONEWAY))
qqnorm(residuals(X.ONEWAY)) # qqplot
leveneTest(X.ONEWAY) 
shapiro.test(residuals(X.ONEWAY))
boxplot(av_HB ~ Treatment, data = X.initial)

# MBC - ONE-WAY ANOVA Initial exoposure (1 hour mean)
MBC.ONEWAY <-aov(av_HB ~ Treatment, data = MBC.initial)
summary(MBC.ONEWAY)
hist(residuals(MBC.ONEWAY)) #plot histogram of residuals
boxplot(residuals(MBC.ONEWAY)) #plot boxplot of residuals
plot(fitted(MBC.ONEWAY),residuals(MBC.ONEWAY))
qqnorm(residuals(MBC.ONEWAY)) # qqplot
leveneTest(MBC.ONEWAY) 
shapiro.test(residuals(MBC.ONEWAY))
boxplot(av_HB ~ Treatment, data = MBC.initial)

# NB - ONE-WAY ANOVA Initial exoposure (1 hour mean)
NB.ONEWAY <-aov(av_HB ~ Treatment, data = NB.initial)
summary(NB.ONEWAY)
hist(residuals(NB.ONEWAY)) #plot histogram of residuals
boxplot(residuals(NB.ONEWAY)) #plot boxplot of residuals
plot(fitted(NB.ONEWAY),residuals(NB.ONEWAY))
qqnorm(residuals(NB.ONEWAY)) # qqplot
leveneTest(NB.ONEWAY) 
shapiro.test(residuals(NB.ONEWAY))
boxplot(av_HB ~ Treatment, data = NB.initial)

# Q - ONE-WAY ANOVA Initial exoposure (1 hour mean)
Q.ONEWAY <-aov(av_HB ~ Treatment, data = Q.initial)
summary(Q.ONEWAY)
hist(residuals(Q.ONEWAY)) #plot histogram of residuals
boxplot(residuals(Q.ONEWAY)) #plot boxplot of residuals
plot(fitted(Q.ONEWAY),residuals(Q.ONEWAY))
qqnorm(residuals(Q.ONEWAY)) # qqplot
leveneTest(Q.ONEWAY) 
shapiro.test(residuals(Q.ONEWAY))
boxplot(av_HB ~ Treatment, data = Q.initial)

####################################### #
# One way ANOVA   RECOVERRY TO NORMOXIA
####################################### #
recovery.hour.ALL <- na.omit(recovery.hour.ALL)
# FI
FI.recovery <- dplyr::filter(recovery.hour.ALL, site == "FI")
# X
X.recovery <- dplyr::filter(recovery.hour.ALL, site == "X")
# MBC
MBC.recovery <- dplyr::filter(recovery.hour.ALL, site == "MB")
# NB
NB.recovery <- dplyr::filter(recovery.hour.ALL, site == "NB")
# Q
Q.recovery <- dplyr::filter(recovery.hour.ALL, site == "Q")

# ANOVAS######################################################## #

# FI - ONE-WAY ANOVA Recovery (1 hour mean)
FI.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = FI.recovery)
summary(FI.ONEWAY.recovery)
hist(residuals(FI.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(FI.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(FI.ONEWAY.recovery),residuals(FI.ONEWAY.recovery))
qqnorm(residuals(FI.ONEWAY.recovery)) # qqplot
leveneTest(FI.ONEWAY.recovery) 
shapiro.test(residuals(FI.ONEWAY.recovery))

# X - ONE-WAY ANOVA Recovery (1 hour mean)
X.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = X.recovery)
summary(X.ONEWAY.recovery)
hist(residuals(X.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(X.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(X.ONEWAY.recovery),residuals(X.ONEWAY.recovery))
qqnorm(residuals(X.ONEWAY.recovery)) # qqplot
leveneTest(X.ONEWAY.recovery) 
shapiro.test(residuals(X.ONEWAY.recovery))

# MBC - ONE-WAY ANOVA Recovery (1 hour mean)
MBC.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = MBC.recovery)
summary(MBC.ONEWAY.recovery)
hist(residuals(MBC.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(MBC.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(MBC.ONEWAY.recovery),residuals(MBC.ONEWAY.recovery))
qqnorm(residuals(MBC.ONEWAY.recovery)) # qqplot
leveneTest(MBC.ONEWAY.recovery) 
shapiro.test(residuals(MBC.ONEWAY.recovery))

# NB - ONE-WAY ANOVA Recovery (1 hour mean)
NB.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = NB.recovery)
summary(NB.ONEWAY.recovery)
hist(residuals(NB.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(NB.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(NB.ONEWAY.recovery),residuals(NB.ONEWAY.recovery))
qqnorm(residuals(NB.ONEWAY.recovery)) # qqplot
leveneTest(NB.ONEWAY.recovery) 
shapiro.test(residuals(NB.ONEWAY.recovery))

# Q - ONE-WAY ANOVA Recovery (1 hour mean)
Q.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = Q.recovery)
summary(Q.ONEWAY.recovery)
hist(residuals(Q.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(Q.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(Q.ONEWAY.recovery),residuals(Q.ONEWAY.recovery))
qqnorm(residuals(Q.ONEWAY.recovery)) # qqplot
leveneTest(Q.ONEWAY.recovery) 
shapiro.test(residuals(Q.ONEWAY.recovery))

############################################### #
# One Way ANOVA - 24 means 
############################################## #
DAT_24_hr

library(dplyr)


# FI
FI_24_hr <- DAT_24_hr %>% 
  mutate(av_HB = rowMeans(select(DAT_24_hr, starts_with("hr")), na.rm = TRUE)) %>% 
  dplyr::filter(Site == "FI")
# X
X_24_hr <- DAT_24_hr %>% 
  mutate(av_HB = rowMeans(select(DAT_24_hr, starts_with("hr")), na.rm = TRUE)) %>% 
  dplyr::filter(Site == "Sag Harbor")
# MBC
MBC_24_hr <- DAT_24_hr %>% 
  mutate(av_HB = rowMeans(select(DAT_24_hr, starts_with("hr")), na.rm = TRUE)) %>% 
  dplyr::filter(Site == "MB")
# NB
NB_24_hr<- DAT_24_hr %>% 
  mutate(av_HB = rowMeans(select(DAT_24_hr, starts_with("hr")), na.rm = TRUE)) %>% 
  dplyr::filter(Site == "Nicoll Bay")
# Q
Q_24_hr <- DAT_24_hr %>% 
  mutate(av_HB = rowMeans(select(DAT_24_hr, starts_with("hr")), na.rm = TRUE)) %>% 
  dplyr::filter(Site == "Quantuck")

# ANOVAS ######################################################## #

# FI - ONE-WAY ANOVA Recovery (1 hour mean)
FI.ONEWAY.24_hr <-aov(av_HB ~ treatment, data = FI_24_hr)
summary(FI.ONEWAY.24_hr)
hist(residuals(FI.ONEWAY.24_hr)) #plot histogram of residuals
boxplot(residuals(FI.ONEWAY.24_hr)) #plot boxplot of residuals
plot(fitted(FI.ONEWAY.24_hr),residuals(FI.ONEWAY.24_hr))
qqnorm(residuals(FI.ONEWAY.24_hr)) # qqplot
leveneTest(FI.ONEWAY.24_hr) 
shapiro.test(residuals(FI.ONEWAY.24_hr))

# X - ONE-WAY ANOVA Recovery (1 hour mean)
X.ONEWAY.24_hr<-aov(av_HB ~ treatment, data = X_24_hr)
summary(X.ONEWAY.24_hr)
hist(residuals(X.ONEWAY.24_hr)) #plot histogram of residuals
boxplot(residuals(X.ONEWAY.24_hr)) #plot boxplot of residuals
plot(fitted(X.ONEWAY.24_hr),residuals(X.ONEWAY.24_hr))
qqnorm(residuals(X.ONEWAY.24_hr)) # qqplot
leveneTest(X.ONEWAY.24_hr) 
shapiro.test(residuals(X.ONEWAY.24_hr))

# MBC - ONE-WAY ANOVA Recovery (1 hour mean)
MBC.ONEWAY.24_hr <-aov(av_HB ~ treatment, data = MBC_24_hr)
summary(MBC.ONEWAY.24_hr)
hist(residuals(MBC.ONEWAY.24_hr)) #plot histogram of residuals
boxplot(residuals(MBC.ONEWAY.24_hr)) #plot boxplot of residuals
plot(fitted(MBC.ONEWAY.24_hr),residuals(MBC.ONEWAY.24_hr))
qqnorm(residuals(MBC.ONEWAY.24_hr)) # qqplot
leveneTest(MBC.ONEWAY.24_hr) 
shapiro.test(residuals(MBC.ONEWAY.24_hr))

# NB - ONE-WAY ANOVA Recovery (1 hour mean)
NB.ONEWAY.24_hr <-aov(av_HB ~ treatment, data = NB_24_hr)
summary(NB.ONEWAY.24_hr)
hist(residuals(NB.ONEWAY.24_hr)) #plot histogram of residuals
boxplot(residuals(NB.ONEWAY.24_hr)) #plot boxplot of residuals
plot(fitted(NB.ONEWAY.24_hr),residuals(NB.ONEWAY.24_hr))
qqnorm(residuals(NB.ONEWAY.24_hr)) # qqplot
leveneTest(NB.ONEWAY.24_hr) 
shapiro.test(residuals(NB.ONEWAY.24_hr))

# Q - ONE-WAY ANOVA Recovery (1 hour mean)
Q.ONEWAY.24_hr <-aov(av_HB ~ treatment, data = Q_24_hr)
summary(Q.ONEWAY.24_hr)
hist(residuals(Q.ONEWAY.24_hr)) #plot histogram of residuals
boxplot(residuals(Q.ONEWAY.24_hr)) #plot boxplot of residuals
plot(fitted(Q.ONEWAY.24_hr),residuals(Q.ONEWAY.24_hr))
qqnorm(residuals(Q.ONEWAY.24_hr)) # qqplot
leveneTest(Q.ONEWAY.24_hr) 
shapiro.test(residuals(Q.ONEWAY.24_hr))

##########
# standard error plots by site for heartbeat rate
# pre and post recovery
# and normoxia
#######################
library(Rmisc)
library(ggplot2)

DAT.initialOM <- na.omit(DAT.initial)
DAT.initialOM
names(DAT.initialOM)

table= summarySE(DAT.initialOM, 
                 measurevar="av_HB", 
                 groupvars=c("Treatment","site"
                 ))
table

# make this to stagger points and see differences
dodge <- position_dodge(width=0.8)  
# plot syntax with STANDARD ERROR
plot4 <- ggplot(table, aes(x=Treatment, 
                          y=av_HB)) + 
  geom_errorbar(aes(ymin=av_HB-se, 
                    ymax=av_HB+se),position = dodge)  +
  geom_point(position=position_dodge(0.8), fill = "white", shape=15, size=4)+
  theme_bw() +
  theme(
    axis.title.y = element_text(vjust= 1.8),
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  coord_cartesian(ylim=c(20,55), xlim=c(0,5)) +
  xlab("treatment") + 
  ylab("heatbeat rate") + 
  facet_wrap(~site)

plot4 # view plot

# plot syntax WITH STANDARD DEVIATION
plot5 <- ggplot(table, aes(x=Treatment, 
                            y=av_HB)) + 
  geom_errorbar(aes(ymin=av_HB-sd, 
                    ymax=av_HB+sd),position = dodge)  +
  geom_point(position=position_dodge(0.8), fill = "white", shape=15, size=4)+
  theme_bw() +
  theme(
    axis.title.y = element_text(vjust= 1.8),
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  coord_cartesian(ylim=c(20,55), xlim=c(0,5)) +
  xlab("treatment") + 
  ylab("heatbeat rate") + 
  facet_wrap(~site)

plot5 # view plot


