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
#subset DAT.recovery
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
#######################################
# Two way anova between site and hb  ##
#######################################

twoway.anova.model <- aov(av_HB ~ Treatment *  site, data = initial.hour.ALL)
shapiro.test(residuals(twoway.anova.model)) # shaprio wilk test of model residuals p = 0.5429;normal distribution
hist((residuals(twoway.anova.model)))# histogram of model - looks normal
boxplot(residuals(twoway.anova.model)) #plot boxplot of residuals - some outliers present
plot(fitted(twoway.anova.model),residuals(twoway.anova.model)) # plot residuals
qqnorm(residuals(twoway.anova.model)) # qqplot - looks normal
summary(twoway.anova.model) # FIVE significant interaction between DAT.initiale and treatment
TukeyHSD(twoway.anova.model, 'Treatment', conf.level=0.95) # tukey test on the effect of treatment with 95% confidence

#######################################
#test sites seperately with T-TESTS ###   INITIAL EXPOSURE HYP VS. NORMOXIA
#######################################

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

# ONE-WAY ANOVA######################################################## #
par(mfrow=c(1,3)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots

# FI
t.test(av_HB ~ Treatment, data = FI.initial) # t.test 0.01848
FI.ONEWAY <-aov(av_HB ~ Treatment, data = FI.initial)
summary(FI.ONEWAY)
hist(residuals(FI.ONEWAY)) #plot histogram of residuals
boxplot(residuals(FI.ONEWAY)) #plot boxplot of residuals
plot(fitted(FI.ONEWAY),residuals(FI.ONEWAY))
qqnorm(residuals(FI.ONEWAY)) # qqplot
leveneTest(FI.ONEWAY) 
shapiro.test(residuals(FI.ONEWAY))
# X
t.test(av_HB ~ Treatment, data = X.initial) # t.test 0.01928
X.ONEWAY <-aov(av_HB ~ Treatment, data = X.initial)
summary(X.ONEWAY)
hist(residuals(X.ONEWAY)) #plot histogram of residuals
boxplot(residuals(X.ONEWAY)) #plot boxplot of residuals
plot(fitted(X.ONEWAY),residuals(X.ONEWAY))
qqnorm(residuals(X.ONEWAY)) # qqplot
leveneTest(X.ONEWAY) 
shapiro.test(residuals(X.ONEWAY))
# MBC
t.test(av_HB ~ Treatment, data = MBC.initial) # t.test 0.3864
MBC.ONEWAY <-aov(av_HB ~ Treatment, data = MBC.initial)
summary(MBC.ONEWAY)
hist(residuals(MBC.ONEWAY)) #plot histogram of residuals
boxplot(residuals(MBC.ONEWAY)) #plot boxplot of residuals
plot(fitted(MBC.ONEWAY),residuals(MBC.ONEWAY))
qqnorm(residuals(MBC.ONEWAY)) # qqplot
leveneTest(MBC.ONEWAY) 
shapiro.test(residuals(MBC.ONEWAY))
# NB
t.test(av_HB ~ Treatment, data = NB.initial) # t.test 0.05637
NB.ONEWAY <-aov(av_HB ~ Treatment, data = NB.initial)
summary(NB.ONEWAY)
hist(residuals(NB.ONEWAY)) #plot histogram of residuals
boxplot(residuals(NB.ONEWAY)) #plot boxplot of residuals
plot(fitted(NB.ONEWAY),residuals(NB.ONEWAY))
qqnorm(residuals(NB.ONEWAY)) # qqplot
leveneTest(NB.ONEWAY) 
shapiro.test(residuals(NB.ONEWAY))
# Q
t.test(av_HB ~ Treatment, data = Q.initial) # t.test 0.08077
Q.ONEWAY <-aov(av_HB ~ Treatment, data = Q.initial)
summary(Q.ONEWAY)
hist(residuals(Q.ONEWAY)) #plot histogram of residuals
boxplot(residuals(Q.ONEWAY)) #plot boxplot of residuals
plot(fitted(Q.ONEWAY),residuals(Q.ONEWAY))
qqnorm(residuals(Q.ONEWAY)) # qqplot
leveneTest(Q.ONEWAY) 
shapiro.test(residuals(Q.ONEWAY))


#######################################
#test sites seperately with T-TESTS ###   RECOVERRY TO NORMOXIA
#######################################
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

# FI
t.test(av_HB ~ Treatment, data = FI.recovery) # t.test 0.01848
FI.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = FI.recovery)
summary(FI.ONEWAY.recovery)
hist(residuals(FI.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(FI.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(FI.ONEWAY.recovery),residuals(FI.ONEWAY.recovery))
qqnorm(residuals(FI.ONEWAY.recovery)) # qqplot
leveneTest(FI.ONEWAY.recovery) 
shapiro.test(residuals(FI.ONEWAY.recovery))
# X
t.test(av_HB ~ Treatment, data = X.recovery) # t.test 0.01928
X.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = X.recovery)
summary(X.ONEWAY.recovery)
hist(residuals(X.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(X.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(X.ONEWAY.recovery),residuals(X.ONEWAY.recovery))
qqnorm(residuals(X.ONEWAY.recovery)) # qqplot
leveneTest(X.ONEWAY.recovery) 
shapiro.test(residuals(X.ONEWAY.recovery))
# MBC
t.test(av_HB ~ Treatment, data = MBC.recovery) # t.test 0.3864
MBC.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = MBC.recovery)
summary(MBC.ONEWAY.recovery)
hist(residuals(MBC.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(MBC.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(MBC.ONEWAY.recovery),residuals(MBC.ONEWAY.recovery))
qqnorm(residuals(MBC.ONEWAY.recovery)) # qqplot
leveneTest(MBC.ONEWAY.recovery) 
shapiro.test(residuals(MBC.ONEWAY.recovery))
# NB
t.test(av_HB ~ Treatment, data = NB.recovery) # t.test 0.05637
NB.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = NB.recovery)
summary(NB.ONEWAY.recovery)
hist(residuals(NB.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(NB.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(NB.ONEWAY.recovery),residuals(NB.ONEWAY.recovery))
qqnorm(residuals(NB.ONEWAY.recovery)) # qqplot
leveneTest(NB.ONEWAY.recovery) 
shapiro.test(residuals(NB.ONEWAY.recovery))
# Q
t.test(av_HB ~ Treatment, data = Q.recovery) # t.test 0.08077
Q.ONEWAY.recovery <-aov(av_HB ~ Treatment, data = Q.recovery)
summary(Q.ONEWAY.recovery)
hist(residuals(Q.ONEWAY.recovery)) #plot histogram of residuals
boxplot(residuals(Q.ONEWAY.recovery)) #plot boxplot of residuals
plot(fitted(Q.ONEWAY.recovery),residuals(Q.ONEWAY.recovery))
qqnorm(residuals(Q.ONEWAY.recovery)) # qqplot
leveneTest(Q.ONEWAY.recovery) 
shapiro.test(residuals(Q.ONEWAY.recovery))


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
                          y=av_HB, 
                          color=site)) + 
  geom_errorbar(aes(ymin=av_HB-se, 
                    ymax=av_HB+se),position = dodge)  +
  geom_point(position=position_dodge(0.8), fill = "white", shape=15, size=4)+
  theme_classic() +
  theme(
    axis.title.y = element_text(vjust= 1.8),
    axis.title.x = element_text(vjust= -0.5),
    axis.title = element_text(face = "bold")) +
  coord_cartesian(ylim=c(20,55), xlim=c(0,5)) +
  xlab("treatment") + 
  ylab("heatbeat rate") + 
  ggtitle("heartbeat rate response to nomroxia (hypoxia vs. control) STANDARD ERROR") +
  scale_color_manual(values = c("black", "blue", "red",
                                "green", "orange"))

plot4

# plot syntax WITH STANDARD DEVIATION
plot5  <- ggplot(table, aes(x=Treatment, 
                            y=av_HB, 
                            color=site)) + 
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
  ggtitle("heartbeat rate response to nomroxia (hypoxia vs. control) STANDARD DEVIATION") +
  scale_color_manual(values = c("black", "blue", "red",
                                "green", "orange"))

plot5

plot4 # stand error
plot5 # stand dev

