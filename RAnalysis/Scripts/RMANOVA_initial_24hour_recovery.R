# script written by Sam J Gurr
# for MANOVA tests of Argopecten irradians heartbeat measurements in a 2016 experiment at SBU

# INTIAL  -
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

# test #### 
# Sig diff anovas - Sag Harbor and Fire Island sig diff at hour 1 of hypoxia and Nicoll Bay at hour 2
# Moneyboue and Wuantuck are delayed in diffference between treatment
HB.24hr.hypoxia.ommit<-na.omit(HB.24hr.hypoxia) # ommit NAs
hourly.average <- HB.24hr.hypoxia.ommit %>%  # get the average data for sites 
  dplyr::select(treatment, variable, Site, value) %>% # call column to summarize 
  dplyr::filter(variable %in% 'hr_17') # NOTE "hr_9" starts sig diffs for Quantuck and hr_17 start sig diff for MB
hourly.average_MB <- hourly.average %>% filter(Site %in% 'MB') # set above to hour 17
summary(aov(lm(value~treatment, data = hourly.average_MB)))
hourly.average_Q <- hourly.average %>% filter(Site %in% 'Quantuck') # set above to hour 9
summary(aov(lm(value~treatment, data = hourly.average_Q)))

#####
# RECOVERY -------------------------------------------------------------- #
# call the first hour of hypoxia
HB.hr.O2.increase <- HB.last.TWOhours %>%  dplyr::select(Site, hb_sensor, treatment, block, min_10, min_20, min_30, min_40, min_50, min_60) # data from the termination of 24 hour hypoxia; hour of O2 INCREASE to 6.0 mg L-1 (every 10 minutes)
HB.hr.O2.increase.OM <-  na.omit(HB.hr.O2.increase)
HB.hr.O2.increase.OM <-  melt(HB.hr.O2.increase.OM, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.O2.increase.OM  <- HB.hr.O2.increase.OM  %>% dplyr::select(-block) # ommit sensor and block
HB.hr.O2.increase <- melt(HB.hr.O2.increase, id.vars=c("Site", "hb_sensor", "treatment", "block")) # melt into a vertical dataset 'variable' = time; 'value' = heartbeat rate
HB.hr.O2.increase.2 <- HB.hr.O2.increase%>% dplyr::select(-hb_sensor, -block) # ommit sensor and block
# call the first hour of hypoxia
HB.hr.recovery.hypoxia <- HB.last.TWOhours %>% dplyr::select(Site, hb_sensor, treatment, block, min_70, min_80, min_90, min_100, min_110, min_120) # initial hour of recovery after 24 hr hypoxia; O2 >= 6.0 mg L-1 (every 10 minutes) 
HB.hr.recovery.hypoxia.OM <-  na.omit(HB.hr.recovery.hypoxia)
HB.hr.recovery.hypoxia.OM <-  melt(HB.hr.recovery.hypoxia.OM, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.recovery.hypoxia.OM  <- HB.hr.recovery.hypoxia.OM  %>% dplyr::select(-block) # ommit sensor and block
HB.hr.recovery.hypoxia <- melt(HB.hr.recovery.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block")) # melt into a vertical dataset 'variable' = time; 'value' = heartbeat rate
HB.hr.recovery.hypoxia.2 <- HB.hr.recovery.hypoxia%>% dplyr::select(-hb_sensor, -block) # ommit sensor and block

################################################################################## #
################################################################################## #
# RM ANOVA tests ----------------------------------------------------------- #
################################################################################## 
################################################################################## #

# ....response of heartbeat rate with
# one WIHTIN = time (as 'variable' in melted tables)
# two BETWEEN = Site and Treatment
# subject variable  = unique ID as Site_hbsensor

# INITIAL RESPONSE TO HYPOXIA - 6 YIMES EVERY 10 MINUTES OVER 1 HOUR
HB.hr.initial.hypoxia.RMANOVA.OM$Subject_ID <- paste(HB.hr.initial.hypoxia.RMANOVA.OM$Site, 
                                                     HB.hr.initial.hypoxia.RMANOVA.OM$hb_sensor, sep ="_")
HB.hr.initial.hypoxia.RMANOVA.OM$treatment <- factor(HB.hr.initial.hypoxia.RMANOVA.OM$treatment)
HB.hr.initial.hypoxia.RMANOVA.OM$Site <- as.factor(HB.hr.initial.hypoxia.RMANOVA.OM$Site)
HB.hr.initial.hypoxia.RMANOVA.OM$Subject_ID <- as.factor(HB.hr.initial.hypoxia.RMANOVA.OM$Subject_ID)
HB.hr.initial.hypoxia.RMANOVA.OM$Subject_ID <- as.numeric(HB.hr.initial.hypoxia.RMANOVA.OM$Subject_ID)
HB.hr.initial.hypoxia.RMANOVA.OM$variable <- as.factor(HB.hr.initial.hypoxia.RMANOVA.OM$variable)
HB.hr.initial.hypoxia.RMANOVA.OM$value <- as.numeric(HB.hr.initial.hypoxia.RMANOVA.OM$value)
HB.hr.initial.hypoxia.RMANOVA.OM$treatment <- as.numeric(HB.hr.initial.hypoxia.RMANOVA.OM$treatment)

mod.INITIAL.RMANOVA <- ezANOVA(data=HB.hr.initial.hypoxia.RMANOVA.OM, 
    dv=value, 
      wid=.(Subject_ID), 
        within=.(variable), 
         between=.(treatment,Site)) 

print(mod.INITIAL.RMANOVA)

# 24 HOURS OF HYPOXIA - 24 TIMES HOURLY MEANS OVER ONE DAY
HB.24hr.hypoxia.OM$Subject_ID <- paste(HB.24hr.hypoxia.OM$Site, 
                                       HB.24hr.hypoxia.OM$hb_sensor, sep ="_")
HB.24hr.hypoxia.OM$Site <- as.factor(HB.24hr.hypoxia.OM$Site)
HB.24hr.hypoxia.OM$Subject_ID <- as.factor(HB.24hr.hypoxia.OM$Subject_ID)
HB.24hr.hypoxia.OM$Subject_ID <- as.numeric(HB.24hr.hypoxia.OM$Subject_ID)
HB.24hr.hypoxia.OM$Subject_ID <- abs(HB.24hr.hypoxia.OM$Subject_ID)
HB.24hr.hypoxia.OM$variable <- as.factor(HB.24hr.hypoxia.OM$variable)
HB.24hr.hypoxia.OM$treatment <- as.factor(HB.24hr.hypoxia.OM$treatment)
HB.24hr.hypoxia.OM$value <- as.factor(HB.24hr.hypoxia.OM$value)
HB.24hr.hypoxia.OM$value <- as.numeric(HB.24hr.hypoxia.OM$value)

#which(is.nan(log(HB.24hr.hypoxia.OM$Subject_ID)))
typeof(HB.24hr.hypoxia.OM$value)

mod.24HOUR.RMANOVA <- ezANOVA(data=HB.24hr.hypoxia.OM, 
                              dv=value, 
                              wid=.(Subject_ID), 
                              within=.(variable), 
                              between=.(treatment,Site))
print(mod.24HOUR.RMANOVA)

?ezANOVA

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
