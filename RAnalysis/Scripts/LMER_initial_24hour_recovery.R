#Title: LMER models Argopecten
#Author: Sam Gurr
#Edited by: Sam Gurr
#Date Last Modified: 20200914
#See Readme file for details
# script written by Sam J Gurr

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
library(lme4)
library(ggplot2)
library(nlme)
library(stargazer)
library(sjPlot)


#set working directory -------------------------------------------------------------- #
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")


# load datasets
HB.first.TWOhours <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/initial_hb_response/Summary_initial_hypoxia_exposure.csv", header=T) # first two hours O2 decline and initial hour of hypoxia 
HB.24hr.hypoxia <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/24-hour_hb_response/Summary_24_hour_exposure.csv", header=T) # full 24 hour mean data (every 10 minutes)  of O2 <= 2.0 mg L-1
HB.last.TWOhours <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/recovery_hb_response/Summary_recovery_to_normoxia.csv", header=T) # last two hours of O2 increase and initial hour of recovery to normoxia 
# ------------------------------------------------------------- 

# first hour of hypoxia
HB.hr.initial.hypoxia <- HB.first.TWOhours %>% dplyr::select(Site, hb_sensor, treatment, block, min_70, min_80, min_90, min_100, min_110, min_120) # initial hour of O2 <= 2.0 mg L-1 (every 10 minutes) 
HB.hr.initial_2<-  melt(HB.hr.initial.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.initial.OM <-  na.omit(HB.hr.initial_2)
HB.hr.initial.OM  <- HB.hr.initial.OM  %>% dplyr::select(-block) # ommit sensor and block

# 24 PROLONGED HYPOXIA -------------------------------------------------------------- #
HB.24hr.hypoxia # view data
HB.24hr.hypoxia_2 <-  melt(HB.24hr.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.24hr.hypoxia.OM <-  na.omit(HB.24hr.hypoxia_2)
HB.24hr.hypoxia.OM  <- HB.24hr.hypoxia.OM  %>% dplyr::select(-block) # ommit sensor and block


# RECOVERY -------------------------------------------------------------- #
HB.hr.recovery.hypoxia <- HB.last.TWOhours %>% dplyr::select(Site, hb_sensor, treatment, block, min_70, min_80, min_90, min_100, min_110, min_120) # initial hour of recovery after 24 hr hypoxia; O2 >= 6.0 mg L-1 (every 10 minutes) 
HB.hr.recovery.hypoxia_2 <-  melt(HB.hr.recovery.hypoxia, id.vars=c("Site", "hb_sensor", "treatment", "block"))
HB.hr.recovery.hypoxia.OM <-  na.omit(HB.hr.recovery.hypoxia_2)
HB.hr.recovery.hypoxia.OM  <- HB.hr.recovery.hypoxia.OM  %>% dplyr::select(-block) # ommit sensor and block

############# #
# PREP DATASETS FOR MODELS ####
############# #

# # INITIAL RESPONSE TO HYPOXIA - 6 YIMES EVERY 10 MINUTES OVER 1 HOUR
HB.hr.initial.OM$Subject_ID <- paste(HB.hr.initial.OM$Site, 
                                     HB.hr.initial.OM$hb_sensor, sep ="_")
HB.hr.initial.OM$treatment <- factor(HB.hr.initial.OM$treatment)
HB.hr.initial.OM$Site <- as.factor(HB.hr.initial.OM$Site)
HB.hr.initial.OM$Subject_ID <- as.factor(HB.hr.initial.OM$Subject_ID)
HB.hr.initial.OM$variable <- as.factor(HB.hr.initial.OM$variable)
HB.hr.initial.OM$value <- as.numeric(HB.hr.initial.OM$value)
HB.hr.initial.OM$treatment <- as.factor(HB.hr.initial.OM$treatment)
sapply(HB.hr.initial.OM, class)

# 24 HOURS OF HYPOXIA - 24 TIMES HOURLY MEANS OVER ONE DAY - Prepare variables
HB.24hr.hypoxia.OM$Subject_ID <- paste(HB.24hr.hypoxia.OM$Site, 
                                       HB.24hr.hypoxia.OM$hb_sensor, sep ="_")
HB.24hr.hypoxia.OM$Site <- as.factor(HB.24hr.hypoxia.OM$Site)
HB.24hr.hypoxia.OM$Subject_ID <- as.factor(HB.24hr.hypoxia.OM$Subject_ID)
HB.24hr.hypoxia.OM$Subject_ID <- abs(HB.24hr.hypoxia.OM$Subject_ID)
HB.24hr.hypoxia.OM$variable <- as.factor(HB.24hr.hypoxia.OM$variable)
HB.24hr.hypoxia.OM$treatment <- as.factor(HB.24hr.hypoxia.OM$treatment)
HB.24hr.hypoxia.OM$value <- as.numeric(HB.24hr.hypoxia.OM$value)
sapply(HB.24hr.hypoxia.OM, class)

# ONE HOUR OF RECPOVERY FROM HYPOXIA - 6 YIMES EVERY 10 MINUTES OVER 1 HOUR
HB.hr.recovery.hypoxia.OM$Subject_ID <- paste(HB.hr.recovery.hypoxia.OM$Site, 
                                              HB.hr.recovery.hypoxia.OM$hb_sensor, sep ="_")
HB.hr.recovery.hypoxia.OM$Site <- as.factor(HB.hr.recovery.hypoxia.OM$Site)
HB.hr.recovery.hypoxia.OM$Subject_ID <- as.factor(HB.hr.recovery.hypoxia.OM$Subject_ID)
HB.hr.recovery.hypoxia.OM$variable <- as.factor(HB.hr.recovery.hypoxia.OM$variable)
HB.hr.recovery.hypoxia.OM$value <- as.numeric(HB.hr.recovery.hypoxia.OM$value)
HB.hr.recovery.hypoxia.OM$treatment <- factor(HB.hr.recovery.hypoxia.OM$treatment)
sapply(HB.hr.recovery.hypoxia.OM, class)





################################################################################## #
################################################################################## #
# LME test ----------------------------------------------------------- #
################################################################################## 
################################################################################## #




# INITIAL DO DECLINE LME  -------------------------------------------------------- #
################################################################################## #




sapply(HB.hr.initial.OM,class)

# REML = F because we are comparing the models with different FIXED effects 

# # only the effect of random SITE as 'restricted' model for AIC comparison
# DOdecline_lme_restricted <- lmer(value ~ (1|Site), REML = F, data=HB.hr.initial.OM)
# anova(DOdecline_lme_restricted)

# FIXED effects of treatment and time  (mod1)
DOdecline_lme_mod1 <- lme(value~treatment*variable,random=~1|Site,data=HB.hr.initial.OM)
anova(DOdecline_lme_mod1)

leveneTest(residuals(DOdecline_lme_mod1) ~ HB.hr.initial.OM$treatment*HB.hr.initial.OM$variable) # pass homoscedasticity
boxplot(residuals(DOdecline_lme_mod1) ~ HB.hr.initial.OM$treatment*HB.hr.initial.OM$variable) # visual of residuals grouped
qqnorm(resid(DOdecline_lme_mod1)) # check for normality of residuals


stargazer(DOdecline_lme_mod1, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

library(emmeans)
emmeans(DOdecline_lme_mod1, list(pairwise ~ treatment), adjust = "tukey")
library(lmerTest)
DOdecline_lmer_mod1 <- lmer(value~treatment*variable + (1|Site),data=HB.hr.initial.OM)
anova(DOdecline_lmer_mod1)
difflsmeans(DOdecline_lmer_mod1, test.effs = "treatment")

# # FIXED effects of treatment ONLY (mod2)
# 
# DOdecline_lme_mod2 <- lmer(value~treatment + (1|Site), REML = F, data=HB.hr.initial.OM)
# anova(DOdecline_lme_mod2)
# 
# leveneTest(residuals(DOdecline_lme_mod2) ~ HB.hr.initial.OM$treatment) # pass homoscedasticity 
# boxplot(residuals(DOdecline_lme_mod2) ~ HB.hr.initial.OM$treatment) # visual of residuals grouped
# qqnorm(resid(DOdecline_lme_mod2)) # check for normality of residuals
# 
# # FIXED effects of time ONLY (mod3)
# DOdecline_lme_mod3 <- lmer(value~variable + (1|Site), REML = F, data=HB.hr.initial.OM)
# anova(DOdecline_lme_mod3)
# 
# leveneTest(residuals(DOdecline_lme_mod3) ~ HB.hr.initial.OM$variable) # pass homoscedasticity
# boxplot(residuals(DOdecline_lme_mod3) ~ HB.hr.initial.OM$variable) # visual of residuals grouped
# qqnorm(resid(DOdecline_lme_mod3)) # check for normality of residuals
# 
# # AIC comparisons
# #compute the AIC-corrected log-base-2 likelihood ratio (a.k.a. "bits" of evidence) - smaller the AIC the better the fit
# (AIC(DOdecline_lme_restricted)-AIC(DOdecline_lme_mod1))*log2(exp(1)) # mod1 treatment and time   == 182.9748
# (AIC(DOdecline_lme_restricted)-AIC(DOdecline_lme_mod2))*log2(exp(1)) # mod2 treatment            == 206.51
# (AIC(DOdecline_lme_restricted)-AIC(DOdecline_lme_mod3))*log2(exp(1)) # mod3 time                 == -12.7127







# 24 HOUR HYPOXIA ---------------------------------------------------------------- #
################################################################################## #




# REML = F because we are comparing the models with different FIXED effects 

# only the effect of random SITE as 'restricted' model for AIC comparison
# DO24HR_lme_restricted <- lmer(value ~ (1|Site), REML = F, data=HB.24hr.hypoxia.OM)
# anova(DO24HR_lme_restricted)

# FIXED effects of treatment and time  (mod1)
#DOdecline_lme_mod1 <- lme(value~treatment*variable,random=~1|Site,data=HB.hr.initial.hypoxia.RMANOVA.OM)
#anova(DOdecline_lme_mod1)

# DO24HR_lme_mod1 <- lmer(value~treatment*variable + (1|Site),  REML = F, data=HB.24hr.hypoxia.OM)
# anova(DO24HR_lme_mod1)

DO24HR_lme_mod1 <- lme(value~treatment*variable,random=~1|Site,data=HB.24hr.hypoxia.OM)
anova(DO24HR_lme_mod1)

leveneTest(residuals(DO24HR_lme_mod1) ~ HB.24hr.hypoxia.OM$treatment*HB.24hr.hypoxia.OM$variable) # pass homoscedasticity
boxplot(residuals(DO24HR_lme_mod1) ~ HB.24hr.hypoxia.OM$treatment*HB.24hr.hypoxia.OM$variable) # visual of residuals grouped
qqnorm(resid(DO24HR_lme_mod1)) # check for normality of residuals

stargazer(DO24HR_lme_mod1, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

emmeans(DO24HR_lme_mod1, list(pairwise ~ variable), adjust = "tukey") #posthoc


# # FIXED effects of treatment ONLY (mod2)
# DO24HR_lme_mod2 <- lmer(value~treatment + (1|Site), REML = F, data=HB.24hr.hypoxia.OM)
# anova(DO24HR_lme_mod2)
# 
# leveneTest(residuals(DO24HR_lme_mod2) ~ HB.24hr.hypoxia.OM$treatment) # does not pass homoscedasticity 
# boxplot(residuals(DO24HR_lme_mod2) ~ HB.24hr.hypoxia.OM$treatment) # visual of residuals grouped
# qqnorm(resid(DO24HR_lme_mod2)) # check for normality of residuals
# 
# # FIXED effects of time ONLY (mod3)
# DO24HR_lme_mod3 <- lmer(value~variable + (1|Site), REML = F, data=HB.24hr.hypoxia.OM)
# anova(DO24HR_lme_mod3)
# 
# leveneTest(residuals(DO24HR_lme_mod3) ~ HB.24hr.hypoxia.OM$variable) # pass homoscedasticity
# boxplot(residuals(DO24HR_lme_mod3) ~ HB.24hr.hypoxia.OM$variable) # visual of residuals grouped
# qqnorm(resid(DO24HR_lme_mod3)) # check for normality of residuals
# 
# # AIC comparisons
# #compute the AIC-corrected log-base-2 likelihood ratio (a.k.a. "bits" of evidence) - smaller the AIC the better the fit
# (AIC(DO24HR_lme_restricted)-AIC(DO24HR_lme_mod1))*log2(exp(1)) # mod1 treatment and timeB  == 825.9074
# (AIC(DO24HR_lme_restricted)-AIC(DO24HR_lme_mod2))*log2(exp(1)) # mod2 treatment            == 746.9893
# (AIC(DO24HR_lme_restricted)-AIC(DO24HR_lme_mod3))*log2(exp(1)) # mod3 time                 == 17.620951




# normoxia    RECOVERY       ----------------------------------------------------- #
################################################################################## #




# # REML = F because we are comparing the models with different FIXED effects 
# 
# # only the effect of random SITE as 'restricted' model for AIC comparison
# RECOVERY_lme_restricted <- lmer(value ~ (1|Site), REML = F, data=HB.hr.recovery.hypoxia.OM)
# anova(RECOVERY_lme_restricted)

# FIXED effects of treatment and time  (mod1)
#DOdecline_lme_mod1 <- lme(value~treatment*variable,random=~1|Site,data=HB.hr.initial.hypoxia.RMANOVA.OM)
#anova(DOdecline_lme_mod1)

RECOVERY_lme_mod1 <-  lme(value~treatment*variable,random=~1|Site,data=HB.hr.recovery.hypoxia.OM)
anova(RECOVERY_lme_mod1)

leveneTest(residuals(RECOVERY_lme_mod1) ~ HB.hr.recovery.hypoxia.OM$treatment*HB.hr.recovery.hypoxia.OM$variable) # pass homoscedasticity
boxplot(residuals(RECOVERY_lme_mod1) ~ HB.hr.recovery.hypoxia.OM$variable) # visual of residuals grouped
qqnorm(resid(RECOVERY_lme_mod1)) # check for normality of residuals

stargazer(RECOVERY_lme_mod1, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

emmeans(RECOVERY_lme_mod1, list(pairwise ~ treatment), adjust = "tukey") #posthoc

# # FIXED effects of treatment ONLY (mod2)
# RECOVERY_lme_mod2 <- lmer(value~treatment + (1|Site), REML = F, data=HB.hr.recovery.hypoxia.OM)
# anova(RECOVERY_lme_mod2)
# 
# leveneTest(residuals(RECOVERY_lme_mod2) ~ HB.hr.recovery.hypoxia.OM$treatment) # pass homoscedasticity
# boxplot(residuals(RECOVERY_lme_mod2) ~ HB.hr.recovery.hypoxia.OM$treatment) # visual of residuals grouped
# qqnorm(resid(RECOVERY_lme_mod2)) # check for normality of residuals
# 
# # FIXED effects of time ONLY (mod3)
# RECOVERY_lme_mod3 <- lmer(value~variable + (1|Site), REML = F, data=HB.hr.recovery.hypoxia.OM)
# anova(RECOVERY_lme_mod3)
# 
# leveneTest(residuals(RECOVERY_lme_mod3) ~ HB.hr.recovery.hypoxia.OM$variable) # pass homoscedasticity
# boxplot(residuals(RECOVERY_lme_mod3) ~ HB.hr.recovery.hypoxia.OM$variable) # visual of residuals grouped
# qqnorm(resid(RECOVERY_lme_mod3)) # check for normality of residuals
# 
# # AIC comparisons
# #compute the AIC-corrected log-base-2 likelihood ratio (a.k.a. "bits" of evidence) - smaller the AIC the better the fit
# (AIC(RECOVERY_lme_restricted)-AIC(RECOVERY_lme_mod1))*log2(exp(1)) # mod1 treatment and time == 121.3628
# (AIC(RECOVERY_lme_restricted)-AIC(RECOVERY_lme_mod2))*log2(exp(1)) # mod2 treatment          == 94.62552
# (AIC(RECOVERY_lme_restricted)-AIC(RECOVERY_lme_mod3))*log2(exp(1)) # mod3 time               == -14.01801








RECOVERY_lme_mod3.plot <- ggplot(HB.hr.recovery.hypoxia.OM, aes(x = treatment, y = value, colour = Site)) +
  facet_wrap(~variable, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(HB.hr.recovery.hypoxia.OM, pred = predict(RECOVERY_lme_mod1)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels
RECOVERY_lme_mod3.plot



library(sjPlot)

# 24 HOUR EXPOSURE LME
mod.24HOUR.LME <- lmer(value ~ treatment + (1|variable) + (1|Site), data = HB.24hr.hypoxia.OM) 
summary(mod.24HOUR.LME)
isSingular(mod.24HOUR.LME, tol = 1e-4)

install.packages("stargazer")
library(stargazer)
stargazer(mod.24HOUR.LME, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


plot_24hr_LME <- ggplot(HB.24hr.hypoxia.OM, aes(x = treatment, y = value, colour = Site)) +
  facet_wrap(~variable, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(HB.24hr.hypoxia.OM, pred = predict(mod.24HOUR.LME)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels
plot_24hr_LME

(re.effects <- plot_model(mod.24HOUR.LME, type = "re", show.values = TRUE)) # Visualise random effects 


# RECOVERY - DO INCRESAE
mod.RECOVERY.LME <- lmer(value ~ treatment + (1|variable) + (1|Site), data = HB.hr.recovery.hypoxia.OM) 
summary(mod.RECOVERY.LME)

stargazer(mod.RECOVERY.LME, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

plot_RECOVERY_LME <- ggplot(HB.hr.recovery.hypoxia.OM, aes(x = treatment, y = value, colour = Site)) +
  facet_wrap(~variable, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(HB.hr.recovery.hypoxia.OM, pred = predict(mod.RECOVERY.LME)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"))  # adding space between panels
plot_RECOVERY_LME
