#Title: Cardiac_response (Figure 5 data and Supplementary Table 1)
#Author:Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20200125
#See Readme file for details

rm(list=ls()) #clears workspace 

# upload libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Set Working Directory:
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")

# laod table
heart <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/Heartbeat_rate_intial_exposure_recovery_MEAN.csv", header=T) #read Size.info data


DF_means <- heart %>% dplyr::select(-c(1:2, 11:16)) # select the subtracted (corrected) hreat rate values and the insitu mean characteristics by site
heart.rates <- DF_means %>% dplyr::select(c(1:8)) # exposure hrs 1,2,3,8,9,15,16; recovery hr 1
vars <- DF_means %>% dplyr::select(c(9:37)) # 30 total variables (5 DO thresholds * 6 characteristics; mag, freq, % time, count, EQ)

Summary_TABLE <- data.frame() # run this before the loop
for(i in 1:ncol(heart.rates)){
  for(m in 1:ncol(vars)) {
    # loop linear regressions for rquared and p values
    rsq <- summary(lm(heart.rates[,i] ~ vars[,m]))$r.squared
    pval <- summary(lm(heart.rates[,i] ~ vars[,m]))$coef[2,"Pr(>|t|)"]
    # assign the data table 
    RSQ.loop <- data.frame(matrix(nrow = 1, ncol = 4)) # create a new data table
    colnames(RSQ.loop) <- c('heart.meas', 'hypoxia.characteristic', 'rsq', 'pval') # assign headers
    RSQ.loop$heart.meas <- colnames(heart.rates[i])
    RSQ.loop$hypoxia.characteristic <- colnames(vars[m])
    RSQ.loop$rsq <- rsq
    RSQ.loop$pval <- pval
    # loop additions 
    df <- data.frame(RSQ.loop) # name dataframe for this single row
    Summary_TABLE <- rbind(Summary_TABLE,df) # bind to a cumulative list dataframe
    }# inside loop'
print(Summary_TABLE) # show loop progress in the console
}# outside loo
Summary_TABLE # loop product 

### - create new summary table to filter out unwanted data (salininty, secchi depth, chl a, site depth)
#### Why? This data was not measured continuously during field acclimation and is not related to the 5 hypoxia thresholds for the heatmap
### - use tidyr::separate to divide the 'hypoxia.characteristic' to 'Descriptor' & 'DO_threshold' by 
Summary_TABLE2 <- Summary_TABLE %>% 
  dplyr::filter(!hypoxia.characteristic  %in% c("Salinity", "secchi_depth", "site_depth", "chl_a_ug_l")) %>% 
  separate(hypoxia.characteristic, c("Descriptor", "DO_threshold"), "_")

# Create Heatmap
Summary_TABLE2$rsq <- format(round(Summary_TABLE2$rsq, 1)) # sig figs for heat map
Summary_TABLE2$rsq <-as.numeric(Summary_TABLE2$rsq) # convert rsq to numeric
Summary_TABLE2$heart.meas = factor(Summary_TABLE2$heart.meas, 
                                   levels=c("SUB_1_hr_init", "SUB_2_hr_init", "SUB_3_hr_init", "SUB_8_hr", 
                                            "SUB_9_hr", "SUB_15_hr","SUB_16_hr", "Recovery_1_hr")) # order the hearbeat rate variables
RSQ_Heatmap_Plots <- Summary_TABLE2 %>% # un heat map
  #format(round(rsq, 2), nsmall = 2) %>% 
  #dplyr::filter(heart.meas %in% c("SUB_1_hr_init", "SUB_2_hr_init", "SUB_3_hr_init", "Recovery_1_hr")) %>% 
  #dplyr::filter(!Descriptor %in% c("MeanConcentration", "EQ")) %>% 
  ggplot(aes(Descriptor,DO_threshold, fill = rsq)) + # the descritor and DO threhsolds and fill by the rsq values
  geom_tile() + # tile command calls the data as a heat map
  scale_fill_gradient(low = "white", high = "blue") + # heatmap as white to grey
  facet_grid(. ~ heart.meas, scales = "free") + # heat map spectrum is determined by the data
  geom_text(size=2, aes(label = rsq)) + # add rsq text to each heatmap segment
  theme(axis.text.x = element_text(angle = 90)) # rotate x axis labels
RSQ_Heatmap_Plots # view plots
# save the heat map to the output folder
ggsave(RSQ_Heatmap_Plots, file="C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/Output/Supplementary_heatmap.pdf",
       width=35, height=20, units = "cm", dpi=500)

################################################################# #
###### BELOW IS A MANUAL METHOD TO VIEW ALL LINEAR REGRESSIONS ## #
###################################################################

################################################### #
#### SITE CHARACTERISTICS ######################### #
################################################### #

site<-heart$site
ch<-heart$chl_a_ug_l
sal<-heart$Salinity
t<-heart$secchi_depth

##HyPOXIA Frequency (count = f) 
###Duration (d)
###mean concentration (magnitude) = m
###  frequency * mean duration / mean magnitude = EQ

f4.8<-heart$count_4.8
m4.8<-heart$MeanConcentration_4.8
d4.8<-heart$MeanDuration_4.8
EQ4.8<-heart$EQ_4.8

f3.5<-heart$count_3.5
m3.5<-heart$MeanConcentration_3.5
d3.5<-heart$MeanDuration_3.5
EQ3.5<-heart$EQ_3.5

f3.0<-heart$count_3.0
m3.0<-heart$MeanConcentration_3.0
d3.0<-heart$MeanDuration_3.0
EQ3.0<-heart$EQ_3.0

f2.0<-heart$count_2.0
m2.0<-heart$MeanConcentration_2.0
d2.0<-heart$MeanDuration_2.0
EQ2.0<-heart$EQ_2.0

f1.0<-heart$count_1.0
m1.0<-heart$MeanConcentration_1.0
d1.0<-heart$MeanDuration_1.0
EQ1.0<-heart$EQ_1.0

## PERCENT TIME OF HYPOheartX PERIODS
PERC4.8<- heart$Percent_4.8
PERC3.5<-heart$Percent_3.5
PERC3.0<-heart$Percent_3.0
PERC2.0<-heart$Percent_2.0
PERC1.0<-heart$Percent_1.0

### LAST WEEK of hypoxia ####
### Note this is only the PERCENT TIME SEVERE (<2.0)
##  and  PERCENT TIME MODERATE (<4.8)

week_sev<-heart$Percent_2.0_last_week
week_mod<-heart$Percent_4.8_last_week

################################################### #
#### Heartbeat data input ######################### #
################################################### #

names(heart)

# mean SUB value of
# (mean prebasal HB - mean HB response)
SUB_1<-heart$SUB_1_hr_init
SUB_2<-heart$SUB_2_hr_init
SUB_3<-heart$SUB_3_hr_init

# raw heartbeat response 
# NOT corrected for basal heartbeat rate
one_hr<-heart$X1_hr_init
two_hr<-heart$X2_hr_init
three_hr<-heart$X3_hr_init

#reccovery (only at 1 hour after normoxia reestablished)
# normoxia meaning >6.0 mg L-1
hr_recov<-heart$Recovery_1_hr


################################################### #
#### 1 hour post hypoxia ########################## #
################################################### #

# CARDIAC RESPONSE  ABOSOLUTE VALUE FROM PREBASAL CORRECTION ####### #

par(mfrow=c(3,2))
p4_1A<-lm(SUB_1~PERC4.8) #Multiple R-squared:  0.5039,	Adjusted R-squared:  0.3385 
summary(p4_1A)
plot(PERC4.8, SUB_1)
abline(p4_1A)

p3.5_1A<-lm(SUB_1~PERC3.5)  # Multiple R-squared:  0.7214,	Adjusted R-squared:  0.6285    LOOK AT THIS ONE
summary(p3.5_1A)
plot(PERC3.5, SUB_1)
abline(p3.5_1A)

p3_1A<-lm(SUB_1~PERC3.0) # Multiple R-squared:  0.5825,	Adjusted R-squared:  0.4434 
summary(p3_1A)
plot(PERC3.0, SUB_1)
abline(p3_1A)

p2_1A<-lm(SUB_1~PERC2.0) #Multiple R-squared:  0.4847,	Adjusted R-squared:  0.3129
summary(p2_1A)
plot(PERC2.0, SUB_1)
abline(p2_1A)

p1_1A<-lm(SUB_1~PERC1.0) #Multiple R-squared:  0.4283,	Adjusted R-squared:  0.2377
summary(p1_1A)
plot(PERC1.0, SUB_1)
abline(p1_1A)


#####  PERCENTTIME OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(PERC3.5, SUB_1)
abline(p3.5_1A)
summary(p3.5_1A)
mylabelp3.5_1A <- bquote(italic(R)^2 == .(format(0.7214, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_1A)


############################################################################## #
############################################################################## #
############# TEST ON THE FREQUENCY OR COUNT OF DATA  ######################## #
###############################################################################

# CARDIAC RESPONSE - SUBTRACTED FROM BASAL RATE(S)

par(mfrow=c(3,2))
c4_1A<-lm(SUB_1~f4.8) #Multiple R-squared:  0.3621,	Adjusted R-squared:  0.1495 
plot(f4.8, SUB_1)
summary(c4_1A)
abline(c4_1A)

c3.5_1A<-lm(SUB_1~f3.5) #Multiple R-squared:  0.6615,	Adjusted R-squared:  0.5487 Look at this one
plot(f3.5, SUB_1)
summary(c3.5_1A)
abline(c3.5_1A)

c3_1A<-lm(SUB_1~f3.0) #Multiple R-squared:  0.6394,	Adjusted R-squared:  0.5192
plot(f3.0, SUB_1)
summary(c3_1A)
abline(c3_1A)

c2_1A<-lm(SUB_1~f2.0) #Multiple R-squared:  0.5891,	Adjusted R-squared:  0.4522 
plot(f2.0, SUB_1)
summary(c2_1A)
abline(c2_1A)

c1_1A<-lm(SUB_1~f1.0) #Multiple R-squared:  0.5638,	Adjusted R-squared:  0.4184
plot(f1.0, SUB_1)
summary(c1_1A)
abline(c1_1A)

#####  Count OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(f3.5, SUB_1)
abline(c3.5_1A)
summary(c3.5_1A)
mylabelDUR2_1A <- bquote(italic(R)^2 == .(format(0.6615, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_1A)


###############################################################################
################     TEST OF THE DURATION DATA       ##########################
###############################################################################

### subtracted value 

par(mfrow=c(3,2))

DUR4.8_1A<-lm(SUB_1~d4.8) #Multiple R-squared:  0.5154,	Adjusted R-squared:  0.3538
plot(d4.8, SUB_1)
summary(DUR4.8_1A)
abline(DUR4.8_1A)

DUR3.5_1A<-lm(SUB_1~d3.5) #Multiple R-squared:  0.7262,	Adjusted R-squared:  0.6349 
plot(d3.5, SUB_1)
summary(DUR3.5_1A)
abline(DUR3.5_1A)

DUR3_1A<-lm(SUB_1~d3.0) #Multiple R-squared:  0.7475,	Adjusted R-squared:  0.6634
plot(d3.0, SUB_1)
summary(DUR3_1A)
abline(DUR3_1A)

DUR2_1A<-lm(SUB_1~d2.0)  # Multiple R-squared:  0.7878,	Adjusted R-squared:  0.717 Highest - look at this one
plot(d2.0, SUB_1)
summary(DUR2_1A)
abline(DUR2_1A)

DUR1_1A<-lm(SUB_1~d1.0) #Multiple R-squared:  0.7363,	Adjusted R-squared:  0.6484 
plot(d1.0, SUB_1)
summary(DUR1_1A)
abline(DUR1_1A)

#####  DURAITON OF EVENTS BELOW 2.0 HAS THE STRONGEST R SQUARED OF 0.8995

par(mfrow=c(1,1))

plot(d2.0, SUB_1)
abline(DUR2_1A)
summary(DUR2_1A)
mylabelDUR2_1A <- bquote(italic(R)^2 == .(format(0.7878, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR2_1A)


##### tEST OF THE mAGNITUDE dATA###########################################

##  subtract from basal

par(mfrow=c(3,2))

MAG4.8_1A<-lm(SUB_1~m4.8) # Multiple R-squared:  0.6897,	Adjusted R-squared:  0.5862
plot(m4.8, SUB_1)
summary(MAG4.8_1A)
abline(MAG4.8_1A)

MAG3.5_1A<-lm(SUB_1~m3.5) #Multiple R-squared:  0.8325,	Adjusted R-squared:  0.7766 
plot(m3.5, SUB_1)
summary(MAG3.5_1A)
abline(MAG3.5_1A)

MAG3_1A<-lm(SUB_1~m3.0) #Multiple R-squared:  0.8421,	Adjusted R-squared:  0.7895 
plot(m3.0, SUB_1)
summary(MAG3_1A)
abline(MAG3_1A)

MAG2_1A<-lm(SUB_1~m2.0)  # Multiple R-squared:  0.9902,	Adjusted R-squared:  0.9869
plot(m2.0, SUB_1)
summary(MAG2_1A)
abline(MAG2_1A)

MAG1_1A<-lm(SUB_1~m1.0) # Multiple R-squared:  0.9524,	Adjusted R-squared:  0.9366
plot(m1.0, SUB_1)
summary(MAG1_1A)
abline(MAG1_1A)


#### MAGNITUDE OF HYPOXIA AT 4.8 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(m2.0, SUB_1)
abline(MAG2_1A)
summary(MAG2_1A)
mylabelMAG2_1A<- bquote(italic(R)^2 == .(format(0.9902 , digits = 4)))
text(x = 0.5, y = 8, labels = mylabelMAG2_1A)


#############        EQUATION       ###################################
#######################################################################

# subtracted from basal
par(mfrow=c(3,2))

EQ4_1A<-lm(SUB_1~EQ4.8) #Multiple R-squared:  0.521,	Adjusted R-squared:  0.3614
plot(EQ4.8, SUB_1)
summary(EQ4_1A)

EQ3.5_1A<-lm(SUB_1~EQ3.5)  #Multiple R-squared:  0.5653,	Adjusted R-squared:  0.4204  Look at this one
plot(EQ3.5, SUB_1)
summary(EQ3.5_1A)
abline(EQ3.5_1A)

EQ3_1A<-lm(SUB_1~EQ3.0) #Multiple R-squared:  0.563,	Adjusted R-squared:  0.4174
plot(EQ3.0, SUB_1)
summary(EQ3_1A)

EQ2_1A<-lm(SUB_1~EQ2.0) #Multiple R-squared:  0.4802,	Adjusted R-squared:  0.3069 
plot(EQ2.0, SUB_1)
summary(EQ2_1A)

EQ1_1A<-lm(SUB_1~EQ1.0) #Multiple R-squared:  0.4208,	Adjusted R-squared:  0.2277 
plot(EQ1.0, SUB_1)
summary(EQ1_1A)


#### EQUATION  OF HYPOXIA AT 3.5 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(EQ3.5, SUB_1)
abline(EQ3.5_1A)
summary(EQ3.5_1A)
mylabelEQ3.5_1A<- bquote(italic(R)^2 == .(format(0.5653 , digits = 4)))
text(x = 0.5, y = 8, labels = mylabelEQ3.5_1A)


#############       finaL FIGURE ONE HOUR STRONGEST REGRESSIONS  ###################################
####################################################################################################


#####  PERCENTTIME OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(3,2))

#PERCENT TIME 

plot(PERC3.5, SUB_1)
abline(p3.5_1A)
summary(p3.5_1A)
mylabelp3.5_1A <- bquote(italic(R)^2 == .(format(0.7214, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_1A)

#  EQUATION 

plot(EQ3.5, SUB_1)
abline(EQ3.5_1A)
summary(EQ3.5_1A)
mylabelEQ3.5_1A<- bquote(italic(R)^2 == .(format(0.5653 , digits = 4)))
text(x = 0.5, y = 8, labels = mylabelEQ3.5_1A)

# MAGNITUDE

plot(m2.0, SUB_1)
abline(MAG2_1A)
summary(MAG2_1A)
mylabelMAG2_1A<- bquote(italic(R)^2 == .(format(0.9902 , digits = 4)))
text(x = 0.5, y = 8, labels = mylabelMAG2_1A)

# DURATION

plot(d2.0, SUB_1)
abline(DUR2_1A)
summary(DUR2_1A)
mylabelDUR2_1A <- bquote(italic(R)^2 == .(format(0.7878, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR2_1A)

# COUNT

plot(f3.5, SUB_1)
abline(c3.5_1A)
summary(c3.5_1A)
mylabelDUR2_1A <- bquote(italic(R)^2 == .(format(0.6615, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_1A)


######################################################################################################################################## #
######################################################################################################################################## #
# TWO HOUR
######################################################################################################################################## #
#########################################################################################################################################

####TEST THE PERCENT FACTOR ON heartbeat

#######################################   ABOSOLUTE VALUE FROM PREBASAL CORRECTION #######
par(mfrow=c(3,2))
p4_2A<-lm(SUB_2~PERC4.8) #Multiple R-squared:  0.8317,	Adjusted R-squared:  0.7757
summary(p4_2A)
plot(PERC4.8, SUB_2)
abline(p4_2A)

p3.5_2A<-lm(SUB_2~PERC3.5)  # Multiple R-squared:  0.9582,	Adjusted R-squared:  0.9442 
summary(p3.5_2A)
plot(PERC3.5, SUB_2)
abline(p3.5_2A)

p3_2A<-lm(SUB_2~PERC3.0) # Multiple R-squared:  0.8802,	Adjusted R-squared:  0.8403 
summary(p3_2A)
plot(PERC3.0, SUB_2)
abline(p3_2A)

p2_2A<-lm(SUB_2~PERC2.0) #Multiple R-squared:  0.778,	Adjusted R-squared:  0.704
summary(p2_2A)
plot(PERC2.0, SUB_2)
abline(p2_2A)

p1_2A<-lm(SUB_2~PERC1.0) #Multiple R-squared:  0.6903,	Adjusted R-squared:  0.5871 
summary(p1_2A)
plot(PERC1.0, SUB_2)
abline(p1_2A)


#####  PERCENTTIME OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(PERC3.5, SUB_2)
abline(p3.5_2A)
summary(p3.5_2A)
mylabelp3.5_2A <- bquote(italic(R)^2 == .(format(0.9582, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_2A)


############################################################################## #
############################################################################## #
############# TEST ON THE FREQUENCY OR COUNT OF DATA  ######################## #
###############################################################################


par(mfrow=c(3,2))
c4_2A<-lm(SUB_2~f4.8) #Multiple R-squared:  0.6874,	Adjusted R-squared:  0.5832 
plot(f4.8, SUB_2)
summary(c4_2A)
abline(c4_2A)

c3.5_2A<-lm(SUB_2~f3.5) #Multiple R-squared:  0.9307,	Adjusted R-squared:  0.9076  Look at this one
plot(f3.5, SUB_2)
summary(c3.5_2A)
abline(c3.5_2A)

c3_2A<-lm(SUB_2~f3.0) #Multiple R-squared:  0.9128,	Adjusted R-squared:  0.8837
plot(f3.0, SUB_2)
summary(c3_2A)
abline(c3_2A)

c2_2A<-lm(SUB_2~f2.0) #Multiple R-squared:  0.8736,	Adjusted R-squared:  0.8315 
plot(f2.0, SUB_2)
summary(c2_2A)
abline(c2_2A)

c1_2A<-lm(SUB_2~f1.0) #Multiple R-squared:  0.8802,	Adjusted R-squared:  0.8403 
plot(f1.0, SUB_2)
summary(c1_2A)
abline(c1_2A)

#####  Count OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(f3.5, SUB_2)
abline(c3.5_2A)
summary(c3.5_2A)
mylabelDUR2_2A <- bquote(italic(R)^2 == .(format(0.9307, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_2A)


############################################################################## #
################     TEST OF THE DURATION DATA       ######################### #
###############################################################################

par(mfrow=c(3,2))

DUR4.8_2A<-lm(SUB_2~d4.8) #Multiple R-squared:  0.8494,	Adjusted R-squared:  0.7992 
plot(d4.8, SUB_2)
summary(DUR4.8_2A)
abline(DUR4.8_2A)

DUR3.5_2A<-lm(SUB_2~d3.5) #Multiple R-squared:  0.9726,	Adjusted R-squared:  0.9635 
plot(d3.5, SUB_2)
summary(DUR3.5_2A)
abline(DUR3.5_2A)

DUR3_2A<-lm(SUB_2~d3.0) #Multiple R-squared:  0.9793,	Adjusted R-squared:  0.9724
plot(d3.0, SUB_2)
summary(DUR3_2A)
abline(DUR3_2A)

DUR2_2A<-lm(SUB_2~d2.0)  # Multiple R-squared:  0.9352,	Adjusted R-squared:  0.9136 
plot(d2.0, SUB_2)
summary(DUR2_2A)
abline(DUR2_2A)

DUR1_2A<-lm(SUB_2~d1.0) #Multiple R-squared:  0.9061,	Adjusted R-squared:  0.8748 
plot(d1.0, SUB_2)
summary(DUR1_2A)
abline(DUR1_2A)

#####  DURAITON OF EVENTS BELOW 2.0 HAS THE STRONGEST R SQUARED OF 0.8995

par(mfrow=c(1,1))

plot(d3.0, SUB_2)
abline(DUR3_2A)
summary(DUR3_2A)
mylabelDUR3_2A <- bquote(italic(R)^2 == .(format(0.9793, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR3_2A)


##### tEST OF THE mAGNITUDE dATA###########################################

## asbosulte subtract from basal

par(mfrow=c(3,2))

MAG4.8_2A<-lm(SUB_2~m4.8) #Multiple R-squared:  0.9591,	Adjusted R-squared:  0.9455 
plot(m4.8, SUB_2)
summary(MAG4.8_2A)
abline(MAG4.8_2A)

MAG3.5_2A<-lm(SUB_2~m3.5) #Multiple R-squared:  0.5513,	Adjusted R-squared:  0.4017 
plot(m3.5, SUB_2)
summary(MAG3.5_2A)
abline(MAG3.5_2A)

MAG3_2A<-lm(SUB_2~m3.0) #Multiple R-squared:  0.5544,	Adjusted R-squared:  0.4059 
plot(m3.0, SUB_2)
summary(MAG3_2A)
abline(MAG3_2A)

MAG2_2A<-lm(SUB_2~m2.0)  # Multiple R-squared:  0.8497,	Adjusted R-squared:  0.7996 
plot(m2.0, SUB_2)
summary(MAG2_2A)
abline(MAG2_2A)

MAG1_2A<-lm(SUB_2~m1.0) # Multiple R-squared:  0.7444,	Adjusted R-squared:  0.6593
plot(m1.0, SUB_2)
summary(MAG1_2A)
abline(MAG1_2A)


#### MAGNITUDE OF HYPOXIA AT 4.8 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(m4.8, SUB_2)
abline(MAG4.8_2A)
summary(MAG4.8_2A)
mylabelMAG4.8_2A<- bquote(italic(R)^2 == .(format(0.9591 , digits = 4)))
text(x = 4, y = 8, labels = mylabelMAG4.8_2A)


#############        EQUATION       ###################################
#######################################################################

par(mfrow=c(3,2))

EQ4_2A<-lm(SUB_2~EQ4.8) #Multiple R-squared:  0.8483,	Adjusted R-squared:  0.7978 
plot(EQ4.8, SUB_2)
summary(EQ4_2A)
abline(EQ4_2A)

EQ3.5_2A<-lm(SUB_2~EQ3.5)  #Multiple R-squared:  0.8766,	Adjusted R-squared:  0.8354 
plot(EQ3.5, SUB_2)
summary(EQ3.5_2A)
abline(EQ3.5_2A)

EQ3_2A<-lm(SUB_2~EQ3.0) #Multiple R-squared:  0.8752,	Adjusted R-squared:  0.8336 
plot(EQ3.0, SUB_2)
summary(EQ3_2A)
abline(EQ3_2A)

EQ2_2A<-lm(SUB_2~EQ2.0) #Multiple R-squared:  0.7988,	Adjusted R-squared:  0.7317  
plot(EQ2.0, SUB_2)
summary(EQ2_2A)
abline(EQ2_2A)

EQ1_2A<-lm(SUB_2~EQ1.0) #Multiple R-squared:  0.6952,	Adjusted R-squared:  0.5936
plot(EQ1.0, SUB_2)
summary(EQ1_2A)
abline(EQ1_2A)


#### EQUATION  OF HYPOXIA AT 3.5 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(EQ3.5, SUB_2)
abline(EQ3.5_2A)
summary(EQ3.5_2A)
mylabelEQ3.5_2A<- bquote(italic(R)^2 == .(format(0.8766 , digits = 4)))
text(x = 20, y = 5, labels = mylabelEQ3.5_2A)

#############       finaL FIGURE ONE HOUR STRONGEST REGRESSIONS  ###################################
####################################################################################################


par(mfrow=c(3,2))

#PERCENT TIME 

plot(PERC3.5, SUB_2)
abline(p3.5_2A)
summary(p3.5_2A)
mylabelp3.5_2A <- bquote(italic(R)^2 == .(format(0.9582, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_2A)

#  EQUATION 

plot(EQ3.5, SUB_2)
abline(EQ3.5_2A)
summary(EQ3.5_2A)
mylabelEQ3.5_2A<- bquote(italic(R)^2 == .(format(0.8766 , digits = 4)))
text(x = 20, y = 5, labels = mylabelEQ3.5_2A)

# MAGNITUDE

plot(m4.8, SUB_2)
abline(MAG4.8_2A)
summary(MAG4.8_2A)
mylabelMAG4.8_2A<- bquote(italic(R)^2 == .(format(0.9591 , digits = 4)))
text(x = 4, y = 6, labels = mylabelMAG4.8_2A)

# DURATION

plot(d3.0, SUB_2)
abline(DUR3_2A)
summary(DUR3_2A)
mylabelDUR3_2A <- bquote(italic(R)^2 == .(format(0.9793, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR3_2A)

# COUNT

plot(f3.5, SUB_2)
abline(c3.5_2A)
summary(c3.5_2A)
mylabelDUR2_2A <- bquote(italic(R)^2 == .(format(0.9307, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_2A)


######################################################################################################################################## #
######################################################################################################################################## #
# THREE HOUR
######################################################################################################################################## # 
#########################################################################################################################################


par(mfrow=c(3,2))
p4_3A<-lm(SUB_3~PERC4.8) #Multiple R-squared:  0.9022,	Adjusted R-squared:  0.8695 
summary(p4_3A)
plot(PERC4.8, SUB_3)
abline(p4_3A)

p3.5_3A<-lm(SUB_3~PERC3.5)  # Multiple R-squared:  0.8814,	Adjusted R-squared:  0.8419 
summary(p3.5_3A)
plot(PERC3.5, SUB_3)
abline(p3.5_3A)

p3_3A<-lm(SUB_3~PERC3.0) # Multiple R-squared:  0.8672,	Adjusted R-squared:  0.8229 
summary(p3_3A)
plot(PERC3.0, SUB_3)
abline(p3_3A)

p2_3A<-lm(SUB_3~PERC2.0) #Multiple R-squared:  0.6939,	Adjusted R-squared:  0.5919 
summary(p2_3A)
plot(PERC2.0, SUB_3)
abline(p2_3A)

p1_3A<-lm(SUB_3~PERC1.0) #Multiple R-squared:  0.5704,	Adjusted R-squared:  0.4271
summary(p1_3A)
plot(PERC1.0, SUB_3)
abline(p1_3A)


#####  PERCENTTIME OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(PERC4.8, SUB_3)
abline(p4_3A)
summary(p4_3A)
mylabelp4_3A <- bquote(italic(R)^2 == .(format(0.9022, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_3A)


############################################################################## #
############################################################################## #
############# TEST ON THE FREQUENCY OR COUNT OF DATA  ######################## #
###############################################################################


par(mfrow=c(3,2))
c4_3A<-lm(SUB_3~f4.8) #Multiple R-squared:  0.7906,	Adjusted R-squared:  0.7209 
plot(f4.8, SUB_3)
summary(c4_3A)
abline(c4_3A)

c3.5_3A<-lm(SUB_3~f3.5) #Multiple R-squared:  0.9394,	Adjusted R-squared:  0.9192 
plot(f3.5, SUB_3)
summary(c3.5_3A)
abline(c3.5_3A)

c3_3A<-lm(SUB_3~f3.0) #Multiple R-squared:  0.9254,	Adjusted R-squared:  0.9005 
plot(f3.0, SUB_3)
summary(c3_3A)
abline(c3_3A)

c2_3A<-lm(SUB_3~f2.0) #Multiple R-squared:  0.8305,	Adjusted R-squared:  0.774
plot(f2.0, SUB_3)
summary(c2_3A)
abline(c2_3A)

c1_3A<-lm(SUB_3~f1.0) #Multiple R-squared:  0.8606,	Adjusted R-squared:  0.8141 
plot(f1.0, SUB_3)
summary(c1_3A)
abline(c1_3A)

#####  Count OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(f3.5, SUB_3)
abline(c3.5_3A)
summary(c3.5_3A)
mylabelDUR2_3A <- bquote(italic(R)^2 == .(format(0.9394, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_3A)


############################################################################## #  
################     TEST OF THE DURATION DATA       ######################### #
###############################################################################
par(mfrow=c(3,2))

DUR4.8_3A<-lm(SUB_3~d4.8) #Multiple R-squared:  0.9401,	Adjusted R-squared:  0.9202 
plot(d4.8, SUB_3)
summary(DUR4.8_3A)
abline(DUR4.8_3A)

DUR3.5_3A<-lm(SUB_3~d3.5) #Multiple R-squared:  0.961,	Adjusted R-squared:  0.948
plot(d3.5, SUB_3)
summary(DUR3.5_3A)
abline(DUR3.5_3A)

DUR3_3A<-lm(SUB_3~d3.0) #Multiple R-squared:  0.9634,	Adjusted R-squared:  0.9512
plot(d3.0, SUB_3)
summary(DUR3_3A)
abline(DUR3_3A)

DUR2_3A<-lm(SUB_3~d2.0)  # Multiple R-squared:  0.8309,	Adjusted R-squared:  0.7745
plot(d2.0, SUB_3)
summary(DUR2_3A)
abline(DUR2_3A)

DUR1_3A<-lm(SUB_3~d1.0) #Multiple R-squared:  0.7369,	Adjusted R-squared:  0.6492 
plot(d1.0, SUB_3)
summary(DUR1_3A)
abline(DUR1_3A)

#####  DURAITON OF EVENTS BELOW 2.0 HAS THE STRONGEST R SQUARED OF 0.8995

par(mfrow=c(1,1))

plot(d3.0, SUB_3)
abline(DUR3_3A)
summary(DUR3_3A)
mylabelDUR3_3A <- bquote(italic(R)^2 == .(format(0.9634, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR3_3A)


##### tEST OF THE mAGNITUDE dATA###########################################

par(mfrow=c(3,2))

MAG4.8_3A<-lm(SUB_3~m4.8) #Multiple R-squared:  0.9659,	Adjusted R-squared:  0.9545 
plot(m4.8, SUB_3)
summary(MAG4.8_3A)
abline(MAG4.8_3A)

MAG3.5_3A<-lm(SUB_3~m3.5) #Multiple R-squared:  0.3799,	Adjusted R-squared:  0.1732 
plot(m3.5, SUB_3)
summary(MAG3.5_3A)
abline(MAG3.5_3A)

MAG3_3A<-lm(SUB_3~m3.0) #Multiple R-squared:  0.3808,	Adjusted R-squared:  0.1743 
plot(m3.0, SUB_3)
summary(MAG3_3A)
abline(MAG3_3A)

MAG2_3A<-lm(SUB_3~m2.0)  #Multiple R-squared:  0.6661,	Adjusted R-squared:  0.5548
plot(m2.0, SUB_3)
summary(MAG2_3A)
abline(MAG2_3A)

MAG1_3A<-lm(SUB_3~m1.0) #Multiple R-squared:  0.5657,	Adjusted R-squared:  0.4209
plot(m1.0, SUB_3)
summary(MAG1_3A)
abline(MAG1_3A)


#### MAGNITUDE OF HYPOXIA AT 4.8 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(m4.8, SUB_3)
abline(MAG4.8_3A)
summary(MAG4.8_3A)
mylabelMAG4.8_3A<- bquote(italic(R)^2 == .(format(0.9659 , digits = 4)))
text(x = 4, y = 8, labels = mylabelMAG4.8_3A)


#############        EQUATION       ################################## #
#######################################################################

par(mfrow=c(3,2))

EQ4_3A<-lm(SUB_3~EQ4.8) #Multiple R-squared:  0.8957,	Adjusted R-squared:  0.8609 
plot(EQ4.8, SUB_3)
summary(EQ4_3A)
abline(EQ4_3A)

EQ3.5_3A<-lm(SUB_3~EQ3.5)  #Multiple R-squared:  0.8633,	Adjusted R-squared:  0.8177 
plot(EQ3.5, SUB_3)
summary(EQ3.5_3A)
abline(EQ3.5_3A)

EQ3_3A<-lm(SUB_3~EQ3.0) #Multiple R-squared:  0.8654,	Adjusted R-squared:  0.8205
plot(EQ3.0, SUB_3)
summary(EQ3_3A)
abline(EQ3_3A)

EQ2_3A<-lm(SUB_3~EQ2.0) #Multiple R-squared:  0.7523,	Adjusted R-squared:  0.6698  
plot(EQ2.0, SUB_3)
summary(EQ2_3A)
abline(EQ2_3A)

EQ1_3A<-lm(SUB_3~EQ1.0) #Multiple R-squared:  0.5889,	Adjusted R-squared:  0.4519
plot(EQ1.0, SUB_3)
summary(EQ1_3A)
abline(EQ1_3A)


#### EQUATION  OF HYPOXIA AT 3.5 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(EQ4.8, SUB_3)
abline(EQ4_3A)
summary(EQ4_3A)
mylabelEQ4_3A<- bquote(italic(R)^2 == .(format(0.8957 , digits = 4)))
text(x = 20, y = 5, labels = mylabelEQ4_3A)

#############       finaL FIGURE ONE HOUR STRONGEST REGRESSIONS  ################################## #
####################################################################################################

par(mfrow=c(3,2))

#PERCENT TIME 

plot(PERC3.5, SUB_3)
abline(p3.5_3A)
summary(p3.5_3A)
mylabelp3.5_3A <- bquote(italic(R)^2 == .(format(0.8814, digits = 4)))
text(x = 10, y =5, labels = mylabelp3.5_3A)

#  EQUATION 

plot(EQ3.5, SUB_3)
abline(EQ3.5_3A)
summary(EQ3.5_3A)
mylabelEQ3.5_3A<- bquote(italic(R)^2 == .(format(0.8354 , digits = 4)))
text(x = 20, y = 5, labels = mylabelEQ3.5_3A)

# MAGNITUDE

plot(m4.8, SUB_3)
abline(MAG4.8_3A)
summary(MAG4.8_3A)
mylabelMAG4.8_3A<- bquote(italic(R)^2 == .(format(0.9545 , digits = 4)))
text(x = 4, y = 5, labels = mylabelMAG4.8_3A)

# DURATION

plot(d3.0, SUB_3)
abline(DUR3_3A)
summary(DUR3_3A)
mylabelDUR3_3A <- bquote(italic(R)^2 == .(format(0.9512, digits = 4)))
text(x = 2, y = 5, labels = mylabelDUR3_3A)

# COUNT

plot(f3.5, SUB_3)
abline(c3.5_3A)
summary(c3.5_3A)
mylabelDUR2_3A <- bquote(italic(R)^2 == .(format(0.9192, digits = 4)))
text(x = 10, y =5, labels = mylabelDUR2_3A)

######################################################################## #
######### RECOVERY 1 HR DATA ########################################### #
##########################################################################
hr_recov # VIEW DATA
# recovery 1 hr inconditions >6.0 after 24 hrs of hypoxic exposure

####TEST THE PERCENT FACTOR ON heartbeat

par(mfrow=c(3,2))
p4_recovA<-lm(hr_recov~PERC4.8) # Multiple R-squared:  0.8385,	Adjusted R-squared:  0.7846 
summary(p4_recovA)
plot(PERC4.8, hr_recov)
abline(p4_recovA)

p3.5_recovA<-lm(hr_recov~PERC3.5)  # Multiple R-squared:  0.7311,	Adjusted R-squared:  0.6415 
summary(p3.5_recovA)
plot(PERC3.5, hr_recov)
abline(p3.5_recovA)

p3_recovA<-lm(hr_recov~PERC3.0) # Multiple R-squared:  0.7511,	Adjusted R-squared:  0.6682 
summary(p3_recovA)
plot(PERC3.0, hr_recov)
abline(p3_recovA)

p2_recovA<-lm(hr_recov~PERC2.0) # Multiple R-squared:  0.5429,	Adjusted R-squared:  0.3905 
summary(p2_recovA)
plot(PERC2.0, hr_recov)
abline(p2_recovA)

p1_recovA<-lm(hr_recov~PERC1.0) # Multiple R-squared:  0.4123,	Adjusted R-squared:  0.2164 
summary(p1_recovA)
plot(PERC1.0, hr_recov)
abline(p1_recovA)


#####  PERCENTTIME OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(PERC4.8, hr_recov)
abline(p4_recovA)
summary(p4_recovA)
mylabelp4_recovA<- bquote(italic(R)^2 == .(format(0.8385, digits = 4)))
text(x = 8, y =2, labels = mylabelp4_recovA)


############################################################################## #
############################################################################## #
############# TEST ON THE FREQUENCY OR COUNT OF DATA  ######################## #
###############################################################################

par(mfrow=c(3,2))
c4_recovA<-lm(hr_recov~f4.8) #Multiple R-squared:  0.7584,	Adjusted R-squared:  0.6779 
plot(f4.8, hr_recov)
summary(c4_recovA)
abline(c4_recovA)

c3.5_recovA<-lm(hr_recov~f3.5) #Multiple R-squared:  0.8393,	Adjusted R-squared:  0.7857 
plot(f3.5, hr_recov)
summary(c3.5_recovA)
abline(c3.5_recovA)

c3_recovA<-lm(hr_recov~f3.0) #Multiple R-squared:  0.8283,	Adjusted R-squared:  0.771 
plot(f3.0, hr_recov)
summary(c3_recovA)
abline(c3_recovA)

c2_recovA<-lm(hr_recov~f2.0) #Multiple R-squared:  0.6986,	Adjusted R-squared:  0.5981 
plot(f2.0, hr_recov)
summary(c2_recovA)
abline(c2_recovA)

c1_recovA<-lm(hr_recov~f1.0) #Multiple R-squared:  0.7368,	Adjusted R-squared:  0.6491 
plot(f1.0, hr_recov)
summary(c1_recovA)
abline(c1_recovA)

#####  Count OF EVENTS BELOW 3.5 HAS THE STRONGEST R SQUARED OF 0.6188

par(mfrow=c(1,1))

plot(f3.5, hr_recov)
abline(c3.5_recovA)
summary(c3.5_recovA)
mylabelDUR2_recovA <- bquote(italic(R)^2 == .(format(0.8393, digits = 4)))
text(x = 5, y =2, labels = mylabelDUR2_recovA)


############################################################################## #
################     TEST OF THE DURATION DATA       ######################### #
###############################################################################

par(mfrow=c(3,2))

DUR4.8_recovA<-lm(hr_recov~d4.8) #Multiple R-squared:  0.8882,	Adjusted R-squared:  0.8509 
plot(d4.8, hr_recov)
summary(DUR4.8_recovA)
abline(DUR4.8_recovA)

DUR3.5_recovA<-lm(hr_recov~d3.5) #Multiple R-squared:  0.8487,	Adjusted R-squared:  0.7983 
plot(d3.5, hr_recov)
summary(DUR3.5_recovA)
abline(DUR3.5_recovA)

DUR3_recovA<-lm(hr_recov~d3.0) #Multiple R-squared:  0.8503,	Adjusted R-squared:  0.8004 
plot(d3.0, hr_recov)
summary(DUR3_recovA)
abline(DUR3_recovA)

DUR2_recovA<-lm(hr_recov~d2.0)  # Multiple R-squared:  0.6788,	Adjusted R-squared:  0.5718 
plot(d2.0, hr_recov)
summary(DUR2_recovA)
abline(DUR2_recovA)

DUR1_recovA<-lm(hr_recov~d1.0) #Multiple R-squared:  0.5486,	Adjusted R-squared:  0.3981 
plot(d1.0, hr_recov)
summary(DUR1_recovA)
abline(DUR1_recovA)

#####  DURAITON OF EVENTS BELOW 4.8 HAS THE STRONGEST R SQUARED OF 0.8882

par(mfrow=c(1,1))

plot(d4.8, hr_recov)
abline(DUR4.8_recovA)
summary(DUR4.8_recovA)
mylabelDUR4.8_recovA<- bquote(italic(R)^2 == .(format(0.8882, digits = 4)))
text(x = 5, y = 2, labels = mylabelDUR4.8_recovA)


##### tEST OF THE mAGNITUDE dATA########################################### #
#############################################################################

par(mfrow=c(3,2))

MAG4.8_recovA<-lm(hr_recov~m4.8) #Multiple R-squared:  0.9659,	Adjusted R-squared:  0.9545 
plot(m4.8, hr_recov)
summary(MAG4.8_recovA)
abline(MAG4.8_recovA)

MAG3.5_recovA<-lm(hr_recov~m3.5) #Multiple R-squared:  0.3799,	Adjusted R-squared:  0.1732 
plot(m3.5, hr_recov)
summary(MAG3.5_recovA)
abline(MAG3.5_recovA)

MAG3_recovA<-lm(hr_recov~m3.0) #Multiple R-squared:  0.3808,	Adjusted R-squared:  0.1743 
plot(m3.0, hr_recov)
summary(MAG3_recovA)
abline(MAG3_recovA)

MAG2_recovA<-lm(hr_recov~m2.0)  #Multiple R-squared:  0.6661,	Adjusted R-squared:  0.5548
plot(m2.0, hr_recov)
summary(MAG2_recovA)
abline(MAG2_recovA)

MAG1_recovA<-lm(hr_recov~m1.0) #Multiple R-squared:  0.5657,	Adjusted R-squared:  0.4209
plot(m1.0, hr_recov)
summary(MAG1_recovA)
abline(MAG1_recovA)


#### MAGNITUDE OF HYPOXIA AT 4.8 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(m4.8, hr_recov)
abline(MAG4.8_recovA)
summary(MAG4.8_recovA)
mylabelMAG4.8_recovA<- bquote(italic(R)^2 == .(format(0.9545 , digits = 4)))
text(x = 4, y = 8, labels = mylabelMAG4.8_recovA)


#############        EQUATION       ################################## #
#######################################################################

par(mfrow=c(3,2))

EQ4_recovA<-lm(hr_recov~EQ4.8) #Multiple R-squared:  0.8957,	Adjusted R-squared:  0.8609 
plot(EQ4.8, hr_recov)
summary(EQ4_recovA)
abline(EQ4_recovA)

EQ3.5_recovA<-lm(hr_recov~EQ3.5)  #Multiple R-squared:  0.8633,	Adjusted R-squared:  0.8177 
plot(EQ3.5, hr_recov)
summary(EQ3.5_recovA)
abline(EQ3.5_recovA)

EQ3_recovA<-lm(hr_recov~EQ3.0) #Multiple R-squared:  0.8654,	Adjusted R-squared:  0.8205
plot(EQ3.0, hr_recov)
summary(EQ3_recovA)
abline(EQ3_recovA)

EQ2_recovA<-lm(hr_recov~EQ2.0) #Multiple R-squared:  0.7523,	Adjusted R-squared:  0.6698  
plot(EQ2.0, hr_recov)
summary(EQ2_recovA)
abline(EQ2_recovA)

EQ1_recovA<-lm(hr_recov~EQ1.0) #Multiple R-squared:  0.5889,	Adjusted R-squared:  0.4519
plot(EQ1.0, hr_recov)
summary(EQ1_recovA)
abline(EQ1_recovA)


#### EQUATION  OF HYPOXIA AT 3.5 HAS STRONGEST R SQUARED OF 0.8686 
par(mfrow=c(1,1))

plot(EQ4.8, hr_recov)
abline(EQ4_recovA)
summary(EQ4_recovA)
mylabelEQ4_recovA<- bquote(italic(R)^2 == .(format(0.8609 , digits = 4)))
text(x = 20, y = 5, labels = mylabelEQ4_recovA)

#############       finaL FIGURE ONE HOUR STRONGEST REGRESSIONS  ################################## #
####################################################################################################

par(mfrow=c(3,2))

#PERCENT TIME 

plot(PERC4.8, hr_recov, xlim=c(0, 30), ylim=c(-12, 16))
abline(p4_recovA)
summary(p4_recovA)
mylabelp4_recovA<- bquote(italic(R)^2 == .(format(0.8385, digits = 4)))
text(x = 8, y =2, labels = mylabelp4_recovA)

#  EQUATION 

#plot(EQ3.5, hr_recov)
#abline(EQ3.5_recovA)
#summary(EQ3.5_recovA)
#mylabelEQ3.5_recovA<- bquote(italic(R)^2 == .(format(0.8354 , digits = 4)))
#text(x = 20, y = 5, labels = mylabelEQ3.5_recovA)

# MAGNITUDE

#plot(m4.8, hr_recov)
#abline(MAG4.8_recovA)
#summary(MAG4.8_recovA)
#mylabelMAG4.8_recovA<- bquote(italic(R)^2 == .(format(0.9545 , digits = 4)))
#text(x = 4, y = 5, labels = mylabelMAG4.8_recovA)

# DURATION

plot(d4.8, hr_recov)
abline(DUR4.8_recovA)
summary(DUR4.8_recovA)
mylabelDUR4.8_recovA<- bquote(italic(R)^2 == .(format(0.8882, digits = 4)))
text(x = 5, y = 2, labels = mylabelDUR4.8_recovA)

# COUNT

plot(f3.5, hr_recov)
abline(c3.5_recovA)
summary(c3.5_recovA)
mylabelDUR2_recovA <- bquote(italic(R)^2 == .(format(0.8393, digits = 4)))
text(x = 5, y =2, labels = mylabelDUR2_recovA)


