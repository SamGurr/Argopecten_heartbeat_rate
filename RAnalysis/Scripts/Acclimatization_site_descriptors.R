#Title: Insitu tempeature 
#Author:Sam Gurr 
#Edited by: Sam Gurr
#Date Last Modified: 20200123
#See Readme file for details


rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('doBY')
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
#Read in required libraries
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(ggpubr)         # Version: 0.1.8 Date: 2018-08-30, Depends: R (>= 3.1.0), ggplot2, magrittrImports: ggrepel, grid, ggsci, stats, utils, tidyr, purrr, dplyr(>=0.7.1), cowplot, ggsignif, scales, gridExtra, glue, polynom
library(Rmisc)          # Version: 1.5 Packaged: 2013-10-21, Depends: lattice, plyr
library(plotrix)        # Version: 3.7-4, Date/Publication: 2018-10-03
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(gridExtra)      # Version: 2.3, Date/Publication: 2017-09-09, Imports: gtable, grid, grDevices, graphics, utils
library(reshape)        # Version: 0.8.7, Date/Publication: 2017-08-06, Depends: R (>= 2.6.1) Imports: plyr
library(multcompView)   # Version: 0.1-7, Date/Publication: 2015-07-31, Imports: grid
library(Rmisc)
library(lmtest)
library(car)
library(ggpubr)

# Set Working Directory:
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")
# load spreadsheets
# all raw insitu HOBO logger data
DO_TEMP.data <- read.csv(file="Data/LIMMN_2016_data/HOBO_logger_data/DO_temp_raw_insitu.csv", header=T) #read Size.info data
DO_TEMP.data # view data
unique(DO_TEMP.data$SITE)
DO_TEMP.data <- DO_TEMP.data %>% dplyr::filter(SITE %in% c("B", "FI", "M", "N", "Q", "X"))
unique(DO_TEMP.data$SITE) # sites used
# weekly site data for chlorophyll, secchi depth, salinity
weekly_data <- read.csv(file="Data/LIMMN_2016_data/Weekly_data/chlorophyl_secchi_sal_minmaxDO.csv", header=T) #read Size.info data
weekly_data <- weekly_data %>% dplyr::filter(LIMMN_site_initals %in% c("B", "F", "M", "N", "Q", "X"))
names(weekly_data) # view column names
# DO change data 
DO_change <- read.csv(file="Data/LIMMN_2016_data/HOBO_logger_data/DO_change.csv", header=T) #read Size.info data
DO_change # view data

######################################### #
############ TABLE 1 #################### #
###########################################

# WATER CHARACTERISTICS TABLE 1
# TEMPERATURE
TEMP.table <- DO_TEMP.data %>% 
  dplyr::select(SITE, TEMP) %>% 
  dplyr::group_by(SITE) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,sd))
TEMP.table # view table
#SALINITY
salinity.table <- ddply(weekly_data, c("LIMMN_site_initals"), summarise,
               N    = length(Salinity),
               mean = mean(Salinity),
               sd   = sd(Salinity),
               se   = sd / sqrt(N))
salinity.table  # view table
#SECCHI DEPTH
weekly_data.secchi_depth.OM <- weekly_data %>% dplyr::select(c(LIMMN_site_initals, secchi_depth)) # select targets
weekly_data.secchi_depth.OM <- na.omit(weekly_data.secchi_depth.OM) # ommit NAs
secchi_depth.table <- ddply(weekly_data.secchi_depth.OM, c("LIMMN_site_initals"), summarise,
                        N    = length(secchi_depth),
                        mean = mean(secchi_depth),
                        sd   = sd(secchi_depth),
                        se   = sd / sqrt(N))
secchi_depth.table  # view table

#CHLOROPHYLL
weekly_data.CHL.OM <- weekly_data %>% dplyr::select(c(LIMMN_site_initals, chl_a_ug_l)) # select targets
weekly_data.CHL.OM <- na.omit(weekly_data.CHL.OM)  # ommit NAs
chl_a_ug_l.table <- ddply(weekly_data.CHL.OM, c("LIMMN_site_initals"), summarise,
                            N    = length(chl_a_ug_l),
                            mean = mean(chl_a_ug_l),
                            sd   = sd(chl_a_ug_l),
                            se   = sd / sqrt(N))
chl_a_ug_l.table  # view table

# DO Dynamics during in-situ acclimatization
total_rows.DO_change <- nrow(DO_change) # 1077 total rows of hourly change data (1077/24 == ~45 days)
total_rows.DO_RAW <- nrow(DO_TEMP.data) # 40169 total rows of raw data

# Fire Island
# NOTE: Percent time of severe and moderate hypoxia in Table 1 was calculated in the R script HypoxiaDetectionWithOutputs.R
# view this data the the 'Outputs' folder

FI_DO_change <- DO_change %>% dplyr::select(Date_Time, FI) # select date and target site
# rate of change higher than 0.5
FI_DO_change.greater_0.5 <- FI_DO_change %>% dplyr::filter(FI > 0.5) # greater than 0.5
FI_DO_change.lower_neg0.5 <- FI_DO_change %>% dplyr::filter(FI < -0.5) # less than 0.5 - sum of both are rates greater than 0.5 whether decreasing or increasing
((nrow(FI_DO_change.greater_0.5) + nrow(FI_DO_change.lower_neg0.5))/total_rows.DO_change)*100 # 6.220984 % time > 0.5 mg L hr
# DO increasde greater than 1
FI_DO_change.INCREASE <- FI_DO_change %>% dplyr::filter(FI > 1) # greater than 1 DO increase
(nrow(FI_DO_change.INCREASE)/total_rows.DO_change)*100 # 0.2785515 % time
# DO decline greater than -1
FI_DO_change.DECLINE <- FI_DO_change %>% dplyr::filter(FI < -1) # less than 1 DO decline
(nrow(FI_DO_change.DECLINE)/total_rows.DO_change)*100 # 0 % time
# minimum rate change
min(FI_DO_change$FI) # minimum rate change -0.8866667
# maximum rate change
max(FI_DO_change$FI) # max rate change 1.131667

# Sag Harbor
# NOTE: Percent time of severe and moderate hypoxia in Table 1 was calculated in the R script HypoxiaDetectionWithOutputs.R
# view this data the the 'Outputs' folder

Sag_DO_change <- DO_change %>% dplyr::select(Date_Time, SAG) # select date and target site
# rate of change higher than 0.5
Sag_DO_change.greater_0.5 <- Sag_DO_change %>% dplyr::filter(SAG > 0.5) # greater than 0.5
Sag_DO_change.lower_neg0.5 <- Sag_DO_change %>% dplyr::filter(SAG < -0.5) # less than 0.5 - sum of both are rates greater than 0.5 whether decreasing or increasing
((nrow(Sag_DO_change.greater_0.5) + nrow(Sag_DO_change.lower_neg0.5))/total_rows.DO_change)*100 # 25.71959 % time > 0.5 mg L hr
# DO increasde greater than 1
Sag_DO_change.INCREASE <- Sag_DO_change %>% dplyr::filter(SAG > 1) # greater than 1 DO increase
(nrow(Sag_DO_change.INCREASE)/total_rows.DO_change)*100 # 3.249768 % time
# DO decline greater than -1
Sag_DO_change.DECLINE <- Sag_DO_change %>% dplyr::filter(SAG < -1) # less than 1 DO decline
(nrow(Sag_DO_change.DECLINE)/total_rows.DO_change)*100 # 2.135562 % time
# minimum rate change
min(Sag_DO_change$SAG) # minimum rate change -1.8
# maximum rate change
max(Sag_DO_change$SAG) # max rate change  2.231667

# Quantuck
# NOTE: Percent time of severe and moderate hypoxia in Table 1 was calculated in the R script HypoxiaDetectionWithOutputs.R
# view this data the the 'Outputs' folder

Q_DO_change <- DO_change %>% dplyr::select(Date_Time, QTK ) # select date and target site
# rate of change higher than 0.5
Q_DO_change.greater_0.5 <- Q_DO_change %>% dplyr::filter(QTK > 0.5) # greater than 0.5
Q_DO_change.lower_neg0.5 <- Q_DO_change %>% dplyr::filter(QTK < -0.5) # less than 0.5 - sum of both are rates greater than 0.5 whether decreasing or increasing
((nrow(Q_DO_change.greater_0.5) + nrow(Q_DO_change.lower_neg0.5))/total_rows.DO_change)*100 # 59.61003 % time > 0.5 mg L hr
# DO increasde greater than 1
Q_DO_change.INCREASE <- Q_DO_change %>% dplyr::filter(QTK > 1) # greater than 1 DO increase
(nrow(Q_DO_change.INCREASE)/total_rows.DO_change)*100 # 15.32033 % time
# DO decline greater than -1
Q_DO_change.DECLINE <- Q_DO_change %>% dplyr::filter(QTK < -1) #  less than 1 DO decline
(nrow(Q_DO_change.DECLINE)/total_rows.DO_change)*100 # 14.94893 % time
# minimum rate change
min(Q_DO_change$Q) # minimum rate change 3.881667
# maximum rate change
max(Q_DO_change$Q) # max rate change  3.881667

# Nicoll Bay
# NOTE: Percent time of severe and moderate hypoxia in Table 1 was calculated in the R script HypoxiaDetectionWithOutputs.R
# view this data the the 'Outputs' folder

NCB_DO_change <- DO_change %>% dplyr::select(Date_Time, NCB) # select date and target site
# rate of change higher than 0.5
NCB_DO_change.greater_0.5 <- NCB_DO_change %>% dplyr::filter(NCB > 0.5) # greater than 0.5
NCB_DO_change.lower_neg0.5 <- NCB_DO_change %>% dplyr::filter(NCB < -0.5) # less than 0.5 - sum of both are rates greater than 0.5 whether decreasing or increasing
((nrow(NCB_DO_change.greater_0.5) + nrow(NCB_DO_change.lower_neg0.5))/total_rows.DO_change)*100 # 48.09656 % time > 0.5 mg L hr
# DO increasde greater than 1
NCB_DO_change.INCREASE <- NCB_DO_change %>% dplyr::filter(NCB > 1) # greater than 1 DO increase
(nrow(NCB_DO_change.INCREASE)/total_rows.DO_change)*100 # 10.67781 % time
# DO decline greater than -1
NCB_DO_change.DECLINE <- NCB_DO_change %>% dplyr::filter(NCB < -1) #  less than 1 DO decline
(nrow(NCB_DO_change.DECLINE)/total_rows.DO_change)*100 # 9.0065 % time
# minimum rate change
min(NCB_DO_change$NCB) # minimum rate change 3.886667
# maximum rate change
max(NCB_DO_change$NCB) # max rate change  4.365

# Moneybogue
# NOTE: Percent time of severe and moderate hypoxia in Table 1 was calculated in the R script HypoxiaDetectionWithOutputs.R
# view this data the the 'Outputs' folder

MBC_DO_change <- DO_change %>% dplyr::select(Date_Time, MBC) # select date and target site
# rate of change higher than 0.5
MBC_DO_change.greater_0.5 <- MBC_DO_change %>% dplyr::filter(MBC> 0.5) # greater than 0.5
MBC_DO_change.lower_neg0.5 <- MBC_DO_change %>% dplyr::filter(MBC< -0.5) # less than 0.5 - sum of both are rates greater than 0.5 whether decreasing or increasing
((nrow(MBC_DO_change.greater_0.5) + nrow(MBC_DO_change.lower_neg0.5))/total_rows.DO_change)*100 # 55.24605 % time > 0.5 mg L hr
# DO increasde greater than 1
MBC_DO_change.INCREASE <- MBC_DO_change %>% dplyr::filter(MBC> 1) # greater than 1 DO increase
(nrow(MBC_DO_change.INCREASE)/total_rows.DO_change)*100 # 13.74188 % time
# DO decline greater than -1
MBC_DO_change.DECLINE <- MBC_DO_change %>% dplyr::filter(MBC< -1) #  less than 1 DO decline
(nrow(MBC_DO_change.DECLINE)/total_rows.DO_change)*100 # 3.926667 % time
# minimum rate change
min(MBC_DO_change$MBC) # minimum rate change 3.926667
# maximum rate change
max(MBC_DO_change$MBC) # max rate change  4.463333
