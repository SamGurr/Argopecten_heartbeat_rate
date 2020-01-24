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
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/RAnalysis/")
# laod table
DO_TEMP.data <- read.csv(file="Data/final_all_sites_DO_no_NA.csv", header=T) #read Size.info data
DO_TEMP.data # view data
unique(DO_TEMP.data$SITE)
DO_TEMP.data <- DO_TEMP.data %>% dplyr::filter(SITE %in% c("B", "FI", "M", "N", "Q", "X"))
unique(DO_TEMP.data$SITE) # sites used

weekly_data <- read.csv(file="Data/LIMMN_2016_data/Weekly_data/chlorophyl_secchi_sal_minmaxDO.csv", header=T) #read Size.info data
weekly_data <- weekly_data %>% dplyr::filter(LIMMN_site_initals %in% c("B", "F", "M", "N", "Q", "X"))
names(weekly_data) # view column names


# WATER CHARACTERISTICS TABLE 1

# TEMPERATURE
TEMP.table <- DO_TEMP.data %>% 
  dplyr::select(SITE, TEMP) %>% 
  dplyr::group_by(SITE) %>% # call column to summarize 
  dplyr::summarise_each(funs(mean,sd))
TEMP.table

#SALINITY
salinity.table <- ddply(weekly_data, c("LIMMN_site_initals"), summarise,
               N    = length(Salinity),
               mean = mean(Salinity),
               sd   = sd(Salinity),
               se   = sd / sqrt(N))
salinity.table

#SECCHI DEPTH
weekly_data.secchi_depth.OM <- weekly_data %>% dplyr::select(c(LIMMN_site_initals, secchi_depth))
weekly_data.secchi_depth.OM <- na.omit(weekly_data.secchi_depth.OM)
secchi_depth.table <- ddply(weekly_data.secchi_depth.OM, c("LIMMN_site_initals"), summarise,
                        N    = length(secchi_depth),
                        mean = mean(secchi_depth),
                        sd   = sd(secchi_depth),
                        se   = sd / sqrt(N))
secchi_depth.table

#CHLOROPHYLL
weekly_data.CHL.OM <- weekly_data %>% dplyr::select(c(LIMMN_site_initals, chl_a_ug_l))
weekly_data.CHL.OM <- na.omit(weekly_data.CHL.OM)
chl_a_ug_l.table <- ddply(weekly_data.CHL.OM, c("LIMMN_site_initals"), summarise,
                            N    = length(chl_a_ug_l),
                            mean = mean(chl_a_ug_l),
                            sd   = sd(chl_a_ug_l),
                            se   = sd / sqrt(N))
chl_a_ug_l.table

