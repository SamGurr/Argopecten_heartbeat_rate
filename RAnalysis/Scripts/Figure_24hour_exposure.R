# script written by Sam J Gurr
# Last edit: 20200604
# 24 -hour hb plot

rm(list=ls()) # remove all data - start fresh

# Install packages if not already in your library-----------------------------------------------------------------------------------------------
require(dplyr) 
require(ggplot2) 
require(tidyr) 
# Load packages and pacage version/date/import/depends info
library(dplyr) 
library(ggplot2) 
library(tidyr) 

#set working directory -------------------------------------------------------------- #
setwd("C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/")

# load datasets
HB.24hr.hypoxia <- read.csv(file="Data/Post_insitu_experiment/Heartbeat_data/SUMMARY/24-hour_hb_response/Summary_24_hour_exposure.csv", header=T) # full 24 hour mean data (every 10 minutes)  of O2 <= 2.0 mg L-1
HB.24hr.hypoxia # view the data
# create pivot table from lide to longer with tidyr
Pivot_Table_24hour <- HB.24hr.hypoxia %>%  
  pivot_longer(
    cols = starts_with("hr"),
    names_to = "time",
    values_to = "heartrates",
    values_drop_na = TRUE
  )
# substr the hr_# to call just the hour number - we will use this for the x axis labels downstream in the figure
Pivot_Table_24hour$hour <- as.numeric(substr(Pivot_Table_24hour$time, 4,5))
#
Plot <- Pivot_Table_24hour %>% 
  group_by(Site,treatment,hour)%>%
  summarise(mean= mean(heartrates), sd= sd(heartrates)) %>% 
  ggplot(aes(x=hour, 
             y=mean, 
             group=treatment)) +
  geom_line() +
  theme_classic() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_point(aes(shape=treatment, fill=treatment), size = 5)+
  scale_shape_manual(values=c(21,19))+
  scale_fill_manual(values=c("white","black")) +
  facet_wrap(~ Site) +
  #facet_grid(Site ~ .) +
  scale_x_continuous(expand = c(0,0)) +
  labs( y = "heartbeat rate (BPM)",
       x = "time (hours)") +
  theme(legend.position="none")
print(Plot)

ggsave(Plot, file="C:/Users/samjg/Documents/My_Projects/Argopecten_hearbeat_rate/Argopecten_heartbeat_rate/RAnalysis/Output/24_hour_exposure_plots.pdf",
       width=35, height=20, units = "cm", dpi=500)
  
