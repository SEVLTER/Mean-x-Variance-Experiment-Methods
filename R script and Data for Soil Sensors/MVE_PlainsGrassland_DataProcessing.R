###### DATA PROCESSING FOR MVE Blue SOIL SENSORS ######
###### Date updated: Date updated: 20 Nov 2022 J. Rudgers ######

rm(list=ls(all=TRUE)) #give R a blank slate
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)

####### QUICK META DATA ######################################
# VWC UNITS are m^3/m^3
# T UNITS are degrees C
# TIMESTAMP is "%Y-%m-%d %H:%M:%S",tz = "UTC"
# sensor manual and information lives here: http://publications.metergroup.com/Manuals/20424_5TM_Manual_Web.pdf

###### IMPORT AND RESHAPE DATA ##################################################
# download from google drive
# https://drive.google.com/drive/folders/1WsS7mtBol_m050OEONeuGqvDy7OOzi87?usp=sharing
# eventually we can have a script that fetches directly 

# set working directory
setwd("C:/Users/jrudgers/Desktop/SEV Data/Mean-Var Experiment/blue grama mean variance/soil moist temp MVE blue/")

###### import data as DAT file
#quickly look at first few lines of data
readLines("MVE_Blue4oct22.dat", n=10)

#read in data file
MVE.dat<-read.delim("MVE_Blue4oct22.dat", sep=",",header=TRUE, skip=1, as.is=TRUE)
head(MVE.dat)
# on 4 oct 2022: 68597 obs 110 variables

#make sure to only run the next line one time or else you will be removing lines of data
#rows 1-2 have logger info that was deleted during data file import
MVE.dat<-MVE.dat[-(1:2),]
head(MVE.dat)

# **current file has two rows 6575-6, RECORD 5055 that do not import properly
# **missing a comma, so gets smushed into a cell in prior time point
# **because data were not recorded between 2/14/19 and 2/18/19, the timestamp is not a date
# **this deletes those two rows after looking at them
MVE.dat[6574,]
MVE.dat[6575,]
MVE.data<-MVE.dat[-c(6574,6575),]
# can write out to check it
# write.csv(MVE.data,"MVE_Plains_sensors.csv")

# make a timeline variable for x-axis of plot by converting to date data
MVE.data$TIMESTAMP<-as.POSIXct(MVE.data$TIMESTAMP,"%Y-%m-%d %H:%M:%OS",tz = "UTC")
summary(MVE.data$TIMESTAMP)
# convert columns to numeric data
MVE.data<-data.frame(MVE.data %>% mutate_if(is.character,as.numeric))
str(MVE.data)

# before changing shape of the file, remove the column called RECORD
# ONLY RUN THIS ONCE, OR ELSE YOU WILL DELETE ROWS
head(MVE.data)
MVE.data<-MVE.data[,-2]

###### Reshape data
# melt it so each column becomes a variable and measurements are all in just one column
MVE.melt<-melt(MVE.data,id="TIMESTAMP",variable_name="sensor")
head(MVE.melt)

# create a sensor.id column
MVE.melt$sensor.id<-MVE.melt$variable

##### Separate sensor column into individual labels
# this can take some processing time
MVE.data<-MVE.melt %>% tidyr::separate(variable,c("sensor","plot","depth"),sep="_") %>% 
dplyr::mutate(sensor=as.factor(sensor),plot=as.factor(plot),depth.f=as.factor(depth),depth=as.numeric(depth)) 
head(MVE.data)

# check that new sensor information columns have correct levels
levels(MVE.data$sensor) # T=temperature of soil (C), #VWC = volumetric water content (%)
levels(MVE.data$plot) #should be a factor: "P1" "P10" "P11" "P12" "P13" "P14" "P15" "P16" "P17" "P18" "P2"  "P3"  "P4"  "P5"  "P6"  "P7"  "P8"  "P9" 
levels(MVE.data$depth.f) #12, 22, or 37 cm deep, should be factor
summary(MVE.data$depth) #same levels but as integer

# extract components and make date values numeric
MVE.data$year<-as.numeric(format(MVE.data$TIMESTAMP,"%Y"))
MVE.data$month<-as.numeric(format(MVE.data$TIMESTAMP,"%m"))
MVE.data$day<-as.numeric(format(MVE.data$TIMESTAMP,"%d"))
MVE.data$hour<-hour(as.POSIXct(MVE.data$TIMESTAMP,"%Y-%m-%d %H:%M:%OS",tz = "UTC"))
MVE.data$minute<-as.numeric(format(MVE.data$TIMESTAMP,"%M"))

# check it
max(MVE.data$hour) #23
head(MVE.data)

##### Merge SENSOR DATA with SENSOR LABELS
MVE.label<-read.csv("MVE_PlainsGrassland_5TM_Sensor_Labels.csv",stringsAsFactors = T)
# var_2019 ... gives level of watering for variance treatment
# trt_2019 ... give total water treatment in each year, levels: -75, -50, -25, 0 +25, +50)
MVE.label$mean_trt<-recode_factor(as.factor(MVE.label$mean_trt),"-25"="drier","0"="ambient")
summary(MVE.label)
MVE.label$block<-as.factor(MVE.label$block)
MVE<-merge(MVE.data,MVE.label,by.x="sensor.id",by.y="sensor.id")
# should have same number of observations as the MVE.data
str(MVE)

# convert treatments to factors (categories)
MVE$mean.f<-as.factor(MVE$mean_trt)
MVE$var.f<-as.factor(MVE$var_trt)


# check that factors are correct
levels(MVE$mean.f)
levels(MVE$var.f)

# make date values factors to enable summaries
MVE$year.f<-as.factor(MVE$year)
MVE$month.f<-as.factor(MVE$month)
MVE$day.f<-as.factor(MVE$day)

# make sure value column is numeric, check on NAs
summary(MVE$value)
# 230412 /7408044 is missing data at 3.1% due to sensor failures 
str(MVE)
max(MVE$TIMESTAMP)
#hour-min-sec of timestamp not preserved if you don't convert timestamp to character data prior to export
MVE$TIMESTAMP<-as.character(MVE$TIMESTAMP)
max(MVE$TIMESTAMP)
# write out processed data file (note: it will be large)
write.csv(MVE,"MVE_PlainsGrassland_SoilMoistureTemperature.csv")


# reimport data file, if needed
#MVE.import<-read.csv("MVE_PlainsGrassland_SoilMoistureTemperature.csv")
#MVE.import$TIMESTAMP<-as.POSIXct(MVE.import$TIMESTAMP,"%Y-%m-%d %H:%M:%OS",tz = "UTC")
#summary(MVE.import$TIMESTAMP)




