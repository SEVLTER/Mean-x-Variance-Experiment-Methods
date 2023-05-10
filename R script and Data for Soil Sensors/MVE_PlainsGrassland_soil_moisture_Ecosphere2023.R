###### GRAPHICS AND ANALYSIS SCRIPT FOR MVE SOIL SENSORS ######
###### Date updated: 10 January 2023 by J. Rudgers ############

rm(list=ls(all=TRUE)) #give R a blank slate
library(tidyverse)
library(lme4)
library(nlme)
library(car)
library(emmeans)
library(piecewiseSEM)
library(visreg)
library(ggthemes)

####### META DATA ######################################
# VWC UNITS are m^3/m^3 multiply by 100 to get as %
# T UNITS are degrees C
# TIMESTAMP is "%Y-%m-%d %H:%M:%S",tz = "UTC"

###### IMPORT AND REFORMAT SOIL SENSOR DATA ##################################################
# set a working directory on your computer
# optional example
setwd("C:/Users/jrudgers/Desktop/SEV Data/Mean-Var Experiment/blue grama mean variance/soil moist temp MVE blue/")
# this step first requires the script MVE_PlainsGrassland_DataProcessing.R
# read in data file
MVE.Plains<-read.csv("MVE_PlainsGrassland_SoilMoistureTemperature.csv", stringsAsFactors = T)
head(MVE.Plains)
# fix data type issues from import
MVE.Plains$TIMESTAMP<-as.POSIXct(MVE.Plains$TIMESTAMP,"%Y-%m-%d",tz = "UTC")
MVE.Plains$depth.f<-as.factor(MVE.Plains$depth)
MVE.Plains$block<-as.factor(MVE.Plains$block)
MVE.Plains$soil.depth<-recode_factor(MVE.Plains$depth.f,'12'="12 cm",'22'="22 cm",'37'="37 cm")
# these steps can be slow to process

########### GET PLOT SCALE VWC VALUES PRE-TREATMENT ###############
#get pretreatment soil moisture data
MVE.VWC.sub<-subset(MVE.Plains,sensor=="VWC"&!is.na(value)&value<1)
MVE.VWC.sub$date<-as.POSIXct(paste(MVE.VWC.sub$year.f,MVE.VWC.sub$month.f,MVE.VWC.sub$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
MVE.pretrt<-subset(MVE.VWC.sub,date<"2019-04-30")
MVE.pretrt.mean<- data.frame(MVE.pretrt %>% group_by(sensor.id) %>% summarize(pretrtmean=mean(value*100),
                                                          pretrtsd=sd(value*100)))
max(MVE.pretrt.mean$pretrtmean)
min(MVE.pretrt.mean$pretrtmean)
summary(MVE.pretrt.mean)

########### ANALYZE TREATMENT EFFECTS ON DAILY VWC - USES THE SIX TRT LEVELS FOR EACH YEAR ########
MVE.VWC.6<-subset(MVE.Plains,sensor=="VWC"&!is.na(value)&value<1)
MVE.VWC.6$date<-as.POSIXct(paste(MVE.VWC.sub$year.f,MVE.VWC.sub$month.f,MVE.VWC.sub$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
MVE.VWC.6.daily<-data.frame(MVE.VWC.6 %>% group_by(year.f,month.f,day.f,date,block,plot,sensor.id,mean.f,var.f,soil.depth,trt_2019,trt_2020,trt_2021) %>% 
                              summarize(mean=mean(value*100,na.rm = TRUE),max=max(value*100,na.rm = TRUE),min=min(value*100,na.rm = TRUE),n=n()))
#merge with pretreatment VWC for covariate analysis
MVE.VWC.6.daily.cov<-merge(MVE.VWC.6.daily,MVE.pretrt.mean,by="sensor.id")
summary(MVE.VWC.6.daily.cov)

### 2020 WATER YEAR ############################
# uses treatment assigned in 2019
MVE.2020<-subset(MVE.VWC.6.daily.cov,date>"2019-10-01"&date<"2020-10-29") #flipped on 10-30
m.2020<-lme(mean ~ trt_2019*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2020)
Anova(m.2020,type=2)
# slopes for figure 5 in Ecosphere MS
emtrends(m.2020,var="trt_2019",~soil.depth)
# note: interaction with soil depth was n.s.
# check model assumptions
hist(resid(m.2020))
qqnorm(resid(m.2020))
plot(m.2020)
#for graphics, use trt_2019 as a factor
MVE.2020$trt_2019.f<-as.factor(MVE.2020$trt_2019)
m.2020.f<-lme(mean ~ trt_2019.f*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2020)
graph.2020<-data.frame(emmeans(m.2020.f,~trt_2019.f|soil.depth))
graph.2020$trt_2019<-as.numeric(as.character(graph.2020$trt_2019.f))
graph.MVE.Plains.2020<-ggplot(data=graph.2020,aes(x=trt_2019,y=emmean,group=soil.depth))+
  facet_grid(rows=vars(soil.depth))+
  theme_tufte(base_size=12,base_family = "sans")+
  geom_errorbar(aes(ymin=(lower.CL),ymax=(upper.CL),width=0.1),position=position_dodge(width=0.4))+
  geom_smooth(formula=y~x,method="lm",color="black",size=0.1)+
  geom_point(size = 2)+
  ylim(0,15)+
  ylab("Volumetric water content (%)")+
  xlab(element_blank())+
  theme(legend.title=element_blank(),legend.position = "top")
graph.MVE.Plains.2020
ggsave("MVE.Plains.2020.percentage.jpg",graph.MVE.Plains.2020,dpi=400,width=4,height=6)

### 2021 WATER YEAR ############################
MVE.2021<-subset(MVE.VWC.6.daily.cov,date>"2020-10-31"&date<"2021-11-08"&sensor.id!="VWC_P2_12")
# note: one sensor failed in this time period, excluded
m.2021<-lme(mean ~ trt_2020*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2021)
Anova(m.2021,type=2)
# note: interaction with soil depth was n.s.
# check model assumptions
hist(resid(m.2021))
qqnorm(resid(m.2021))
plot(m.2021)
#slope values for MS figure
emtrends(m.2021,var="trt_2020",~soil.depth)
#make treatment a factor for graphics
MVE.2021$trt_2020.f<-as.factor(MVE.2021$trt_2020)
m.2021.f<-lme(mean ~ trt_2020.f*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2021)
Anova(m.2021.f,type=2)
graph.2021<-data.frame(emmeans(m.2021.f,~trt_2020.f|soil.depth))
graph.2021$trt_2020<-as.numeric(as.character(graph.2021$trt_2020))
graph.MVE.Plains.2021<-ggplot(data=graph.2021,aes(x=trt_2020,y=emmean,group=soil.depth))+
  facet_grid(rows=vars(soil.depth))+
  theme_tufte(base_size=12,base_family = "sans")+
  geom_errorbar(aes(ymin=(lower.CL),ymax=(upper.CL),width=0.1),position=position_dodge(width=0.4))+
  geom_smooth(formula=y~x,method="lm",color="black",size=0.1)+
  geom_point(size = 2)+
  ylim(0, 15)+
  ylab("Volumetric water content (%)")+
  xlab(element_blank())+
  theme(legend.title=element_blank(),legend.position = "top")
graph.MVE.Plains.2021
ggsave("MVE.Plains.2021.percentage.jpg",graph.MVE.Plains.2021,dpi=400,width=4,height=6)

### 2022 WATER YEAR ############################
MVE.2022<-subset(MVE.VWC.6.daily.cov,date>"2021-11-10"&date<"2022-11-06")#flipped on 9 Nov in prior year, 7 Nov in current year
m.2022<-lme(mean ~ trt_2021*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2022)
Anova(m.2022,type=2)
# slopes for figure in MS
emtrends(m.2022,var="trt_2021",~soil.depth)
# make treatment a factor for graphics
MVE.2022$trt_2021.f<-as.factor(MVE.2022$trt_2021)
m.2022.f<-lme(mean ~ trt_2021.f*soil.depth + pretrtmean, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.2022)
Anova(m.2022.f,type=2)
graph.2022<-data.frame(emmeans(m.2022.f,~trt_2021.f|soil.depth))
graph.2022$trt_2021<-as.numeric(as.character(graph.2022$trt_2021))
graph.MVE.Plains.2022<-ggplot(data=graph.2022,aes(x=trt_2021,y=emmean,group=soil.depth))+
  facet_grid(rows=vars(soil.depth))+
  theme_tufte(base_size=12,base_family = "sans")+
  geom_errorbar(aes(ymin=(lower.CL),ymax=(upper.CL),width=0.1),position=position_dodge(width=0.4))+
  geom_smooth(formula=y~x,method="lm",color="black",size=0.1)+
  geom_point(size = 2)+
  ylim(0, 15)+
  ylab("Volumetric water content (%)")+
  xlab(element_blank())+
  theme(legend.title=element_blank(),legend.position = "top")
graph.MVE.Plains.2022
ggsave("MVE.Plains.2022.percentage.jpg",graph.MVE.Plains.2022,dpi=400,width=4,height=6)

####################### CALCULATE and GRAPH CV in VWC for MORE VARIANCE TREATMENT ################
# summarize to get means for each sensor.id X water year
# 2020 water year
MVE.2020.sensor<- data.frame(MVE.2020 %>% group_by(trt_2019.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                                    summarize(mean=mean(mean),pre=mean(pretrtmean)))
MVE.2020.sensor$wateryear<-as.factor("2020")
MVE.2020.sensor$effectsize<-(MVE.2020.sensor$mean-MVE.2020.sensor$pre)/(MVE.2020.sensor$pre)*100
MVE.2020.sensor <-MVE.2020.sensor %>% rename(MVE_trt=trt_2019.f)
summary(MVE.2020.sensor)
# 2021 water year
MVE.2021.sensor<- data.frame(MVE.2021 %>% group_by(trt_2020.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                               summarize(mean=mean(mean),pre=mean(pretrtmean)))
MVE.2021.sensor$wateryear<-as.factor("2021")
MVE.2021.sensor$effectsize<-(MVE.2021.sensor$mean-MVE.2021.sensor$pre)/(MVE.2021.sensor$pre)*100
MVE.2021.sensor <-MVE.2021.sensor %>% rename(MVE_trt=trt_2020.f)
summary(MVE.2021.sensor)
# 2022 water year
MVE.2022.sensor<- data.frame(MVE.2022 %>% group_by(trt_2021.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                               summarize(mean=mean(mean),pre=mean(pretrtmean)))
MVE.2022.sensor$wateryear<-as.factor("2022")
MVE.2022.sensor$effectsize<-(MVE.2022.sensor$mean-MVE.2022.sensor$pre)/(MVE.2022.sensor$pre)*100
MVE.2022.sensor <-MVE.2022.sensor %>% rename(MVE_trt=trt_2021.f)
summary(MVE.2022.sensor)
# bind together water years
MVE.sensor<-rbind(MVE.2020.sensor,MVE.2021.sensor,MVE.2022.sensor)
#sd function does not work if variable name is mean
MVE.sensor<- MVE.sensor %>% rename(value = mean)
summary(MVE.sensor)
#determine the CV across years 
MVE.Plains.VWC.CV<- data.frame(MVE.sensor %>% group_by(mean.f,var.f,soil.depth) %>% 
                                 summarize(mean=mean(value),sd=sd(value),n=n()))
MVE.Plains.VWC.CV$CV<-MVE.Plains.VWC.CV$sd/MVE.Plains.VWC.CV$mean
mCV<-lm(CV~mean.f*var.f+soil.depth,data=MVE.Plains.VWC.CV)
Anova(mCV,type=2)
mCV.mean<-lm(CV~mean.f*soil.depth,data=MVE.Plains.VWC.CV)
Anova(mCV.mean,type=2)
MVE_CV<-data.frame(emmeans(mCV,~var.f|soil.depth))
MVE_CV$var.f<-recode_factor(MVE_CV$var.f, 'ambient'="Ambient", 'elevated'="More variance")
MVE_CV
#12 cm
(0.200-0.161)/0.161*100
#22 cm
(0.184-0.145)/0.145*100
#37 cm
(0.182-0.143)/0.143*100

# create graph for Figure 6 Ecosphere MS
graph.CV<-ggplot(data=MVE_CV,aes(x=var.f,y=emmean,group=soil.depth,fill=var.f))+
  facet_grid(rows=vars(soil.depth))+
  theme_tufte(base_size=12,base_family = "sans")+
  geom_errorbar(aes(ymin=(emmean-SE),ymax=(emmean+SE),width=0.1),position=position_dodge(width=0.4))+
  geom_col()+
  scale_fill_manual(values = c("gray45", "turquoise3"))+
  ylim(0, 0.24)+
  ylab("CV of volumetric water content")+
  xlab(element_blank())+
  theme(legend.title=element_blank(),legend.position = "none")
graph.CV
ggsave("MVE_CV_graph.jpg",graph.CV,dpi=400,width=3,height=8)

######### GRAPH MICROENVIRONMENT OF PLOT VS RESPONSE TO TRT #################
#read in data 
MVE.microenv.sensor<-read.csv("MVE.microenv.sensor.csv",stringsAsFactors = T)
MVE.microenv.sensor$trt.f<-as.factor(MVE.microenv.sensor$trt.f)
MVE.microenv.sensor$trt.f<-recode_factor(MVE.microenv.sensor$trt.f,
                                         '-75'= "Precipitation -75%",
                                         '-50'= "Precipitation -50%",
                                         '-25'= "Precipitation -25%",
                                         '0'= "Precipitation ambient",
                                         '25'= "Precipitation +25%",
                                         '50'= "Precipitation +50%")
MVE.microenv.sensor$soil.depth<-as.factor(MVE.microenv.sensor$soil.depth)
summary(MVE.microenv.sensor)
#159 obs 11 variables
MVE.microenv.graph<-ggplot(data=MVE.microenv.sensor,aes(x=pretrt.mean,y=effectsize))+
  facet_wrap(vars(trt.f), dir = "v")+
  theme_tufte(base_size=18,base_family = "sans")+
  geom_point(aes(colour = factor(soil.depth)),size=3)+
  geom_smooth(formula=y~x,method="lm",color="black",size=0.5,se=F)+
  ylab("Effect size (%)")+
  xlab("Pre-treatment volumetric water content (%)")+
  scale_colour_manual(values=c("grey70","grey40","black"),labels=c("12 cm","22 cm","37 cm"))+
  theme(legend.title=element_blank(),legend.position = "top")
MVE.microenv.graph
ggsave("MVE.microenv.graph.jpg",MVE.microenv.graph,dpi=400,width=8,height=12)

# microenv analysis
m.micro<-lmer(effectsize~trt.f*pretrt.mean*soil.depth+(1|sensor.id),data=MVE.microenv.sensor,REML=F)
Anova(m.micro,type=2)
emtrends(m.micro,var="pretrt.mean",~soil.depth)
pairs(emtrends(m.micro,var="pretrt.mean",~soil.depth))
# do same analysis for absolute difference between treatment and pretreatment
m.micro.abs<-lmer(abs_diff~trt.f*pretrt.mean*soil.depth+(1|sensor.id),data=MVE.microenv.sensor,REML=F)
Anova(m.micro.abs,type=2)
# same results
# get slope across all depths
m.micro2<-lmer(effectsize~pretrt.mean+(1|sensor.id),data=MVE.microenv.sensor,REML=F)
summary(m.micro2)



########## VWC DAILY SUMMARY for PLOTS of SINGLE SENSORS ####################################
MVE.Plains.VWC<-subset(MVE.Plains,sensor=="VWC"&!is.na(value)&TIMESTAMP>"2019-11-01"&value<1)
summary(MVE.Plains.VWC$TIMESTAMP)
summary(MVE.Plains.VWC$value)
# looks at data post 2018
# could also do on and before "2018-11-04 19:30:00"
# prior to this date, VWC is in different units (%), so would need to split dataset and divide those early dates by 100
# summarize to daily data
MVE.Plains.VWC.daily<- data.frame(MVE.Plains.VWC %>% group_by(year.f,month.f,day.f,block,plot,sensor.id,mean.f,var.f,depth.f) %>% 
                                    summarize(mean=mean(value,na.rm = TRUE),max=max(value,na.rm = TRUE),min=min(value,na.rm = TRUE),n=n()))
# make a timeline variable using julian day
MVE.Plains.VWC.daily$date<-as.POSIXlt(paste(MVE.Plains.VWC.daily$year.f,MVE.Plains.VWC.daily$month.f,MVE.Plains.VWC.daily$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
MVE.Plains.VWC.daily$julian<-as.numeric(format(MVE.Plains.VWC.daily$date, "%j"))
MVE.Plains.VWC.daily$time.f<-as.factor(paste(MVE.Plains.VWC.daily$year,MVE.Plains.VWC.daily$julian,sep="."))
# recode factors to create better names for graphics
MVE.Plains.VWC.daily$soil.depth<-recode_factor(MVE.Plains.VWC.daily$depth.f,'12'="12 cm",'22'="22 cm",'37'="37 cm")
MVE.Plains.VWC.daily$date<-as.POSIXct(paste(MVE.Plains.VWC.daily$year.f,MVE.Plains.VWC.daily$month.f,MVE.Plains.VWC.daily$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
summary(MVE.Plains.VWC.daily)

########## VWC GRAPHS OF SINGLE SENSORS for SUPPLEMENT ####################################
# subset the data to ambient mean, ambient variance
MVE_ambient<-subset(MVE.Plains.VWC.daily, MVE.Plains.VWC.daily$mean.f=="ambient"&MVE.Plains.VWC.daily$var.f=="ambient")
summary(MVE_ambient)
# graph plot 12 as example
MVE_ambient_12<-subset(MVE_ambient, plot=="P12")
graph.VWC.ambient12<-ggplot(data=MVE_ambient_12,aes(x=date,y=mean*100))+
  theme_bw(base_size=20,base_line_size=0.2)+
  #geom_line(size=0.5,position=position_dodge(width=0.4))+
  geom_point(aes(colour = factor(soil.depth)),position=position_dodge(width=0.4))+
  ylab("Volumetric water content (%)")+
  ylim(0,25)+
  xlab(element_blank())+
  scale_colour_manual(values=c("tan4","lightblue3","blue3"),labels=c("12 cm","22 cm","37 cm"))+
  theme(legend.title=element_blank(),legend.position = "bottom")
graph.VWC.ambient12
# subset the data to one treatment: drier elevated variance
MVE_drier.elevated<-subset(MVE.Plains.VWC.daily, MVE.Plains.VWC.daily$mean.f=="drier"&MVE.Plains.VWC.daily$var.f=="elevated")
# graph plot 7
MVE_P7<-subset(MVE_drier.elevated,plot=="P7")
graph.VWC.P7<-ggplot(data=MVE_P7,aes(x=date,y=mean*100))+
  theme_bw(base_size=20,base_line_size=0.2)+
  #geom_line(size=0.5,position=position_dodge(width=0.4))+
  geom_point(aes(colour = factor(soil.depth)),position=position_dodge(width=0.4))+
  ylab("Volumetric water content (%)")+
  ylim(0,25)+
  xlab(element_blank())+
  scale_colour_manual(values=c("tan4","lightblue3","blue3"),labels=c("12 cm","22 cm","37 cm"))+
  theme(legend.title=element_blank(),legend.position = "bottom")
graph.VWC.P7
# graph plot 8
MVE_P8<-subset(MVE_drier.elevated,plot=="P8")
graph.VWC.P8<-ggplot(data=MVE_P8,aes(x=date,y=mean*100))+
  theme_bw(base_size=20,base_line_size=0.2)+
  #geom_line(size=0.5,position=position_dodge(width=0.4))+
  geom_point(aes(colour = factor(soil.depth)),position=position_dodge(width=0.4))+
  ylab("Volumetric water content (%)")+
  ylim(0,25)+
  xlab(element_blank())+
  scale_colour_manual(values=c("tan4","lightblue3","blue3"),labels=c("12 cm","22 cm","37 cm"))+
  theme(legend.title=element_blank(),legend.position = "bottom")
graph.VWC.P8

############### ANALYSIS OF SOIL TEMPERATURE DATA ##########################################
MVE.Temp<-subset(MVE.Plains,sensor=="T"&!is.na(value)&TIMESTAMP>"2018-11-04 19:30:00")
summary(MVE.Temp$TIMESTAMP)
# looks at data post 2018
# on and before "2018-11-04 19:30:00"

########## Temperature DAILY MEANS ####################################
# summarize to daily data
MVE.Plains.Temp.daily<- data.frame(MVE.Temp %>% group_by(year.f,month.f,day.f,block,plot,mean.f,var.f,trt_2020,trt_2020,trt_2021,depth.f,soil.depth,sensor.id) %>% 
                             summarize(mean=mean(value,na.rm = TRUE),max=max(value,na.rm = TRUE),min=min(value,na.rm = TRUE),n=n()))
summary(MVE.Plains.Temp.daily)
# make a timeline variable using julian day
MVE.Plains.Temp.daily$date<-as.POSIXlt(paste(MVE.Plains.Temp.daily$year.f,MVE.Plains.Temp.daily$month.f,MVE.Plains.Temp.daily$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
MVE.Plains.Temp.daily$julian<-as.numeric(format(MVE.Plains.Temp.daily$date, "%j"))
MVE.Plains.Temp.daily$time.f<-as.factor(paste(MVE.Plains.Temp.daily$year,MVE.Plains.Temp.daily$julian,sep="."))
summary(MVE.Plains.Temp.daily)

# recode factors to create better names for graphics
MVE.Plains.Temp.daily$soil.depth<-recode_factor(MVE.Plains.Temp.daily$soil.depth,'12'="12 cm",'22'="22 cm",'37'="37 cm")
MVE.Plains.Temp.daily$date<-as.POSIXct(paste(MVE.Plains.Temp.daily$year.f,MVE.Plains.Temp.daily$month.f,MVE.Plains.Temp.daily$day.f,sep="-"),"%Y-%m-%d",tz = "UTC")
summary(MVE.Plains.Temp.daily$date)
summary(MVE.Plains.Temp.daily$mean)

#export daily data summary from the larger hourly dataset
#write.csv(MVE.daily.Temp,"MVE_PlainsGrassland_Temperature_daily.csv")

########### TEMPERATURE USE SIX TRT LEVELS ###########################
### 2020
MVE.T.2020<-subset(MVE.Plains.Temp.daily, date>"2019-10-01"&date<"2020-10-29")
MVE.T.2020$trt_2019.f<-as.factor(MVE.T.2020$trt_2019)
summary(MVE.T.2020)
m.T.2020<-lme(mean ~ trt_2019.f+soil.depth, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.T.2020)
Anova(m.T.2020,type=2)
emmeans(m.T.2020,~trt_2019.f|soil.depth)
emmeans(m.T.2020,~soil.depth)
### 2021
MVE.T.2021<-subset(MVE.Plains.Temp.daily,date>"2020-10-31"&date<"2021-11-08"&sensor.id!="T_P2_12")
MVE.T.2021$trt_2020.f<-as.factor(MVE.T.2021$trt_2020)
summary(MVE.T.2021$sensor.id)
m.T.2021<-lme(mean ~ trt_2020.f+soil.depth, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.T.2021)
Anova(m.T.2021,type=2)
emmeans(m.T.2021,~trt_2020.f|soil.depth)
### 2022
MVE.T.2022<-subset(MVE.Plains.Temp.daily,date>"2021-11-10"&date<"2022-11-06")
MVE.T.2022$trt_2021.f<-as.factor(MVE.T.2022$trt_2021)
m.T.2022<-lme(mean ~ trt_2021.f+soil.depth, random=~1|block/plot/sensor.id, correlation = corAR1(form = ~ 1 | block/plot/sensor.id), method="ML",data=MVE.T.2022)
Anova(m.T.2022,type=2)
emtrends(m.T.2022,var="trt_2021",~soil.depth)
emmeans(m.T.2020,~soil.depth)

####################### CALCULATE and GRAPH CV in Soil Temp for MORE VARIANCE TREATMENT ################
# summarize to get means for each sensor.id X water year
# 2020 water year
MVE.2020.T.sensor<- data.frame(MVE.T.2020 %>% group_by(trt_2019.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                               summarize(mean=mean(mean)))
MVE.2020.T.sensor$wateryear<-as.factor("2020")
MVE.2020.T.sensor <-MVE.2020.T.sensor %>% rename(MVE_trt=trt_2019.f)
summary(MVE.2020.T.sensor)
# 2021 water year
MVE.2021.T.sensor<- data.frame(MVE.T.2021 %>% group_by(trt_2020.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                               summarize(mean=mean(mean)))
MVE.2021.T.sensor$wateryear<-as.factor("2021")
MVE.2021.T.sensor <-MVE.2021.T.sensor %>% rename(MVE_trt=trt_2020.f)
summary(MVE.2021.T.sensor)
# 2022 water year
MVE.2022.T.sensor<- data.frame(MVE.T.2022 %>% group_by(trt_2021.f,mean.f,var.f,soil.depth,sensor.id) %>% 
                               summarize(mean=mean(mean)))
MVE.2022.T.sensor$wateryear<-as.factor("2022")
MVE.2022.T.sensor <-MVE.2022.T.sensor %>% rename(MVE_trt=trt_2021.f)
summary(MVE.2022.T.sensor)
# bind together water years
MVE.T.sensor<-rbind(MVE.2020.T.sensor,MVE.2021.T.sensor,MVE.2022.T.sensor)
#sd function does not work if variable name is mean
MVE.T.sensor<- MVE.T.sensor %>% rename(value = mean)
summary(MVE.T.sensor)
#determine the CV in temp across years 
MVE.Plains.T.CV<- data.frame(MVE.T.sensor %>% group_by(mean.f,var.f,soil.depth) %>% 
                                 summarize(mean=mean(value),sd=sd(value),n=n()))
MVE.Plains.T.CV$CV<-MVE.Plains.T.CV$sd/MVE.Plains.T.CV$mean
mCV.T<-lm(CV~mean.f*var.f+soil.depth,data=MVE.Plains.T.CV)
Anova(mCV.T,type=2)
mCV.T2<-lm(CV~var.f*soil.depth,data=MVE.Plains.T.CV)
Anova(mCV.T2,type=2)
MVE_CV.T<-data.frame(emmeans(mCV.T,~var.f))
MVE_CV.T$var.f<-recode_factor(MVE_CV.T$var.f, 'ambient'="Ambient", 'elevated'="More variance")
MVE_CV.T
#all depths
(0.06-0.033)/0.033*100





