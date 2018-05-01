#A code to calibrate Tris calibrate your pH sensor on the total scale and calculate pH from mV and Temperature
#pH data is exported as a csv file with the name "Data Today's Date" into your working directory
#Modified by Nyssa Silbiger 6/25/14
#Edited 1/7/2016
#Edited by Yaamini Venkataraman on 4/30/2018

#### SET WORKING DIRECTORY ####
main<-setwd(getwd())

#### LOAD LIBRARIES AND SET FUNCTIONALITY ####
install.packages("seacarb") #Install seacarb package
library(seacarb)
rm(list=ls())

#### LOAD DATA FILES (ONLY CHANGE THESE TWO LINES OF CODE) ####
filename.calib <- "2018-04-30-Calibration-Measurements/2017-04-08-Calibration-Measurements.csv" #Calibration file
filename.data <- "2018-04-30-Grab-Samples/2017-04-08-Grab-Samples.csv" #Data file

#### READ IN CALIBRATION DATA ####
z<-read.table(filename.calib, header = T, sep=",", na.string="NA") 
#mVTris is milivolts of the tris calibration
#TTris is temperature of the tris calibration 
#there needs to be at least three points of Tris values at different temperature

attach(z) #this lets me search the variable names within "z"

#load your own data
Data <- read.table(filename.data, header=T, sep=",", na.string="NA") #make sure to change "your data"
#mV is milivolts from your data
#Tin is temperature from your data (measured at the time of collection)

#### RUN THE CALIBRATION ####
#fit a linear model between temperature and mV
mVTris_t<-lm(mVTris~TTris)

print(summary(mVTris_t))

#plot temperature vs millivolt -- make sure it looks like a straight line
plot(TTris, mVTris, xlab="Temperature", ylab="mV")
abline(mVTris_t)

#### FIT PH FOR DATA ####
attach(Data)
#variables
#Tin is your temperature in situ
#mV is your milivolt

#constants
R<- 8.31451 #gas constant
Far<-96485.309 #Faraday constant

#Tin<-T#Temp measured in situ

#calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
mvTris<-Tin*mVTris_t$coefficients[2]+mVTris_t$coefficients[1]

STris<-34.5 #salinity of the Tris

phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(Tin+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(Tin+273.15)-0.11149858*(Tin+273.15)

# calculate pherror
Tris<-lm(phTris~Tin) #Linear model of temperate by pH of tris
TrisCalc<-25*Tris$coefficients[2]+Tris$coefficients[1] #calculate the pH at 25 deg C
pHError<-((TrisCalc-8.0835)/8.0835)*100 #percent error of your probe measurement


#calculate the pH of your samples
pH<-phTris+(mvTris/1000-mV/1000)/(R*(Data$Tin+273.15)*log(10)/Far)

Data<-cbind(Data,pH)

date<-as.Date(Data$Date[1], format = "%m/%d/%y")
write.table(Data,paste("Data",date[1],".csv"),sep=",", col.names = NA)