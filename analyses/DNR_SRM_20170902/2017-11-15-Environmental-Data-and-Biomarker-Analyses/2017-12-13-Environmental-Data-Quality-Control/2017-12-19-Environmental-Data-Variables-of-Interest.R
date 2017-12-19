#In this script, I'll analyze variables of interest for my temperature, pH, DO and salinity data.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT AND FORMAT DATA ####

#Temperature
temperatureData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA") #Import temperature data
head(temperatureData) #Confirm import
colnames(temperatureData) #Get column names
temperatureData <- temperatureData[,c(2:3, seq(from = 33, to = 41, by = 2))] #Save temperature data from bare outplants as a new dataframe
head(temperatureData) #Confirm subset
colnames(temperatureData) <- c("Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
temperatureData$Date <- as.Date(temperatureData$Date, format = "%m/%d/%y") #Convert entries to dates
temperatureData$DateTime <- paste(temperatureData$Date, temperatureData$Time) #Create new DateTime column to easily merge tide and environmental data
head(temperatureData) #Confirm changes

#pH
pHData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-pH-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA") #Import pH data
head(pHData) #Confirm import
colnames(pHData) #Get column names
pHData <- pHData[,-c(1, 10:14)] #Save only pH data as a new dataframe
colnames(pHData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
head(pHData) #Confirm changes

#DO
DOData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-DO-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA") #Import DO data
head(DOData) #Confirm import
colnames(DOData) #Get column names
DOData <- DOData[,-c(1, 10:14)] #Save DO data as a new dataframe
colnames(DOData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
head(DOData) #Confirm changes

#Salinity
salinityData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Salinity-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA") #Import salinity data
head(salinityData) #Confirm import
colnames(salinityData) #Get column names
salinityData <- salinityData[, -c(1, 10:14)] #Save salinity data as a new dataframe
colnames(salinityData) <- c("DateTime", "Date", "Time", "CIB", "FBB", "PGE", "SKB", "WBB") #Rename columns
head(salinityData) #Confirm changes

#### CALCULATE IMPORTANT VARIABLES ####

#Temperature
max(temperatureData$WBB, na.rm = TRUE) #Calculate maximum
min(temperatureData$WBB, na.rm = TRUE) #Calculate minimum
mean(temperatureData$WBB, na.rm = TRUE) #Calculate mean
var(temperatureData$WBB, na.rm = TRUE) #Calculate variance
sqrt(var(temperatureData$WBB, na.rm = TRUE)) #Calculate SD
((sum((temperatureData$WBB > (mean(temperatureData$WBB, na.rm = TRUE) + (2*sqrt(var(temperatureData$WBB, na.rm = TRUE))))), (temperatureData$WBB < (mean(temperatureData$WBB, na.rm = TRUE) - (2*sqrt(var(temperatureData$WBB, na.rm = TRUE))))), na.rm = TRUE))/(length(temperatureData$WBB)))*100 #Calculate percentage of data more than 2 SDs away from mean
as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[2]) #Calculate first quartile
as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[3]) #Calculate median
as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[4]) #Calculate third quartile
as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[2]) #Calculate IQR
((sum((temperatureData$WBB > as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[4]) + 1.5*(as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[2]))), (temperatureData$WBB < as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[2]) - 1.5*(as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData$WBB, na.rm = TRUE)[2]))), na.rm = TRUE))/(length(temperatureData$WBB)))*100 #Calculate percentage of data more than 1.5*IQR away from the upper or lower quartiles