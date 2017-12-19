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
temperatureData <- data.frame(temperatureData$DateTime,
                              temperatureData$Date,
                              temperatureData$Time,
                              temperatureData$CIB,
                              temperatureData$FBB,
                              temperatureData$PGB,
                              temperatureData$SKB,
                              temperatureData$WBB) #Reorganize columns
colnames(temperatureData) <- c("DateTime", "Date", "Time", "CIB", "FBB", "PGB", "SKB", "WBB") #Rename columns
head(temperatureData) #Confirm changes

#pH
pHData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-pH-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA") #Import pH data
head(pHData) #Confirm import
colnames(pHData) #Get column names
pHData <- pHData[,-c(1, 10:14)] #Save only pH data as a new dataframe
colnames(pHData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
pHData <- data.frame(pHData$DateTime,
                     pHData$Date,
                     pHData$Time,
                     pHData$CIB,
                     pHData$FBB,
                     pHData$PGB,
                     pHData$SKB,
                     pHData$WBB) #Reorganize columns
colnames(pHData) <- c("DateTime", "Date", "Time", "CIB", "FBB", "PGB", "SKB", "WBB") #Rename columns
head(pHData) #Confirm changes

#DO
DOData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-DO-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA") #Import DO data
head(DOData) #Confirm import
colnames(DOData) #Get column names
DOData <- DOData[,-c(1, 10:14)] #Save DO data as a new dataframe
colnames(DOData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
DOData <- data.frame(DOData$DateTime,
                     DOData$Date,
                     DOData$Time,
                     DOData$CIB,
                     DOData$FBB,
                     DOData$PGB,
                     DOData$SKB,
                     DOData$WBB) #Reorganize columns
colnames(DOData) <- c("DateTime", "Date", "Time", "CIB", "FBB", "PGB", "SKB", "WBB") #Rename columns
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
temperatureVariablesofInterest <- data.frame("blank1" = rep(0, times = 12),
                                             "blank2" = rep(0, times = 12),
                                             "blank3" = rep(0, times = 12),
                                             "CI" = rep(0, times = 12),
                                             "FB" = rep(0, times = 12),
                                             "PG" = rep(0, times = 12),
                                             "SK" = rep(0, times = 12),
                                             "WB" = rep(0, times = 12)) #Create an empty dataframe. Columns 4-8 are for sites, following the pattern of all environmental variable datasets.
rownames(temperatureVariablesofInterest) <- c("Maximum", "Minimum", "Range", "Mean", "Variance", "StandardDeviation", "Percent±2SD", "FirstQuartile", "Median", "ThirdQuartile", "IQR", "Percent±1.5IQR")
head(temperatureVariablesofInterest) #Confirm dataframe creation

nSites <- 8 #Sites are in columns 4-8
for(i in 4:nSites) { #For each site
  temperatureVariablesofInterest[1, i] <- max(temperatureData[, i], na.rm = TRUE) #Calculate maximum
  temperatureVariablesofInterest[2, i] <- min(temperatureData[, i], na.rm = TRUE) #Calculate minimum
  temperatureVariablesofInterest[3, i] <- ((max(temperatureData[, i], na.rm = TRUE)) - (min(temperatureData[, i], na.rm = TRUE))) #Calculate range
  temperatureVariablesofInterest[4, i] <- mean(temperatureData[, i], na.rm = TRUE) #Calculate mean
  temperatureVariablesofInterest[5, i] <- var(temperatureData[, i], na.rm = TRUE) #Calculate variance
  temperatureVariablesofInterest[6, i] <- sqrt(var(temperatureData[, i], na.rm = TRUE)) #Calculate SD
  temperatureVariablesofInterest[7, i] <- ((sum((temperatureData[, i] > (mean(temperatureData[, i], na.rm = TRUE) + (2*sqrt(var(temperatureData[, i], na.rm = TRUE))))), (temperatureData[, i] < (mean(temperatureData[, i], na.rm = TRUE) - (2*sqrt(var(temperatureData[, i], na.rm = TRUE))))), na.rm = TRUE))/(length(temperatureData[, i])))*100 #Calculate percentage of data more than 2 SDs away from mean
  temperatureVariablesofInterest[8, i] <- as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[2]) #Calculate first quartile
  temperatureVariablesofInterest[9, i] <- as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[3]) #Calculate median
  temperatureVariablesofInterest[10, i] <- as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[4]) #Calculate third quartile
  temperatureVariablesofInterest[11, i] <- as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[2]) #Calculate IQR
  temperatureVariablesofInterest[12, i] <- ((sum((temperatureData[, i] > as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[4]) + 1.5*(as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[2]))), (temperatureData[, i] < as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[2]) - 1.5*(as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(temperatureData[, i], na.rm = TRUE)[2]))), na.rm = TRUE))/(length(temperatureData[, i])))*100 #Calculate percentage of data more than 1.5*IQR away from the upper or lower quartiles
} #Calculate variables of interest at each site for temperature data
head(temperatureVariablesofInterest) #Confirm data table is properly filled
temperatureVariablesofInterest <- temperatureVariablesofInterest[, -c(1:3)] #Remove blank columns
head(temperatureVariablesofInterest) #Confirm changes

#pH
pHVariablesofInterest <- data.frame("blank1" = rep(0, times = 12),
                                    "blank2" = rep(0, times = 12),
                                    "blank3" = rep(0, times = 12),
                                    "CI" = rep(0, times = 12),
                                    "FB" = rep(0, times = 12),
                                    "PG" = rep(0, times = 12),
                                    "SK" = rep(0, times = 12),
                                    "WB" = rep(0, times = 12)) #Create an empty dataframe. Columns 4-8 are for sites, following the pattern of all environmental variable datasets.
rownames(pHVariablesofInterest) <- c("Maximum", "Minimum", "Range", "Mean", "Variance", "StandardDeviation", "Percent±2SD", "FirstQuartile", "Median", "ThirdQuartile", "IQR", "Percent±1.5IQR")
head(pHVariablesofInterest) #Confirm dataframe creation

for(i in 4:nSites) { #For each site
  pHVariablesofInterest[1, i] <- max(pHData[, i], na.rm = TRUE) #Calculate maximum
  pHVariablesofInterest[2, i] <- min(pHData[, i], na.rm = TRUE) #Calculate minimum
  pHVariablesofInterest[3, i] <- ((max(pHData[, i], na.rm = TRUE)) - (min(pHData[, i], na.rm = TRUE))) #Calculate range
  pHVariablesofInterest[4, i] <- mean(pHData[, i], na.rm = TRUE) #Calculate mean
  pHVariablesofInterest[5, i] <- var(pHData[, i], na.rm = TRUE) #Calculate variance
  pHVariablesofInterest[6, i] <- sqrt(var(pHData[, i], na.rm = TRUE)) #Calculate SD
  pHVariablesofInterest[7, i] <- ((sum((pHData[, i] > (mean(pHData[, i], na.rm = TRUE) + (2*sqrt(var(pHData[, i], na.rm = TRUE))))), (pHData[, i] < (mean(pHData[, i], na.rm = TRUE) - (2*sqrt(var(pHData[, i], na.rm = TRUE))))), na.rm = TRUE))/(length(pHData[, i])))*100 #Calculate percentage of data more than 2 SDs away from mean
  pHVariablesofInterest[8, i] <- as.numeric(quantile(pHData[, i], na.rm = TRUE)[2]) #Calculate first quartile
  pHVariablesofInterest[9, i] <- as.numeric(quantile(pHData[, i], na.rm = TRUE)[3]) #Calculate median
  pHVariablesofInterest[10, i] <- as.numeric(quantile(pHData[, i], na.rm = TRUE)[4]) #Calculate third quartile
  pHVariablesofInterest[11, i] <- as.numeric(quantile(pHData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(pHData[, i], na.rm = TRUE)[2]) #Calculate IQR
  pHVariablesofInterest[12, i] <- ((sum((pHData[, i] > as.numeric(quantile(pHData[, i], na.rm = TRUE)[4]) + 1.5*(as.numeric(quantile(pHData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(pHData[, i], na.rm = TRUE)[2]))), (pHData[, i] < as.numeric(quantile(pHData[, i], na.rm = TRUE)[2]) - 1.5*(as.numeric(quantile(pHData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(pHData[, i], na.rm = TRUE)[2]))), na.rm = TRUE))/(length(pHData[, i])))*100 #Calculate percentage of data more than 1.5*IQR away from the upper or lower quartiles
} #Calculate variables of interest at each site for pH data
head(pHVariablesofInterest) #Confirm data table is properly filled
pHVariablesofInterest <- pHVariablesofInterest[, -c(1:3)] #Remove blank columns
head(pHVariablesofInterest) #Confirm changes

#DO
DOVariablesofInterest <- data.frame("blank1" = rep(0, times = 12),
                                    "blank2" = rep(0, times = 12),
                                    "blank3" = rep(0, times = 12),
                                    "CI" = rep(0, times = 12),
                                    "FB" = rep(0, times = 12),
                                    "PG" = rep(0, times = 12),
                                    "SK" = rep(0, times = 12),
                                    "WB" = rep(0, times = 12)) #Create an empty dataframe. Columns 4-8 are for sites, following the pattern of all environmental variable datasets.
rownames(DOVariablesofInterest) <- c("Maximum", "Minimum", "Range", "Mean", "Variance", "StandardDeviation", "Percent±2SD", "FirstQuartile", "Median", "ThirdQuartile", "IQR", "Percent±1.5IQR")
head(DOVariablesofInterest) #Confirm dataframe creation

for(i in 4:nSites) { #For each site
  DOVariablesofInterest[1, i] <- max(DOData[, i], na.rm = TRUE) #Calculate maximum
  DOVariablesofInterest[2, i] <- min(DOData[, i], na.rm = TRUE) #Calculate minimum
  DOVariablesofInterest[3, i] <- ((max(DOData[, i], na.rm = TRUE)) - (min(DOData[, i], na.rm = TRUE))) #Calculate range
  DOVariablesofInterest[4, i] <- mean(DOData[, i], na.rm = TRUE) #Calculate mean
  DOVariablesofInterest[5, i] <- var(DOData[, i], na.rm = TRUE) #Calculate variance
  DOVariablesofInterest[6, i] <- sqrt(var(DOData[, i], na.rm = TRUE)) #Calculate SD
  DOVariablesofInterest[7, i] <- ((sum((DOData[, i] > (mean(DOData[, i], na.rm = TRUE) + (2*sqrt(var(DOData[, i], na.rm = TRUE))))), (DOData[, i] < (mean(DOData[, i], na.rm = TRUE) - (2*sqrt(var(DOData[, i], na.rm = TRUE))))), na.rm = TRUE))/(length(DOData[, i])))*100 #Calculate percentage of data more than 2 SDs away from mean
  DOVariablesofInterest[8, i] <- as.numeric(quantile(DOData[, i], na.rm = TRUE)[2]) #Calculate first quartile
  DOVariablesofInterest[9, i] <- as.numeric(quantile(DOData[, i], na.rm = TRUE)[3]) #Calculate median
  DOVariablesofInterest[10, i] <- as.numeric(quantile(DOData[, i], na.rm = TRUE)[4]) #Calculate third quartile
  DOVariablesofInterest[11, i] <- as.numeric(quantile(DOData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(DOData[, i], na.rm = TRUE)[2]) #Calculate IQR
  DOVariablesofInterest[12, i] <- ((sum((DOData[, i] > as.numeric(quantile(DOData[, i], na.rm = TRUE)[4]) + 1.5*(as.numeric(quantile(DOData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(DOData[, i], na.rm = TRUE)[2]))), (DOData[, i] < as.numeric(quantile(DOData[, i], na.rm = TRUE)[2]) - 1.5*(as.numeric(quantile(DOData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(DOData[, i], na.rm = TRUE)[2]))), na.rm = TRUE))/(length(DOData[, i])))*100 #Calculate percentage of data more than 1.5*IQR away from the upper or lower quartiles
} #Calculate variables of interest at each site for DO data
head(DOVariablesofInterest) #Confirm data table is properly filled
DOVariablesofInterest <- DOVariablesofInterest[, -c(1:3)] #Remove blank columns
head(DOVariablesofInterest) #Confirm changes

#Salinity
salinityVariablesofInterest <- data.frame("blank1" = rep(0, times = 12),
                                          "blank2" = rep(0, times = 12),
                                          "blank3" = rep(0, times = 12),
                                          "CI" = rep(0, times = 12),
                                          "FB" = rep(0, times = 12),
                                          "PG" = rep(0, times = 12),
                                          "SK" = rep(0, times = 12),
                                          "WB" = rep(0, times = 12)) #Create an empty dataframe. Columns 4-8 are for sites, following the pattern of all environmental variable datasets.
rownames(salinityVariablesofInterest) <- c("Maximum", "Minimum", "Range", "Mean", "Variance", "StandardDeviation", "Percent±2SD", "FirstQuartile", "Median", "ThirdQuartile", "IQR", "Percent±1.5IQR")
head(salinityVariablesofInterest) #Confirm dataframe creation

for(i in 4:nSites) { #For each site
  salinityVariablesofInterest[1, i] <- max(salinityData[, i], na.rm = TRUE) #Calculate maximum
  salinityVariablesofInterest[2, i] <- min(salinityData[, i], na.rm = TRUE) #Calculate minimum
  salinityVariablesofInterest[3, i] <- ((max(salinityData[, i], na.rm = TRUE)) - (min(salinityData[, i], na.rm = TRUE))) #Calculate range
  salinityVariablesofInterest[4, i] <- mean(salinityData[, i], na.rm = TRUE) #Calculate mean
  salinityVariablesofInterest[5, i] <- var(salinityData[, i], na.rm = TRUE) #Calculate variance
  salinityVariablesofInterest[6, i] <- sqrt(var(salinityData[, i], na.rm = TRUE)) #Calculate SD
  salinityVariablesofInterest[7, i] <- ((sum((salinityData[, i] > (mean(salinityData[, i], na.rm = TRUE) + (2*sqrt(var(salinityData[, i], na.rm = TRUE))))), (salinityData[, i] < (mean(salinityData[, i], na.rm = TRUE) - (2*sqrt(var(salinityData[, i], na.rm = TRUE))))), na.rm = TRUE))/(length(salinityData[, i])))*100 #Calculate percentage of data more than 2 SDs away from mean
  salinityVariablesofInterest[8, i] <- as.numeric(quantile(salinityData[, i], na.rm = TRUE)[2]) #Calculate first quartile
  salinityVariablesofInterest[8, i] <- as.numeric(quantile(salinityData[, i], na.rm = TRUE)[3]) #Calculate median
  salinityVariablesofInterest[10, i] <- as.numeric(quantile(salinityData[, i], na.rm = TRUE)[4]) #Calculate third quartile
  salinityVariablesofInterest[11, i] <- as.numeric(quantile(salinityData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(salinityData[, i], na.rm = TRUE)[2]) #Calculate IQR
  salinityVariablesofInterest[12, i] <- ((sum((salinityData[, i] > as.numeric(quantile(salinityData[, i], na.rm = TRUE)[4]) + 1.5*(as.numeric(quantile(salinityData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(salinityData[, i], na.rm = TRUE)[2]))), (salinityData[, i] < as.numeric(quantile(salinityData[, i], na.rm = TRUE)[2]) - 1.5*(as.numeric(quantile(salinityData[, i], na.rm = TRUE)[4]) - as.numeric(quantile(salinityData[, i], na.rm = TRUE)[2]))), na.rm = TRUE))/(length(salinityData[, i])))*100 #Calculate percentage of data more than 1.5*IQR away from the upper or lower quartiles
} #Calculate variables of interest at each site for salinity data
head(salinityVariablesofInterest) #Confirm data table is properly filled
salinityVariablesofInterest <- salinityVariablesofInterest[, -c(1:3)] #Remove blank columns
head(salinityVariablesofInterest) #Confirm changes

#### SAVE TABLES ####

#Temperature
#write.csv(temperatureVariablesofInterest, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-19-Temperature-Data-Variables-of-Interest.csv")

#pH
#write.csv(pHVariablesofInterest, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-19-pH-Data-Variables-of-Interest.csv")

#DO
#write.csv(DOVariablesofInterest, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-19-DO-Data-Variables-of-Interest.csv")

#Salinity
#write.csv(salinityVariablesofInterest, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-19-Salinity-Data-Variables-of-Interest.csv")