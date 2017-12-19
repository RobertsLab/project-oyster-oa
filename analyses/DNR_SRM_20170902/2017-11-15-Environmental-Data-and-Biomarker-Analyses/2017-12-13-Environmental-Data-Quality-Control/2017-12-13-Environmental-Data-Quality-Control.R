#In this script, I'll quality control the environmental data. Using tidal data, I can replace any readings from when the probe was out of water with N/As.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT AND FORMAT DATA ####

tideData <- read.csv("../../data/DNR/2017-12-13-Tidal-Data-by-Site.csv", header = TRUE) #Import the tide data
head(tideData) #Confirm import
tideData$Date <- as.Date(tideData$Date, format = "%m/%d/%y") #Convert entries to dates
tideData$DateTime <- paste(tideData$Date, "", tideData$Time) #Create new DateTime column to easily merge tide and environmental data
colnames(tideData) <- c("Date", "Time", "CI-Tide", "FB-Tide", "PG-Tide", "SK-Tide", "WB-Tide", "DateTime")
head(tideData) #Confirm changes

pHDOData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA") #Import file with pH and DO data
head(pHDOData) #Confirm import
colnames(pHDOData) #View column names

pHData <- pHDOData[,c(2:3, 5, 7, 8, 10, 12)] #Subset only the bare outplant pH data. Do not include DateTime column, since the formatting is off
head(pHData) #Confirm subset
colnames(pHData) <- c("Date", "Time", "WBB-pH", "SKB-pH", "PGB-pH", "CIB-pH", "FBB-pH") #Rename columns
pHData$Date <- as.Date(pHData$Date, format = "%m/%d/%y") #Convert entries to dates
pHData$DateTime <- paste(pHData$Date, pHData$Time) #Create new DateTime column to easily merge tide and environmental data
head(pHData) #Confirm changes

DOData <- pHDOData[,c(2:3, 23, 25, 27, 29, 31)] #Subset only the bare outplant DO data
head(DOData) #Confirm subset
colnames(DOData) <- c("Date", "Time", "WBB-DO", "SKB-DO", "PGB-DO", "CIB-DO", "FBB-DO") #Rename columns
DOData$Date <- as.Date(DOData$Date, format = "%m/%d/%y") #Convert entries to dates
DOData$DateTime <- paste(DOData$Date, DOData$Time) #Create new DateTime column to easily merge tide and environmental data
head(DOData) #Confirm changes

salinityData <- read.csv("../../data/DNR/2017-11-25-Calculated-Salinity-Output-from-Micah.csv", header = TRUE, na.strings = "NA") #Import salinity data
head(salinityData) #Confirm import
colnames(salinityData) #Get column names
salinityData <- salinityData[,c(1:2, 3, 7, 11, 15, 18)] #Subset only the salinity information from bare outplants. Needed to use PGE instead of PGB since PGB has no salinity data.
head(salinityData) #Confirm subset
colnames(salinityData) <- c("Date", "Time", "CIB-Salinity", "FBB-Salinity", "PGE-Salinity", "SKB-Salinity", "WBB-Salinity") #Rename columns
salinityData$Date <- as.Date(salinityData$Date, format = "%d/%m/%Y") #Convert entries to dates
salinityData$DateTime <- paste(salinityData$Date, salinityData$Time) #Create new DateTime column to easily merge tide and environmental data
head(salinityData)

#### MERGE TIDAL DATA WITH ENVIRONMENTAL DATA ####

pHTideData <- merge(x = pHData, y = tideData, by = "DateTime") #Merge pH and tide data
head(pHTideData) #Confirm merge
colnames(pHTideData) #Get column names
pHTideData <- pHTideData[, -c(9:10)] #Remove redundant date and time columns
colnames(pHTideData) <- c("DateTime", "Date", "Time", "WBB-pH", "SKB-pH", "PGB-pH", "CIB-pH", "FBB-pH", "CI-Tide", "FB-Tide", "PG-Tide", "SK-Tide", "WB-Tide") #Change column names
head(pHTideData) #Confirm changes

DOTideData <- merge(x = DOData, y = tideData, by = "DateTime") #Merge DO and tide data
head(DOTideData) #Confirm merge
colnames(DOTideData) #Get column names
DOTideData <- DOTideData[, -c(9:10)] #Remove redundant date and time columns
colnames(DOTideData) <- c("DateTime", "Date", "Time", "WBB-DO", "SKB-DO", "PGB-DO", "CIB-DO", "FBB-DO", "CI-Tide", "FB-Tide", "PG-Tide", "SK-Tide", "WB-Tide") #Change column names
head(DOTideData) #Confirm changes

salinityTideData <- merge(x = salinityData, y = tideData, by = "DateTime")
