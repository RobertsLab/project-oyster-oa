#In this script, I'll quality control the environmental data. Using tidal data, I can replace any readings from when the probe was out of water with N/As.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT AND FORMAT DATA ####

tideData <- read.csv("../../data/DNR/2017-12-13-Tidal-Data-by-Site.csv", header = TRUE, strip.white = TRUE) #Import the tide data
head(tideData) #Confirm import
tideData$Date <- as.Date(tideData$Date, format = "%m/%d/%y") #Convert entries to dates
tideData$DateTime <- paste(tideData$Date, tideData$Time) #Create new DateTime column to easily merge tide and environmental data
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

salinityData <- read.csv("../../data/DNR/2018-05-30-Fixed-Salinity-from-Micah.csv", header = TRUE, na.strings = "NA", strip.white = TRUE) #Import salinity data and remove white space from end of Date and Time columns
head(salinityData) #Confirm import
colnames(salinityData) #Get column names
salinityData <- salinityData[,c(1:2, 4, 10, 12, 14, 20)] #Subset only the salinity information from bare outplants. Needed to use PGE instead of PGB since PGB has no salinity data. Also use FBE and WBE due to probe burial at bare sites.
head(salinityData) #Confirm subset
colnames(salinityData) <- c("Date", "Time", "CIB-Salinity", "FBB-Salinity", "PGE-Salinity", "SKB-Salinity", "WBB-Salinity") #Rename columns
salinityData$Date <- as.Date(salinityData$Date, format = "%d/%m/%Y") #Convert entries to dates
salinityData$DateTime <- paste(salinityData$Date, salinityData$Time) #Create new DateTime column to easily merge tide and environmental data.
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

salinityTideData <- merge(x = salinityData, y = tideData, by = "DateTime") #Merge salinity and tide data
head(salinityTideData) #Confirm merge
colnames(salinityTideData) #Get column names
salinityTideData <- salinityTideData[, -c(9:10)] #Remove redundant date and time columns
colnames(salinityTideData) <- c("DateTime", "Date", "Time", "CIB-Salinity", "FBB-Salinity", "PGE-Salinity", "SKB-Salinity", "WBB-Salinity", "CI-Tide", "FB-Tide", "PG-Tide", "SK-Tide", "WB-Tide") #Change column names
head(salinityTideData) #Confirm changes

#### QUANTIFY EXPOSURE TIMES ####

#Count exposure intervals, multiply by 10 to convert to minutes, and divide by 60 to convert to hours.
((length(which(tideData$`CI-Tide` <= 1))*10)/60) #125.6667
((length(which(tideData$`FB-Tide` <= 1))*10)/60) #188
((length(which(tideData$`PG-Tide` <= 1))*10)/60) #146.8333
((length(which(tideData$`SK-Tide` <= 1))*10)/60) #138.3333
((length(which(tideData$`WB-Tide` <= 1))*10)/60) #113.6667

#### REMOVE EXPOSURE TIMES ####

#pH Data
pHTideData$`CIB-pH`[pHTideData$`CI-Tide` <= 1] <- NA #Replace CIB-pH values with "NA" when tide is less than 1
pHTideData$`FBB-pH`[pHTideData$`FB-Tide` <= 1] <- NA #Replace FBB-pH values with "NA" when tide is less than 1
pHTideData$`PGB-pH`[pHTideData$`PG-Tide` <= 1] <- NA #Replace PGB-pH values with "NA" when tide is less than 1
pHTideData$`SKB-pH`[pHTideData$`SK-Tide` <= 1] <- NA #Replace SKB-pH values with "NA" when tide is less than 1
pHTideData$`WBB-pH`[pHTideData$`WB-Tide` <= 1] <- NA #Replace WBB-pH values with "NA" when tide is less than 1

#Convert to numeric values
pHTideData$`CIB-pH` <- as.numeric(pHTideData$`CIB-pH`)
pHTideData$`FBB-pH` <- as.numeric(pHTideData$`FBB-pH`)
pHTideData$`PGB-pH` <- as.numeric(pHTideData$`PGB-pH`)
pHTideData$`SKB-pH` <- as.numeric(pHTideData$`SKB-pH`)
pHTideData$`WBB-pH` <- as.numeric(pHTideData$`WBB-pH`)

#DO Data
DOTideData$`CIB-DO`[DOTideData$`CI-Tide` <= 1] <- NA #Replace CIB-DO values with "NA" when tide is less than 1
DOTideData$`FBB-DO`[DOTideData$`FB-Tide` <= 1] <- NA #Replace FBB-DO values with "NA" when tide is less than 1
DOTideData$`PGB-DO`[DOTideData$`PG-Tide` <= 1] <- NA #Replace PGB-DO values with "NA" when tide is less than 1
DOTideData$`SKB-DO`[DOTideData$`SK-Tide` <= 1] <- NA #Replace SKB-DO values with "NA" when tide is less than 1
DOTideData$`WBB-DO`[DOTideData$`WB-Tide` <= 1] <- NA #Replace WBB-DO values with "NA" when tide is less than 1

#Convert to numeric values
DOTideData$`CIB-DO` <- as.numeric(DOTideData$`CIB-DO`)
DOTideData$`FBB-DO` <- as.numeric(DOTideData$`FBB-DO`)
DOTideData$`PGB-DO` <- as.numeric(DOTideData$`PGB-DO`)
DOTideData$`SKB-DO` <- as.numeric(DOTideData$`SKB-DO`)
DOTideData$`WBB-DO` <- as.numeric(DOTideData$`WBB-DO`)

#Salinity Data
salinityTideData$`CIB-Salinity`[salinityTideData$`CIB-Salinity` <= 1] <- NA #Replace CIB-Salinity values with "NA" when tide is less than 1
salinityTideData$`FBB-Salinity`[salinityTideData$`FB-Salinity` <= 1] <- NA #Replace FBB-Salinity values with "NA" when tide is less than 1
salinityTideData$`PGE-Salinity`[salinityTideData$`PG-Salinity` <= 1] <- NA #Replace PGE-Salinity values with "NA" when tide is less than 1
salinityTideData$`SKB-Salinity`[salinityTideData$`SK-Salinity` <= 1] <- NA #Replace SKB-Salinity values with "NA" when tide is less than 1
salinityTideData$`WBB-Salinity`[salinityTideData$`WB-Salinity` <= 1] <- NA #Replace WBB-Salinity values with "NA" when tide is less than 1

#Convert to numeric values
salinityTideData$`CIB-Salinity` <- as.numeric(salinityTideData$`CIB-Salinity`)
salinityTideData$`FBB-Salinity` <- as.numeric(salinityTideData$`FBB-Salinity`)
salinityTideData$`PGE-Salinity` <- as.numeric(salinityTideData$`PGE-Salinity`)
salinityTideData$`SKB-Salinity` <- as.numeric(salinityTideData$`SKB-Salinity`)
salinityTideData$`WBB-Salinity` <- as.numeric(salinityTideData$`WBB-Salinity`)

#### REMOVE OUTLIERS ####

#pH Data
nSites <- 8 #Sites are from columns 4 to 8
for(i in 4:nSites) { #For individual site data
  upperBound <- as.numeric((quantile(pHTideData[, i], na.rm = TRUE)[4]) + (1.5*(quantile(pHTideData[, i], na.rm = TRUE)[4] - quantile(pHTideData[, i], na.rm = TRUE)[2]))) #Calculate upper bound
  lowerBound <- as.numeric((quantile(pHTideData[, i], na.rm = TRUE)[2]) - (1.5*(quantile(pHTideData[, i], na.rm = TRUE)[4] - quantile(pHTideData[, i], na.rm = TRUE)[2]))) #Calculate lower bound
  pHTideData[, i][pHTideData[, i] > upperBound] <- NA #Replace any values higher than upper bound with NA
  pHTideData[, i][pHTideData[, i] < lowerBound] <- NA #Replace any values lower than upper bound with NA
} #Replace outliers with NA values

#DO Data
for(i in 4:nSites) { #For individual site data
  upperBound <- as.numeric((quantile(DOTideData[, i], na.rm = TRUE)[4]) + (1.5*(quantile(DOTideData[, i], na.rm = TRUE)[4] - quantile(DOTideData[, i], na.rm = TRUE)[2]))) #Calculate upper bound
  lowerBound <- 0 #Dissolved oxygen content cannot be less than zero
  DOTideData[, i][DOTideData[, i] > upperBound] <- NA #Replace any values higher than upper bound with NA
  DOTideData[, i][DOTideData[, i] < lowerBound] <- NA #Replace any values lower than upper bound with NA
} #Replace outliers with NA values

#Salinity Data
for(i in 4:nSites) { #For individual site data
  upperBound <- as.numeric((quantile(salinityTideData[, i], na.rm = TRUE)[4]) + (1.5*(quantile(salinityTideData[, i], na.rm = TRUE)[4] - quantile(salinityTideData[, i], na.rm = TRUE)[2]))) #Calculate upper bound
  lowerBound <- as.numeric((quantile(salinityTideData[, i], na.rm = TRUE)[2]) - (1.5*(quantile(salinityTideData[, i], na.rm = TRUE)[4] - quantile(salinityTideData[, i], na.rm = TRUE)[2]))) #Calculate lower bound
  salinityTideData[, i][salinityTideData[, i] > upperBound] <- NA #Replace any values higher than upper bound with NA
  salinityTideData[, i][salinityTideData[, i] < lowerBound] <- NA #Replace any values lower than upper bound with NA
} #Replace outliers with NA values

#### WRITE OUT AS NEW DATAFRAMES ####

#write.csv(pHTideData, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality3-Control/2017-12-18-pH-Data-QC-with-Tide-Data.csv") #pH
#write.csv(DOTideData, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-DO-Data-QC-with-Tide-Data.csv") #DO
#write.csv(salinityTideData, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Salinity-Data-QC-with-Tide-Data.csv") #Salinity