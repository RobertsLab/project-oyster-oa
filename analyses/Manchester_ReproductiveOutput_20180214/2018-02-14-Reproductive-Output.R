#In this script, I'll identify any significant differences in reproductive output between oysters exposed to low and ambient pH treatments. Specifically, I'll look at the number of eggs produced and the larval hatch rate.

#### SET WORKING DIRECTORY ####

getwd()
setwd("Documents/project-oyster-oa/") #Set working directory as repository

#### EGG PRODUCTION ####

#### IMPORT DATA ####

eggProduction <- read.csv("data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2018-02-14-Egg-Production-Data.csv", header = TRUE) #Import egg production data
head(eggProduction) #Confirm import
colnames(eggProduction) #Get column names. I don't need anything after Female.Egg.Count
eggProduction <- eggProduction[,c(1:7)] #Retain necessary columns.
colnames(eggProduction) #Confirm change

#### CORRECT FOR NUMBER OF FEMALES THAT CONTRIBUTED ####
#22 females contributed to the low pool, 26 to the ambient, and 6 to the heat shock

correctedLowContribution <- (eggProduction$Female.Egg.Count..eggs.[1])/22 #1190909 eggs per female
correctedAmbientContribution <- (eggProduction$Female.Egg.Count..eggs.[25])/26 #1256410 eggs per female
correctedHeatShockContribution <- (eggProduction$Female.Egg.Count..eggs.[49])/6 #2300000 eggs per female

max(correctedAmbientContribution, correctedHeatShockContribution, correctedLowContribution) #Heat shock females produced more eggs on average
min(correctedAmbientContribution, correctedHeatShockContribution, correctedLowContribution) #Low pH  females produced less eggs on average
