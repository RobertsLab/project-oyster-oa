#In this script, I'll identify any significant differences in reproductive output between oysters exposed to low and ambient pH treatments. Specifically, I'll look at the number of eggs produced the larval hatch rate, and the total number of larvae produced.

#### SET WORKING DIRECTORY ####

getwd()
setwd("Documents/project-oyster-oa/") #Set working directory as repository

##### EGG PRODUCTION #####

#### IMPORT DATA ####

eggProduction <- read.csv("data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2018-02-14-Egg-Production-Data.csv", header = TRUE) #Import egg production data
head(eggProduction) #Confirm import
colnames(eggProduction) #Get column names. I don't need anything after Female.Egg.Count
eggProduction <- eggProduction[,c(1:9)] #Retain necessary columns.
colnames(eggProduction) #Confirm change

#### CORRECT FOR NUMBER OF FEMALES THAT CONTRIBUTED ####
#22 females contributed to the low pool, 26 to the ambient, and 6 to the heat shock

correctedEggProduction <- data.frame("FemalePool" = c("Low", "Ambient", "HeatShock"),
                                     "EggCount1" = c(((eggProduction$Sample.Egg.Count.1[1])/22), ((eggProduction$Sample.Egg.Count.1[25])/26), ((eggProduction$Sample.Egg.Count.1[49])/6)),
                                     "EggCount2" = c(((eggProduction$Sample.Egg.Count.2[1])/22), ((eggProduction$Sample.Egg.Count.2[25])/26), ((eggProduction$Sample.Egg.Count.2[49])/6)),
                                     "EggCount3" = c(((eggProduction$Sample.Egg.Count.3[1])/22), ((eggProduction$Sample.Egg.Count.3[25])/26), ((eggProduction$Sample.Egg.Count.3[49])/6))) #Create a new dataframe with corrected data
head(correctedEggProduction) #Confirm creation

rownames(correctedEggProduction) <- correctedEggProduction$FemalePool #Set pool as row names
correctedEggProduction <- correctedEggProduction[, -1] #Remove Female Pool column
head(correctedEggProduction) #Confirm change

correctedEggProduction <- data.frame(t(correctedEggProduction)) #Transpose dataframe and maintain dataframe structure
head(correctedEggProduction) #Confirm change

#### CALCULATE EGG COUNT STANDARD DEVIATION ####

lowSD <- sqrt(var(correctedEggProduction$Low)) #SD = 278276.9
ambSD <- sqrt(var(correctedEggProduction$Ambient)) #SD = 114210.1
hsSD <- sqrt(var(correctedEggProduction$HeatShock)) #SD = 210323.8

#These standard deviations are large...

#### CALCULATE AVERAGE EGG COUNTS ####

lowAverage <- mean(correctedEggProduction$Low) #1190909
ambAverage <- mean(correctedEggProduction$Ambient) #1256410
hsAverage <- mean(correctedEggProduction$HeatShock) #2300000

#Heat shock animals produced most eggs on average, low pH the least

#### ADD TO TABLE ####

correctedEggProduction <- data.frame(t(correctedEggProduction)) #Transpose dataframe and maintain dataframe structure
head(correctedEggProduction) #Confirm change

correctedEggProduction$AverageEggCount <- c(lowAverage, ambAverage, hsAverage) #Add average counts
correctedEggProduction$StandardDeviation <- c(lowSD, ambSD, hsSD) #Add standard deviations
head(correctedEggProduction) #Confirm changes

#### ANOVA ####

correctedEggProduction <- data.frame(t(correctedEggProduction)) #Transpose dataframe and maintain dataframe structure
head(correctedEggProduction) #Confirm change

eggProductionANOVAData <- data.frame("Treatment" = c(rep("Low", times = 3), rep("Ambient", times = 3), rep("HeatShock", times =3)),
                                     "EggCount" = c(correctedEggProduction$Low[1:3], correctedEggProduction$Ambient[1:3], correctedEggProduction$HeatShock[1:3])) #Create new dataframe with only treatment and egg count columns
head(eggProductionANOVAData) #Confirm dataframe creation

treatmentANOVA <- aov(EggCount ~ Treatment, data = eggProductionANOVAData) #One-way ANOVA by female treatment
summary(treatmentANOVA)[[1]][["F value"]][[1]] #F = 25.87017
summary(treatmentANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.00112206
TukeyHSD(treatmentANOVA) #Significant differences are between Heat Shock and pH treatments (HS-A = 0.0022743; HS-L = 0.0016528)

##### LARVAL HATCH RATE #####

#### IMPORT DATA ####

hatchRate <- read.csv("data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2018-02-14-Hatch-Rate-Data.csv", header = TRUE) #Import hatch rate data
head(hatchRate) #Confirm import

#### ANOVA ####

hatchRateANOVA <- aov(Average.Hatch.Rate ~ Parental.Treatment, data = hatchRate) #One-way ANOVA by parental treatment

##### TOTAL NUMBER OF LARVAE #####