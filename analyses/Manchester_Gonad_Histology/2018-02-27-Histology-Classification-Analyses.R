#In this script I will analyze my histology data.

##### SET WORKING DIRECTORY #####

getwd() #Working directory is set to the Manchester_Gonad_Histology folder

##### IMPORT DATA #####

histologyData <- read.csv("../../data/Manchester/2018-02-27-Gigas-Histology-Classification.csv") #Import histology data
head(histologyData) #Confirm import
histologyData <- histologyData[, -6] #Remove Notes column

##### MATURATION STAGE #####

#### REFORMAT DATA ####
mature.stage <- 3 #Set the maturation stage to be anything that is ripe and spawning
histologyData$Mature <- rep(0, nrow(histologyData)) #Create a new column for Yes/No information on maturity. 0 is immature, 1 is mature
histologyData$Mature[which(histologyData$Stage >= mature.stage)] <- rep(1, length(which(histologyData$Stage >= mature.stage))) #For anything where the stage is greater than or equal to the mature.stage, replace it with a 1
head(histologyData) #Confirm changes

histologyData$Treatment <- c(rep("Ambient", times = 20), rep("Low", times = 10), rep("Ambient", times = 10)) #Create a treatment column
head(histologyData) #Confirm changes

histologyData$modifiedSex <- rep(0, nrow(histologyData)) #Create a new column to modify sex classifications
histologyData$modifiedSex[which(histologyData$Sex == "N/A")] <- rep("unripe", length(which(histologyData$Sex == "N/A"))) #For anything where sex is N/A, replace 0 with "unripe"
histologyData$modifiedSex[which(histologyData$Sex == "female")] <- rep("female", length(which(histologyData$Sex == "female"))) #For anything where sex is female, replace 0 with "female"
histologyData$modifiedSex[which(histologyData$Sex == "male")] <- rep("male", length(which(histologyData$Sex == "male"))) #For anything where sex is male, replace 0 with "male"
head(histologyData) #Confirm changes

#### STEPWISE ADDITION MODEL BUILDING ####

#Find first significant variable using a binomial GLM and cannonical logit link
glm1 <- glm(Mature ~ factor(Treatment), family = binomial(link = "logit"), data = histologyData) #Ambient vs. low pH
anova(glm1)
1-pf(1.7886/(48.658/38), 1, 38) #0.244599
glm2 <- glm(Mature ~ factor(modifiedSex), family = binomial(link = "logit"), data = histologyData) #Female vs. male vs. unripe
anova(glm2)
1-pf(13.9/(36.547/37), 2, 37) #2.85091e-05. modifiedSex is the most significant, so this is the base model
glm3 <- glm(Mature ~ Ferrous.inclusion.presence, family = binomial(link = "logit"), data = histologyData) #Ferrous inclusion vs. no ferrous inclusions
anova(glm3)
1-pf(1.3318/(49.115/38), 1, 38) #0.3164832

#Use add1


#### SEX RATIO ####

#contingency table...? try both with pre and without pre data