#In this script I will analyze my histology data.

#### SET WORKING DIRECTORY ####

setwd()
getwd()

#### IMPORT DATA####

histologyData <- read.csv("Documents/project-oyster-oa/data/Manchester/2018-02-27-Gigas-Histology-Classification.csv")
head(histologyData)

#### MATURATION STAGE ####

mature.stage <- 3

histologyData$Mature <- rep(0, nrow(histologyData))
histologyData$Mature[which(histologyData$Stage>=mature.stage)]<-rep(1, length(which(histologyData$Stage>=mature.stage)))
histologyData$Treatment <- c(rep("Ambient", times = 20), rep("Low", times = 10), rep("Ambient", times = 10))

mature.glm <- glm(Mature ~ Treatment + Pre.or.Post.OA, data = histologyData, family = binomial)
summary(mature.glm)

#### SEX RATIO ####

#contingency table...? try both with pre and without pre data