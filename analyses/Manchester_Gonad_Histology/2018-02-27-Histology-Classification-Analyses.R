#In this script I will analyze my gonad histology data. I will first see if treatment affected maturation state. Then, I will see if the sex ratios between treatments are homogenous.

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
mature.glm1 <- glm(Mature ~ factor(Treatment), family = binomial(link = "logit"), data = histologyData) #Ambient vs. low pH
anova(mature.glm1)
1-pf(3.7926/(30.024/38), 1, 38) #0.03466119 (Use Deviance/(ResDev/ResDF) to find F-value)
mature.glm2 <- glm(Mature ~ factor(modifiedSex), family = binomial(link = "logit"), data = histologyData) #Female vs. male vs. unripe
anova(mature.glm2)
1-pf(13.175/(20.641/37), 2, 37) #2.456105e-07. modifiedSex is the most significant, so this is the base model
mature.glm3 <- glm(Mature ~ Ferrous.inclusion.presence, family = binomial(link = "logit"), data = histologyData) #Ferrous inclusion vs. no ferrous inclusions
anova(mature.glm3)
1-pf(0.24821/(33.569/38), 1, 38) #0.5991481
mature.glm4 <- glm(Mature ~ Pre.or.Post.OA, family = binomial(link = "logit"), data = histologyData) #Pre vs. post sampling
anova(mature.glm4)
1-pf(0.79731/(33.019/38), 1, 38) #0.3441647

#Use add1 to find next significant variable
add1(mature.glm2, ~. + factor(Treatment) + Ferrous.inclusion.presence, test = "F", data = histologyData) #Neither variable is significant, so none will be included. modifiedSex is the only significant predictor

#See where differences are in mature.glm2

summary(mature.glm2) #Males are more mature than females

#### VISUALIZE MATURITY DATA ####
immatureLow <- length(which(histologyData$Treatment == "Low" & histologyData$Mature == 0)) #Count the number of oysters that are immature and were exposed to low pH
matureLow <- length(which(histologyData$Treatment == "Low" & histologyData$Mature == 1)) #Count the number of oysters that are mature and were exposed to low pH
immatureAmb <- length(which(histologyData$Treatment == "Ambient" & histologyData$Mature == 0 & histologyData$Pre.or.Post.OA == "Post")) #Count the number of oysters post-sampling that are immature and were exposed to ambient pH
matureAmb <- length(which(histologyData$Treatment == "Ambient" & histologyData$Mature == 1 & histologyData$Pre.or.Post.OA == "Post")) #Count the number of oysters post-sampling that are immature and were exposed to ambient pH

maturityFrequency <- data.frame("treatment" = c("Low", "Low", "Ambient", "Ambient"),
           "maturity" = c("Immature", "Mature", "Immature", "Mature"),
           "count" = c(immatureLow, matureLow, immatureAmb, matureAmb)) #Create a new dataframe
is.numeric(maturityFrequency$count)

#jpeg("2018-03-07-Maturity-by-Treatment-Barchart.jpeg", height = 1000, width = 1000)
barplot(height = maturityFrequency$count, xlab = "Treatment", ylab = "Number of Individuals", main = "Maturity by Treatment", col = "white") #Preliminary bar chart
#dev.off()

##### SEX RATIO #####
#Chi-squared test of homogeneity using only post-treatment sex classifications

#### CALCULATE CONTINGENCY TABLE VALUES ####
maleLow <- length(which(histologyData$Treatment == "Low" & histologyData$modifiedSex == "male")) #Count the number of oysters that are male and were exposed to low pH
femaleLow <- length(which(histologyData$Treatment == "Low" & histologyData$modifiedSex == "female"))
unripeLow <- length(which(histologyData$Treatment == "Low" & histologyData$modifiedSex == "unripe"))
maleAmb <- length(which(histologyData$Treatment == "Ambient" & histologyData$modifiedSex == "male" & histologyData$Pre.or.Post.OA == "Post")) #Include specification for post-treamtent samples only
femaleAmb <- length(which(histologyData$Treatment == "Ambient" & histologyData$modifiedSex == "female" & histologyData$Pre.or.Post.OA == "Post")) #Include specification for post-treamtent samples only
unripeAmb <- length(which(histologyData$Treatment == "Ambient" & histologyData$modifiedSex == "unripe" & histologyData$Pre.or.Post.OA == "Post")) #Include specification for post-treamtent samples only

#### CREATE CONTINGENCY TABLE ####
sexRatioContingencyTable <- data.frame("count" = c(maleLow, femaleLow, unripeLow, maleAmb, femaleAmb, unripeAmb),
                                       "row" = c(rep(1, times = 3), rep(2, times = 3)),
                                       "column" = c(1, 2, 3, 1, 2, 3)) #Make contingency table where rows specify pH treatment and columns specify sex classification
head(sexRatioContingencyTable) #Confirm table creation
r <- as.factor(sexRatioContingencyTable$row) #Recognize as factor
c <- as.factor(sexRatioContingencyTable$column) #Recognize as factor

#### TEST INDEPENDENT VARIABLES ####
ratio.glm1 <- glm(count ~ r + c, family = poisson(link = "log"), data = sexRatioContingencyTable) #Create a poisson GLM with a log link
anova(ratio.glm1)
1-pchisq(3.2779, 2) #0.1941838, Insignificant, so model fits. Testing interaction will leave df = 0. Sex ratios are homogenous between low and ambient pH treatments.