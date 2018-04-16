#In this script, I'll identify any significant differences in reproductive output between oysters exposed to low and ambient pH treatments. Specifically, I'll look at the number of eggs produced and the larval hatch rate.

#### SET WORKING DIRECTORY ####

getwd()
setwd("../../") #Set working directory as repository

##### EGG PRODUCTION #####
#If there were differences in gonad maturation between females exposed to stressors and those kept at ambient temperatures, then a different amount of eggs would be produced.

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

#### VISUALIZE DATA ####

correctedEggProduction <- data.frame(t(correctedEggProduction)) #Transpose dataframe and maintain dataframe structure
head(correctedEggProduction) #Confirm change
eggProductionANOVAData <- data.frame("Treatment" = c(rep("Low", times = 3), rep("Ambient", times = 3), rep("HeatShock", times =3)),
                                     "EggCount" = c(correctedEggProduction$Low[1:3], correctedEggProduction$Ambient[1:3], correctedEggProduction$HeatShock[1:3])) #Create new dataframe with only treatment and egg count columns
head(eggProductionANOVAData) #Confirm dataframe creation

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-02-14-Egg-Production-by-Treatment.jpeg", width = 1500, height = 1000)
plot(x = eggProductionANOVAData$Treatment, y = eggProductionANOVAData$EggCount, xlab = "Treatment", ylab = "Egg Count", main = "Egg Production by Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
#dev.off()

#### VISUALIZE JUST THE PH TREATMENT DATA ####

tail(eggProductionANOVAData) #Need to remove heat shock data, or rows 7-9
eggProductionpHOnly <- eggProductionANOVAData[-c(7:9),] #Remove rows 7-9
tail(eggProductionpHOnly) #Confirm removal
eggProductionpHOnly$Treatment <- factor(eggProductionpHOnly$Treatment) #Make sure residual factors are no longer present

pHTreatmentANOVA <- aov(EggCount ~ Treatment, data = eggProductionpHOnly) #One-way ANOVA by female treatment
sqrt(summary(pHTreatmentANOVA)[[1]][["F value"]][[1]]) #t = 0.3771626
summary(pHTreatmentANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.7252096

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-03-07-Egg-Production-by-pH-Treatment.jpeg", width = 1500, height = 1000)
plot(x = eggProductionpHOnly$Treatment, y = eggProductionpHOnly$EggCount, xlab = "Treatment", ylab = "Egg Count", main = "Egg Production by Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
legend("topleft", bty = "n", legend = paste("t =", format(sqrt(summary(pHTreatmentANOVA)[[1]][["F value"]][[1]]), digits = 4), "p =", format(summary(pHTreatmentANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add t and p-value
#dev.off()

#### ANOVA ####

treatmentANOVA <- aov(EggCount ~ Treatment, data = eggProductionANOVAData) #One-way ANOVA by female treatment
summary(treatmentANOVA)[[1]][["F value"]][[1]] #F = 25.87017
summary(treatmentANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.00112206
TukeyHSD(treatmentANOVA) #Significant differences are between Heat Shock and pH treatments (HS-A = 0.0022743; L-HS = 0.0016528)

##### LARVAL HATCH RATE #####
#Hatch rate is my proxy for the number of larvae produced.

#### IMPORT DATA ####

hatchRate <- read.csv("data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2018-02-14-Hatch-Rate-Data.csv", header = TRUE) #Import hatch rate data
head(hatchRate) #Confirm import

#### VISUALIZE DATA ####

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-02-14-Hatch-Rate-by-Treatment.jpeg", width = 1500, height = 1000)
plot(x = hatchRate$Parental.Treatment, y = hatchRate$Average.Hatch.Rate, xlab = "Treatment", ylab = "Hatch Rate", main = "Hatch Rate by Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
#dev.off()

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-02-14-Hatch-Rate-by-Female-Treatment.jpeg", width = 1500, height = 1000)
plot(x = hatchRate$Female.Treatment, y = hatchRate$Average.Hatch.Rate, xlab = "Treatment", ylab = "Hatch Rate", main = "Hatch Rate by Female Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
#dev.off()

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-02-14-Hatch-Rate-by-Male-Treatment.jpeg", width = 1500, height = 1000)
plot(x = hatchRate$Male.Treatment, y = hatchRate$Average.Hatch.Rate, xlab = "Treatment", ylab = "Hatch Rate", main = "Hatch Rate by Male Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
#dev.off()

#### VISUALIZE JUST THE PH TREATMENT DATA ####

hatchRatepHOnly <- hatchRate[-c(25:26),] #Remove heat shock data
tail(hatchRatepHOnly) #Confirm removal
hatchRatepHOnly$Parental.Treatment <- factor(hatchRatepHOnly$Parental.Treatment) #Make sure residual factors are no longer present
hatchRatepHOnly$Female.Treatment <- factor(hatchRatepHOnly$Female.Treatment) #Make sure residual factors are no longer present
hatchRatepHOnly$Male.Treatment <- factor(hatchRatepHOnly$Male.Treatment) #Make sure residual factors are no longer present

hatchRatepHOnly <- hatchRatepHOnly[-24,] #Remove outlier from Ambient-Ambient group
tail(hatchRatepHOnly) #Confirm removal

hatchRatepHTreatmentANOVA <- aov(Average.Hatch.Rate ~ Parental.Treatment, data = hatchRatepHOnly) #One-way ANOVA by parental treatment
summary(hatchRatepHTreatmentANOVA)[[1]][["F value"]][[1]] #F = 3.10953
summary(hatchRatepHTreatmentANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.05082859

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-03-07-Hatch-Rate-by-pH-Treatment-All-Groups.jpeg", width = 1500, height = 1000)
plot(x = hatchRatepHOnly$Parental.Treatment, y = hatchRatepHOnly$Average.Hatch.Rate, xlab = "Treatment", ylab = "Hatch Rate", main = "Hatch Rate by Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
legend("topright", bty = "n", legend = paste("t =", format(sqrt(summary(hatchRatepHTreatmentANOVA)[[1]][["F value"]][[1]]), digits = 4), "p =", format(summary(hatchRatepHTreatmentANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add t and p-value
#dev.off()

hatchRatepHTreatmentFemaleANOVA <- aov(Average.Hatch.Rate ~ Female.Treatment, data = hatchRatepHOnly) #One-way ANOVA by female treatment
sqrt(summary(hatchRatepHTreatmentFemaleANOVA)[[1]][["F value"]][[1]]) #t = 2.99445
summary(hatchRatepHTreatmentFemaleANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.006908886

#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-03-07-Hatch-Rate-by-pH-Treatment-Female-Groups.jpeg", width = 1500, height = 1000)
plot(x = hatchRatepHOnly$Female.Treatment, y = hatchRatepHOnly$Average.Hatch.Rate, xlab = "Treatment", ylab = "Hatch Rate", main = "Hatch Rate by Female Treatment", cex.main = 4, cex.axis = 1, cex.lab = 1.4) #Preliminary plot
legend("topright", bty = "n", legend = paste("t =", format(sqrt(summary(hatchRatepHTreatmentFemaleANOVA)[[1]][["F value"]][[1]]), digits = 4), "p =", format(summary(hatchRatepHTreatmentFemaleANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add t and p-value
#dev.off()

#### ANOVA ####

hatchRateANOVA <- aov(Average.Hatch.Rate ~ Parental.Treatment, data = hatchRate) #One-way ANOVA by parental treatment
summary(hatchRateANOVA)[[1]][["F value"]][[1]] #F = 2.559579
summary(hatchRateANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.06864339

hatchRateFemaleANOVA <- aov(Average.Hatch.Rate ~ Female.Treatment, data = hatchRate) #One-way ANOVA by female treatment
summary(hatchRateFemaleANOVA)[[1]][["F value"]][[1]] #F = 5.606537
summary(hatchRateFemaleANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.01039109
TukeyHSD(hatchRateFemaleANOVA, method = "BH") #Significant difference between low and ambient pH hatch rates, with low pH hatch rates less than ambient pH hatch rates (L-A = 0.0111872)

hatchRateMaleANOVA <- aov(Average.Hatch.Rate ~ Male.Treatment, data = hatchRate) #One-way ANOVA by male treatment
summary(hatchRateMaleANOVA)[[1]][["F value"]][[1]] #F = 0.3819449
summary(hatchRateMaleANOVA)[[1]][["Pr(>F)"]][[1]] #p = 0.6867814

#### MANCHESTER PAPER FIGURES ####
#jpeg(filename = "analyses/Manchester_ReproductiveOutput_20180214/2018-04-16-Manchester-Paper-Figure.jpeg", width = 1500, height = 1000)
plot(x = hatchRatepHOnly$Parental.Treatment, y = hatchRatepHOnly$Average.Hatch.Rate, cex.axis = 2) #Preliminary plot. Will modify in InDesign for publication
#dev.off()