#In this script, I'll make individual NMDS plots for proteins based on functions: oxidative stress (4 proteins), heat shock (1 protein), acid-base balance (1 protein), drug resistance (1 protein), fatty acid metabolism (1 protein), carbohydrate metabolism (2 proteins), cell growth and maintenance (5 proteins) 

#### SET WORKING DIRECTOR ####

setwd("../..") #Set directory to master SRM folder
getwd()

#### IMPORT DATA ####

SRMAveragedAreas <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv") #Import averaged areas
rownames(SRMAveragedAreas) <- SRMAveragedAreas$X #Set rownames
SRMAveragedAreas <- SRMAveragedAreas[,-1] #Removed rownames column
head(SRMAveragedAreas) #Confirm import

sampleColumnNames <- c("O01", "O04", "O08", "O10", "O100", "O101", "O102", "O106", "O118", "O121", "O124", "O131", "O137", "O140", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O51", "O56", "O60", "O64", "O66", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector

biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
colnames(biologicalReplicates) <- c("Sample.Number", "Site", "Eelgrass.Condition") #Rename columns
head(biologicalReplicates) #Confirm import

#### SEPARATE DATA BY PROTEIN FUNCTION ####

oxidativeStressAreas <- SRMAveragedAreas[c(4:5, 6:8, 14:16, 17),] #Subset oxidative stress areas
heatShockAreas <- SRMAveragedAreas[c(20:21), ] #Subset heat shock areas
acidBaseAreas <- SRMAveragedAreas[c(11:13),] #Subset acid base balance areas
drugResistanceAreas <- SRMAveragedAreas[c(27:28),] #Subset drug resistance areas
fattyAcidAreas <- SRMAveragedAreas[c(1:3),] #Subset fatty acid metabolism areas
carbohydrateAreas <- SRMAveragedAreas[c(18:19, 24:26),] #Subset carbohydrate metabolism areas
growthMaintenanceAreas <- SRMAveragedAreas[c(9:10, 22:23, 29:31, 32:34, 35:37),] #Subset cell growth and maintenance areas

#### LOAD NMDS DEPENDENCIES ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS. It can be found the project-oyster-oa repo >> analyses >> DNR_Preliminary_Analyses_20170321. It is also at the bottom of this script.
install.packages("vegan") #Install vegan package
library(vegan)

#### ASSIGN COLORS AND SHAPES ####
#Because colors and shapes are related to sample IDs, I can use the same color scheme regardless of which protein function is displayed

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(Sample.Number = sampleColumnNames,
                            y = rep(x = 0, times = length(sampleColumnNames))) #Create a temporary dataframe with sample  names
head(temporaryData) #Confirm dataframe creation
NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "Sample.Number") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[,-2] #Remove empty column
head(NMDSColorShapeCustomization) #Confirm removal

#Create a color and shape palette
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Site),] #Reorder so sites are sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Colors <- c(rep(x = "red", times = sum(NMDSColorShapeCustomization$Site == "CI")),
                 rep(x = "blue", times = sum(NMDSColorShapeCustomization$Site == "FB")),
                 rep(x = "magenta", times = sum(NMDSColorShapeCustomization$Site == "PG")),
                 rep(x = "green", times = sum(NMDSColorShapeCustomization$Site == "SK")),
                 rep(x = "black", times = sum(NMDSColorShapeCustomization$Site == "WB"))) #Create a color vector
NMDSColorShapeCustomization[,4] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Eelgrass.Condition),] #Reorder so eelgrass condition is sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Shapes <- c(rep(x = 16, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Bare")),
                 rep(x = 17, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Eelgrass"))) #Make a shape vector
NMDSColorShapeCustomization[,5] <- NMDS.Shapes #Add the shape vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Sample.Number),] #Resort by sample number
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
colnames(NMDSColorShapeCustomization) <- c("Sample.Number", "Site", "Eelgrass.Condition", "Color", "Shape") #Change column names
head(NMDSColorShapeCustomization) #Confirm change

#### PREPARE FOR ANOSIM ####
#This code makes the table of factor needed for my ANOSIM

ANOSIMReplicates <- biologicalReplicates #Subset sample numbers used as IDs in ANOSIM
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
ANOSIMReplicates$Site.Eelgrass <- paste(ANOSIMReplicates$Site, ANOSIMReplicates$Eelgrass.Condition)
head(ANOSIMReplicates) #Confirm changes

str(ANOSIMReplicates) #Examine structure
ANOSIMReplicates$Site <- factor(ANOSIMReplicates$Site) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Eelgrass.Condition <- factor(ANOSIMReplicates$Eelgrass.Condition) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Site.Eelgrass <- factor(ANOSIMReplicates$Site.Eelgrass) #Make sure only preesnt factors are recognized
str(ANOSIMReplicates) #Confirm structure

#### OXIDATIVE STRESS ####

#Make NMDS plot

area.protID <- oxidativeStressAreas #Duplicate dataframe
head(area.protID) #Confirm changes
area.t <- t(area.protID) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.oxidativestress.averaged.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.oxidativestress.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.oxidativestress.averaged.euclidean <- envfit(proc.nmds.oxidativestress.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.oxidativestress.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.oxidativestress.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged-Oxidative-Stress.jpeg", width = 1000, height = 750)
fig.nmds.oxidativestress <- ordiplot(proc.nmds.oxidativestress.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds.oxidativestress, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#ANOSIM

dissimArea.t <- vegdist(area.t, "euclidean") #Calculate dissimilarity matrix
siteANOSIMOxidativeStress <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIMOxidativeStress)
plot(siteANOSIMOxidativeStress)
simper(proc.nmds.oxidativestress.averaged.euclidean, ANOSIMReplicates$Site)

eelgrassANOSIMOxidativeStress <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIMOxidativeStress)
plot(eelgrassANOSIMOxidativeStress)
simper(proc.nmds.oxidativestress.averaged.euclidean, ANOSIMReplicates$Eelgrass.Condition)

siteEelgrassANOSIMOxidativeStress <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIMOxidativeStress)
plot(siteEelgrassANOSIMOxidativeStress)
simper(proc.nmds.oxidativestress.averaged.euclidean, ANOSIMReplicates$Site.Eelgrass)

#### HEAT SHOCK ####

#Make NMDS plot

area.protID2 <- heatShockAreas #Duplicate dataframe
head(area.protID2) #Confirm changes
area2.t <- t(area.protID2) #Transpose the file so that rows and columns are switched
head(area2.t) #Confirm transposition
area2.tra <- (area2.t+1) #Add 1 to all values before transforming
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.heatshock.averaged.euclidean <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.heatshock.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.heatshock.averaged.euclidean <- envfit(proc.nmds.heatshock.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.heatshock.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.heatshock.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged-HeatShock.jpeg", width = 1000, height = 750)
fig.nmds.heatshock <- ordiplot(proc.nmds.heatshock.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds.heatshock, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#ANOSIM

dissimArea2.t <- vegdist(area2.t, "euclidean") #Calculate dissimilarity matrix
siteANOSIMHeatShock <- anosim(dat = dissimArea2.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIMHeatShock)
plot(siteANOSIMHeatShock)
simper(proc.nmds.heatshock.averaged.euclidean, ANOSIMReplicates$Site)

eelgrassANOSIMHeatShock <- anosim(dat = dissimArea2.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIMHeatShock)
plot(eelgrassANOSIMHeatShock)
simper(proc.nmds.heatshock.averaged.euclidean, ANOSIMReplicates$Eelgrass.Condition)

siteEelgrassANOSIMHeatShock <- anosim(dat = dissimArea2.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIMHeatShock)
plot(siteEelgrassANOSIMHeatShock)
simper(proc.nmds.heatshock.averaged.euclidean, ANOSIMReplicates$Site.Eelgrass)

#### ACID-BASE BALANCE ####

#Make NMDS plot

area.protID3 <- acidBaseAreas #Duplicate dataframe
head(area.protID3) #Confirm changes
area3.t <- t(area.protID3) #Transpose the file so that rows and columns are switched
head(area3.t) #Confirm transposition
area3.tra <- (area3.t+1) #Add 1 to all values before transforming
area3.tra <- data.trans(area3.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.acidbase.averaged.euclidean <- metaMDS(area3.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.acidbase.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.acidbase.averaged.euclidean <- envfit(proc.nmds.heatshock.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.acidbase.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.acidbase.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged-AcidBase.jpeg", width = 1000, height = 750)
fig.nmds.acidbase <- ordiplot(proc.nmds.acidbase.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds.acidbase, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#ANOSIM

dissimArea3.t <- vegdist(area3.t, "euclidean") #Calculate dissimilarity matrix
siteANOSIMAcidBase <- anosim(dat = dissimArea3.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIMAcidBase)
plot(siteANOSIMAcidBase)

eelgrassANOSIMAcidBase <- anosim(dat = dissimArea3.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIMAcidBase)
plot(eelgrassANOSIMAcidBase)

siteEelgrassANOSIMAcidBase <- anosim(dat = dissimArea3.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIMAcidBase)
plot(siteEelgrassANOSIMAcidBase)

####
### DRUG RESISTANCE ####

#Make NMDS plot

area.protID4 <- drugResistanceAreas #Duplicate dataframe
head(area.protID4) #Confirm changes
area4.t <- t(area.protID4) #Transpose the file so that rows and columns are switched
head(area4.t) #Confirm transposition
area4.tra <- (area4.t+1) #Add 1 to all values before transforming
area4.tra <- data.trans(area4.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.drug.averaged.euclidean <- metaMDS(area4.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.drug.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.drug.averaged.euclidean <- envfit(proc.nmds.heatshock.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.drug.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.drug.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged-DrugResistance.jpeg", width = 1000, height = 750)
fig.nmds.drug <- ordiplot(proc.nmds.drug.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds.drug, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#ANOSIM

dissimArea4.t <- vegdist(area4.t, "euclidean") #Calculate dissimilarity matrix
siteANOSIMADrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIMADrug)
plot(siteANOSIMADrug)

eelgrassANOSIMDrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIMDrug)
plot(eelgrassANOSIMDrug)

siteEelgrassANOSIMADrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIMADrug)
plot(siteEelgrassANOSIMADrug)

#### FATTY ACID METABOLISM ####

#Make NMDS plot

area.protID4 <- fattyAcidAreas #Duplicate dataframe
head(area.protID4) #Confirm changes
area4.t <- t(area.protID4) #Transpose the file so that rows and columns are switched
head(area4.t) #Confirm transposition
area4.tra <- (area4.t+1) #Add 1 to all values before transforming
area4.tra <- data.trans(area4.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.drug.averaged.euclidean <- metaMDS(area4.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.drug.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.drug.averaged.euclidean <- envfit(proc.nmds.heatshock.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.drug.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.drug.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged-DrugResistance.jpeg", width = 1000, height = 750)
fig.nmds.drug <- ordiplot(proc.nmds.drug.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds.drug, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#ANOSIM

dissimArea4.t <- vegdist(area4.t, "euclidean") #Calculate dissimilarity matrix
siteANOSIMADrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIMADrug)
plot(siteANOSIMADrug)

eelgrassANOSIMDrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIMDrug)
plot(eelgrassANOSIMDrug)

siteEelgrassANOSIMADrug <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIMADrug)
plot(siteEelgrassANOSIMADrug)