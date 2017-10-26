#Steven went through my normalized data set of technical replicates and calculated a coefficient of variation for each transition. He removed any transitions that had CVs greater than 20. Using this dataset, filter down further to CV ≤ 15, and then average my technical replicates and proceed with an NMDS and ANOSIM.

#### SET WORKING DIRECTORY ####

setwd("../../..") #Change working directory to the master SRM folder, and not the folder where this script is hosted.
getwd() #Confirm changes

#### IMPORT DATA ####

SRMModifiedAreas <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-24-Norm-modSR.csv", header = TRUE) #Import Steven's modified dataset
head(SRMModifiedAreas) #Confirm import
SRMModifiedAreas <- subset(SRMModifiedAreas, subset = SRMModifiedAreas$CoV <= 15) #Only keep rows with coefficient of variation ≤ 15
max(SRMModifiedAreas$CoV) <= 15 #Statement should be TRUE if maximum does not exceed 15
SRMModifiedAreas <- SRMModifiedAreas[,-c(1, 6)] #Remove unnecessary columns (first column and CoV column)
head(SRMModifiedAreas) #Confirm changes. Now I only have Protein.Name, Sample, Replicate1 and Replicate2 (peak areas from Skyline, which are a proxy for protein abundance)

#### AVERAGE TECHNICAL REPLICATES ####

SRMAveragedAreas <- SRMModifiedAreas #Dupliate dataframe
SRMAveragedAreas$Average.Area <- ((SRMAveragedAreas$Replicate1 + SRMAveragedAreas$Replicate2)/2) #Average peak areas and save as a new column
SRMAveragedAreas <- SRMAveragedAreas[, -c(3:4)] #Remove replicate area columns
head(SRMAveragedAreas) #Confirm changes

#### REFORMAT DATAFRAME FOR NMDS ####

#The dataframe I use needs to have sample numbers as columns, protein name as rownames, and area data in the middle. This is called casting.
library(reshape2) #Install package to pivot table
SRMAveragedAreasPivoted <- dcast(SRMAveragedAreas, Protein.Name ~ Sample) #Cast table!
head(SRMAveragedAreasPivoted) #Confirm table cast!

SRMAveragedAreasPivotedCorrected <- SRMAveragedAreasPivoted #Duplicate dataframe
SRMAveragedAreasPivotedCorrected[is.na(SRMAveragedAreasPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMAveragedAreasPivotedCorrected) #Confirm there are no NAs

rownames(SRMAveragedAreasPivotedCorrected) <- SRMAveragedAreasPivotedCorrected$Protein.Name #Save protein name column as rownames
SRMAveragedAreasPivotedCorrected <- SRMAveragedAreasPivotedCorrected[, -1] #Remove Protein.Name column
head(SRMAveragedAreasPivotedCorrected) #Confirm changes

#### MAKE BASE NMDS PLOT ####

source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

area.protID <- SRMAveragedAreasPivotedCorrected #Save all area data as a new dataframe
head(area.protID) #Confirm change

area.t <- t(area.protID) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.averaged.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
#stressplot(proc.nmds.norm.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.norm.averaged.euclidean <- envfit(proc.nmds.norm.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.norm.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### ASSIGN COLORS AND SHAPES ####

sampleIDs <- c("O01", "O04", "O06", "O08", "O10", "O100", "O101", "O102", "O103", "O106", "O118", "O121", "O122", "O124", "O128", "O131", "O137", "O14", "O140", "O145", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O49", "O51", "O52", "O56", "O60", "O64", "O66", "O71", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector

biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
head(biologicalReplicates) #Confirm import

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(sampleIDs = sampleIDs,
                            y = rep(x = 0, times = length(sampleIDs))) #Create a temporary dataframe with sample IDs used in NMDS
head(temporaryData) #Confirm dataframe creation

NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "sampleIDs") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[, -2] #Remove empty column
tail(NMDSColorShapeCustomization) #Confirm removal
NMDSColorShapeCustomization$sampleIDs #Confirm all sample IDs are there

#Create a color and shape palette
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(site),] #Reorder so sites are sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Colors <- c(rep(x = "red", times = sum(NMDSColorShapeCustomization$site == "CI")),
                 rep(x = "blue", times = sum(NMDSColorShapeCustomization$site == "FB")),
                 rep(x = "black", times = sum(NMDSColorShapeCustomization$site == "PG")),
                 rep(x = "green", times = sum(NMDSColorShapeCustomization$site == "SK")),
                 rep(x = "magenta", times = sum(NMDSColorShapeCustomization$site == "WB"))) #Create a color vector
NMDSColorShapeCustomization[,4] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(eelgrassCondition),] #Reorder so eelgrass condition is sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Shapes <- c(rep(x = 16, times = sum(NMDSColorShapeCustomization$eelgrassCondition == "Bare")),
                 rep(x = 17, times = sum(NMDSColorShapeCustomization$eelgrassCondition == "Eelgrass"))) #Make a shape vector
NMDSColorShapeCustomization[,5] <- NMDS.Shapes #Add the shape vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
colnames(NMDSColorShapeCustomization) <- c("Sample.Number", "Site", "Eelgrass.Condition", "Color", "Shape") #Change column names
head(NMDSColorShapeCustomization) #Confirm change

#### NMDS REFINEMENT ####

#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-NMDS-Norm-Analysis-Averaged-CV15-Filtered.jpeg", width = 1000, height = 1000)
fig.nmds <- ordiplot(proc.nmds.norm.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab= "Axis 1", ylab= "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", cex = .5, pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))
title("Protein Expression Similarities between Sites and Habitats")
#dev.off()

#### ANOSIM ####

dissimArea.t <- vegdist(area.t, "euclidean") #Calculate dissimilarity matrix
ANOSIMReplicates <- biologicalReplicates #Duplicate dataframe
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
head(ANOSIMReplicates) #Confirm changes

ANOSIMReplicates$eelgrassCondition <- factor(ANOSIMReplicates$eelgrassCondition) #Make sure residual factors are no longer present
ANOSIMReplicates$site <- factor(ANOSIMReplicates$site) #Make sure residual factors are no longer present
str(ANOSIMReplicates) #Confirm new factors

siteNormANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site presence
summary(siteNormANOSIM)
plot(siteNormANOSIM)
simper(proc.nmds.norm.averaged.euclidean, ANOSIMReplicates$site) #Error in rowSums(comm, na.rm = TRUE): 'x' must be an array of at least two dimensions

eelgrassNormANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassNormANOSIM)
plot(eelgrassNormANOSIM)
simper(proc.nmds.norm.averaged.euclidean, ANOSIMReplicates$eelgrassCondition) #Error in rowSums(comm, na.rm = TRUE): 'x' must be an array of at least two dimensions