#Steven went through my normalized data set of technical replicates and calculated a coefficient of variation for each transition. He removed any transitions that had CVs greater than 20. Using this dataset, I'll average my technical replicates and proceed with an NMDS and ANOSIM.

#### SET WORKING DIRECTORY ####

setwd("../..") #Change working directory to the master SRM folder, and not the folder where this script is hosted.
getwd() #Confirm changes

#### IMPORT DATA ####

SRMModifiedAreas <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-24-Norm-modSR.csv", header = TRUE) #Import Steven's modified dataset
head(SRMModifiedAreas) #Confirm import. There are six columns. The first column is null, and the last column, coefficient of variation, is not needed.
SRMModifiedAreas <- SRMModifiedAreas[,-c(1, 6)] #Remove unnecessary columns
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

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(Sample.Number = technicalReplicates,
                            y = rep(x = 0, times = length(technicalReplicates))) #Create a temporary dataframe with technical replicate names used in NMDS
head(temporaryData) #Confirm dataframe creation
NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "Sample.Number") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[-c(97:98),-2] #Remove OBLNK2 and empty column
tail(NMDSColorShapeCustomization) #Confirm removal
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[seq(from = 1, to = 95, by = 2),] #Keep only every other row
head(NMDSColorShapeCustomization) #Confirm changes
NMDSColorShapeCustomization$Sample.Number #Confirm changes