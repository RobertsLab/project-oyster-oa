#In this script, I'll use an NMDS plot to see how my technical replication fared after Steven filtered out transitions with a coefficient of variation greater than 20.

#### SET WORKING DIRECTORY ####

setwd("../..") #Change working directory to the master SRM folder, and not the folder where this script is hosted.
getwd() #Confirm changes

#### IMPORT DATA ####

SRMModifiedAreas <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-24-Norm-modSR.csv", header = TRUE) #Import Steven's modified dataset
head(SRMModifiedAreas) #Confirm import. There are six columns. The first column is null, and the last column, coefficient of variation, is not needed.
SRMModifiedAreas <- SRMModifiedAreas[,-c(1, 6)] #Remove unnecessary columns
head(SRMModifiedAreas) #Confirm changes. Now I only have Protein.Name, Sample, Replicate1 and Replicate2 (peak areas from Skyline, which are a proxy for protein abundance)

#### REMFORMAT DATA ####

SRMModifiedAreasReplicate1 <- SRMModifiedAreas[, 1:3] #Duplicate dataframe, keeping only Replicate 1 information
SRMModifiedAreasReplicate1$Sample.Number <- paste(SRMModifiedAreasReplicate1$Sample, "-1") #Create a new column with sample and replicate number
SRMModifiedAreasReplicate1 <- SRMModifiedAreasReplicate1[, -2] #Remove Sample column
colnames(SRMModifiedAreasReplicate1) <- c("Protein.Name", "Area", "Sample") #Change column names
head(SRMModifiedAreasReplicate1) #Confirm changes

SRMModifiedAreasReplicate2 <- SRMModifiedAreas[, c(1, 2, 4)] #Duplicate dataframe, keeping only Replicate 2 information
SRMModifiedAreasReplicate2$Sample.Number <- paste(SRMModifiedAreasReplicate2$Sample, "-2") #Create a new column with sample and replicate number
SRMModifiedAreasReplicate2 <- SRMModifiedAreasReplicate2[, -2] #Remove Sample column
colnames(SRMModifiedAreasReplicate2) <- c("Protein.Name", "Area", "Sample") #Change column names
head(SRMModifiedAreasReplicate2) #Confirm changes

SRMModifiedAreasLong <- rbind(SRMModifiedAreasReplicate1, SRMModifiedAreasReplicate2) #Paste dataframes together
head(SRMModifiedAreasLong) #Confirm changes
tail(SRMModifiedAreasLong) #Confirm changes
transform(SRMModifiedAreasLong, Area = as.numeric(Area)) #Make sure Area is recognized as a numeric variable
is.numeric(SRMModifiedAreasLong$Area) #Confirm change

library(reshape2) #Install package to pivot table
SRMModifiedAreasPivoted <- dcast(SRMModifiedAreasLong, Protein.Name ~ Sample, value.var = "Area") #Cast table! Protein.Name remains as a column with Sample Number as column headers. Area column defined as values.
head(SRMModifiedAreasPivoted) #Confirm cast
rownames(SRMModifiedAreasPivoted) <- SRMModifiedAreasPivoted[, 1] #Save Protein.Name column as rownames
SRMModifiedAreasPivoted <- SRMModifiedAreasPivoted[, -1] #Remove Protein.Name column
head(SRMModifiedAreasPivoted) #Confirm changes

SRMModifiedAreasPivotedCorrected <- SRMModifiedAreasPivoted #Duplicate dataframe
SRMModifiedAreasPivotedCorrected[is.na(SRMModifiedAreasPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMModifiedAreasPivotedCorrected) #Confirm there are no NAs

#### LOAD FUNCTIONS NEEDED FOR NMDS ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste
install.packages("vegan") #Install vegan package
library(vegan)

#### MAKE NMDS PLOT ####

area.protID <- SRMModifiedAreasPivotedCorrected #Save all area data as a new dataframe
head(area.protID) #Confirm changes

area.t <- t(area.protID) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.norm.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.euclidean <- envfit(proc.nmds.nonnorm.euclidean$points, area.t, perm = 1000) #Calculate loadings
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-NMDS-TechnicalReplication-Normalized-after-CV-Filtering.jpeg", width = 1000, height = 1000) #Save plot
ordiplot(proc.nmds.norm.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.euclidean, p.max=.01, col='blue') #Plot eigenvectors
#dev.off() #Turn off plotting mechanism

proc.nmds.norm.euclidean.log <- metaMDS(area.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance
stressplot(proc.nmds.norm.euclidean.log) #Make Shepard plot
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-NMDS-TechnicalReplication-Normalized-LogTransformed-after-CV-Filtering.jpeg", width = 1000, height = 1000) #Save plot
ordiplot(proc.nmds.norm.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#My technical replication looks SO MUCH BETTER! It's still not fantastic, but my samples are clearly lining up. I saved both the nontransformed and log transformed plots. Next, I'm going to use the nontransformed NMDS ordinations to calculate distances between my technical replicates.

#### REFORMAT DATAFRAME FOR NMDS ####

#The goal is to have the row names of my new dataframe be Protein/Peptides/Transitions, with the column names as the sample number

#My first step is to change my dataframe from long to wide (i.e. cast it)
library(reshape2) #Instal package to pivot table
SRMDataNMDSPivoted <- dcast(SRMNormalizedDataNMDS, Protein.Name + Peptide.Sequence + Fragment.Ion ~ Sample.Number) #Cast table! Protein/Peptides/Transitions remain as columns with Sample Number as column headers. Normalized.Area used as value column by default.
head(SRMDataNMDSPivoted) #Confirm cast.
SRMDataNMDSPivoted$RowNames <- paste(SRMDataNMDSPivoted$Protein.Name, SRMDataNMDSPivoted$Peptide.Sequence, SRMDataNMDSPivoted$Fragment.Ion) #Merge Protein, Peptide and Transition information into one column
head(SRMDataNMDSPivoted) #Confirm column merge
SRMDataNMDSPivoted <- SRMDataNMDSPivoted[,-c(1:3)] #Remove unmerged columns
head(SRMDataNMDSPivoted) #Confirm column removal
#write.csv(SRMDataNMDSPivoted, file = "2017-09-11-SRM-Data-Normalized-NMDS-Pivoted.csv") #Wrote out as .csv to make future analyses easier.

#### NMDS PLOT ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste.
install.packages("vegan") #Install vegan package
library(vegan)

SRMDataNMDSPivotedCorrected <- SRMDataNMDSPivoted #Duplicate dataframe
SRMDataNMDSPivotedCorrected[is.na(SRMDataNMDSPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSPivotedCorrected) #Confirm there are no NAs

area.protID2 <- SRMDataNMDSPivotedCorrected[-93] #Save all area data as a new dataframe
rownames(area.protID2) <- SRMDataNMDSPivotedCorrected[,93] #Make sure last column of protein names is recognized as row names instead of values
head(area.protID2) #Confirm changes

area2.t <- t(area.protID2) #Transpose the file so that rows and columns are switched
head(area2.t) #Confirm transposition
area2.tra <- (area2.t+1) #Add 1 to all values before transforming
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.euclidean <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.euclidean) #Make Shepard plot
ordiplot(proc.nmds.euclidean) #Plot basic NMDS
vec.proc.nmds.euclidean <- envfit(proc.nmds.euclidean$points, area2.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
plot(vec.proc.nmds.euclidean, p.max=.01, col='blue') #Plot eigenvectors

proc.nmds.euclidean.log <- metaMDS(area2.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance
#stressplot(proc.nmds.euclidean.log) #Make Shepard plot
#ordiplot(proc.nmds.euclidean.log) #Plot basic NMDS
ordiplot(proc.nmds.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names

proc.nmds.euclidean.autotransform <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = TRUE) #Make MDS dissimilarity matrix using euclidean distance and autotransformation
#stressplot(proc.nmds.euclidean.autotransform) #Make Shepard plot
#ordiplot(proc.nmds.euclidean.autotransform) #Plot basic NMDS
ordiplot(proc.nmds.euclidean.autotransform, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names

#jpeg(filename = "2017-09-08-NMDS-TechnicalReplication-Normalized.jpeg", width = 1000, height = 1000)
#ordiplot(proc.nmds.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#### CALCULATE DISTANCES BETWEEN TECHNICAL REPLICATE ORDINATIONS ####

NMDSCoordinates <- proc.nmds.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinates) #Confirm dataframe creation
nSamples <- length(NMDSCoordinates)/2 #Calculate the number of samples
sampleDistances <- vector(length = nSamples) #Create an empty vector to store distance values
for(i in 1:nSamples) { #For rows in NMDSCoordinates
  sampleDistances[i] <- sqrt((NMDSCoordinates[i,1]-NMDSCoordinates[i,2])^2 + (NMDSCoordinates[i+1,1]-NMDSCoordinates[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistances[i]) #Print the distance value
}
sampleDistances #Confirm vector creation. This vector has all consecutive pairs, including those that are not paris of technical replicates. I need to retain just the odd numbered rows.
technicalReplicates <- rownames(NMDSCoordinates) #Save rownames as a new vector
technicalReplicates #Confirm vector creation
technicalReplicateDistances <- data.frame(Sample = technicalReplicates[seq(from = 1, to = nSamples, by = 2)], 
                                          Distance = sampleDistances[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistances) #Confirm dataframe creation
tail(technicalReplicateDistances) #Confirm dataframe creation

#### PLOT DISTANCES BETWEEN TECHNICAL REPLICATE ORDINATIONS ####

#jpeg(filename = "2017-09-08-NMDS-TechnicalReplication-Ordination-Distances.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistances$Sample, y = technicalReplicateDistances$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()