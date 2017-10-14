#In this script, I'll use two different adjusted R-squared cutoffs to screen transition data and reassess my technical replication. I'll try these cutoffs both with and without normalizing my data.

#### CUTOFF = 0.6, NONNORMALIZED DATA ####

#Import data
SRMDataNMDSNonNormalizedPivoted <- read.csv("2017-09-07-SRM-Data-NMDS-Pivoted.csv", header = TRUE) #Import pivoted data from first technical replication script
head(SRMDataNMDSNonNormalizedPivoted) #Confirm import
rownames(SRMDataNMDSNonNormalizedPivoted) <- SRMDataNMDSNonNormalizedPivoted$RowNames #Assign rownames
SRMDataNMDSNonNormalizedPivoted <- SRMDataNMDSNonNormalizedPivoted[,-1] #Remove extraneous column
SRMDataNMDSNonNormalizedPivoted <- SRMDataNMDSNonNormalizedPivoted[,-93] #Remove RowNames column
head(SRMDataNMDSNonNormalizedPivoted) #Confirm changes

#Remove transitions with R-squared values below 0.6 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111

SRMDataNMDSNonNormalizedPivotedCutoff1 <- SRMDataNMDSNonNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNonNormalizedPivotedCutoff1 <- SRMDataNMDSNonNormalizedPivotedCutoff1[-c(1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNonNormalizedPivotedCutoff1) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining

#I found CHOYP_PSA.1.1|m.27259 had one peptide with only one transition, so I removed that as well.
SRMDataNMDSNonNormalizedPivotedCutoff1 <- SRMDataNMDSNonNormalizedPivotedCutoff1[-88,]
rownames(SRMDataNMDSNonNormalizedPivotedCutoff1) #Confirm changes

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste
install.packages("vegan") #Install vegan package
library(vegan)

#Format data for NMDS
SRMDataNMDSNonNormalizedPivotedCorrectedCutoff1 <- SRMDataNMDSNonNormalizedPivotedCutoff1 #Duplicate dataframe
SRMDataNMDSNonNormalizedPivotedCorrectedCutoff1[is.na(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff1)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff1) #Confirm there are no NAs

area.t <- t(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff1) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

#Rerun NMDS
proc.nmds.nonnorm.cutoff1.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.nonnorm.cutoff1.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.cutoff1.euclidean) #Plot basic NMDS
#vec.proc.nmds.nonnorm.cutoff1.euclidean <- envfit(proc.nmds.nonnorm.cutoff1.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.nonnorm.cutoff1.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.nonnorm.cutoff1.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#proc.nmds.nonnorm.cutoff1.euclidean.log <- metaMDS(area.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance
#stressplot(proc.nmds.nonnorm.cutoff1.euclidean.log) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.cutoff1.euclidean.log) #Plot basic NMDS
#ordiplot(proc.nmds.nonnorm.cutoff1.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names. This is super messy!

#proc.nmds.nonnorm.cutoff1.euclidean.autotransform <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = TRUE) #Make MDS dissimilarity matrix using euclidean distance and autotransformation. Stress is (nearly) zero - you may have insufficient data
#stressplot(proc.nmds.nonnorm.cutoff1.euclidean.autotransform) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.cutoff1.euclidean.autotransform) #Plot basic NMDS
#ordiplot(proc.nmds.nonnorm.cutoff1.euclidean.autotransform, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names. Also not good.

#I'm going to save the NMDS plot from my first run (euclidean distances, untransformed data)
jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-NonNormalized-Cutoff1.jpeg", width = 1000, height = 1000)
ordiplot(proc.nmds.nonnorm.cutoff1.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
dev.off()

#Calculate distance between ordinations
NMDSCoordinatesNonNormalizedCutoff1 <- proc.nmds.nonnorm.cutoff1.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinatesNonNormalizedCutoff1) #Confirm dataframe creation
nSamples <- length(NMDSCoordinatesNonNormalizedCutoff1)/2 #Calculate the number of samples
sampleDistancesNonNormalizedCutoff1 <- vector(length = nSamples) #Create an empty vector to store distance values

for(i in 1:nSamples) { #For rows in NMDSCoordinatesNonNormalizedCutoff1
  sampleDistancesNonNormalizedCutoff1[i] <- sqrt((NMDSCoordinatesNonNormalizedCutoff1[i,1]-NMDSCoordinatesNonNormalizedCutoff1[i,2])^2 + (NMDSCoordinatesNonNormalizedCutoff1[i+1,1]-NMDSCoordinatesNonNormalizedCutoff1[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistancesNonNormalizedCutoff1[i]) #Print the distance value
}
sampleDistancesNonNormalizedCutoff1 #Confirm vector creation. This vector has all consecutive pairs, including those that are not paris of technical replicates. I need to retain just the odd numbered rows.
technicalReplicatesNonNormalizedCuttof1 <- rownames(NMDSCoordinatesNonNormalizedCutoff1) #Save rownames as a new vector
technicalReplicatesNonNormalizedCuttof1 #Confirm vector creation
technicalReplicateDistancesNonNormalizedCuttof1 <- data.frame(Sample = technicalReplicatesNonNormalizedCuttof1[seq(from = 1, to = nSamples, by = 2)], 
                                                             Distance = sampleDistancesNonNormalizedCutoff1[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistancesNonNormalizedCuttof1) #Confirm dataframe creation
tail(technicalReplicateDistancesNonNormalizedCuttof1) #Confirm dataframe creation

#Plot distance between technical replicate ordinations and save
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Ordination-Distances-NonNormalized-Cutoff1.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistancesNonNormalizedCuttof1$Sample, y = technicalReplicateDistancesNonNormalizedCuttof1$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()

#### CUTOFF = 0.6, NORMALIZED DATA ####

#Import data
SRMDataNMDSNormalizedPivoted <- read.csv("2017-09-11-SRM-Data-Normalized-NMDS-Pivoted.csv", header = TRUE) #Import pivoted and normalized data from first technical replication script
head(SRMDataNMDSNormalizedPivoted) #Confirm import
rownames(SRMDataNMDSNormalizedPivoted) <- SRMDataNMDSNormalizedPivoted$RowNames #Assign rownames
SRMDataNMDSNormalizedPivoted <- SRMDataNMDSNormalizedPivoted[,-1] #Remove extraneous column
SRMDataNMDSNormalizedPivoted <- SRMDataNMDSNormalizedPivoted[,-93] #Remove RowNames column
head(SRMDataNMDSNormalizedPivoted) #Confirm changes

#Remove transitions with R-squared values below 0.6 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111

SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1[-c(1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNormalizedPivotedCutoff1) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining

#I found CHOYP_PSA.1.1|m.27259 had one peptide with only one transition, so I removed that as well.
SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1[-88,]
rownames(SRMDataNMDSNormalizedPivotedCutoff1) #Confirm changes

#Load the source file for the biostats package
#source("biostats.R") #Either load the source R script or copy paste
#install.packages("vegan") #Install vegan package
#library(vegan)

#Format data for NMDS
SRMDataNMDSNormalizedPivotedCorrectedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1 #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCorrectedCutoff1[is.na(SRMDataNMDSNormalizedPivotedCorrectedCutoff1)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Confirm there are no NAs

area.t2 <- t(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Transpose the file so that rows and columns are switched
head(area.t2) #Confirm transposition
area.tra2 <- (area.t2+1) #Add 1 to all values before transforming
area.tra2 <- data.trans(area.tra2, method = 'log', plot = FALSE) #log(x+1) transformation

#Rerun NMDS
proc.nmds.norm.cutoff1.euclidean <- metaMDS(area.t2, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.norm.cutoff1.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff1.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.cutoff1.euclidean <- envfit(proc.nmds.norm.cutoff1.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.cutoff1.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.cutoff1.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#proc.nmds.norm.cutoff1.euclidean.log <- metaMDS(area.tra2, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance and log transformed data
#stressplot(proc.nmds.norm.cutoff1.euclidean.log) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff1.euclidean.log) #Plot basic NMDS
#ordiplot(proc.nmds.norm.cutoff1.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names. No difference between transforming and not transforming.

proc.nmds.norm.cutoff1.euclidean.autotransform <- metaMDS(area.t2, distance = 'euclidean', k = 2, trymax = 10000, autotransform = TRUE) #Make MDS dissimilarity matrix using euclidean distance and autotransformation.
#stressplot(proc.nmds.norm.cutoff1.euclidean.autotransform) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff1.euclidean.autotransform) #Plot basic NMDS
ordiplot(proc.nmds.norm.cutoff1.euclidean.autotransform, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names. Also no real difference.

#I'm going to save the NMDS plot from my first run (euclidean distances, untransformed data)
jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Normalized-Cutoff1.jpeg", width = 1000, height = 1000)
ordiplot(proc.nmds.norm.cutoff1.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
dev.off()

#Calculate distance between ordinations
NMDSCoordinatesNormalizedCutoff1 <- proc.nmds.norm.cutoff1.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinatesNormalizedCutoff1) #Confirm dataframe creation
nSamples <- length(NMDSCoordinatesNormalizedCutoff1)/2 #Calculate the number of samples
sampleDistancesNormalizedCutoff1 <- vector(length = nSamples) #Create an empty vector to store distance values

for(i in 1:nSamples) { #For rows in NMDSCoordinatesNormalizedCutoff1
  sampleDistancesNormalizedCutoff1[i] <- sqrt((NMDSCoordinatesNormalizedCutoff1[i,1]-NMDSCoordinatesNormalizedCutoff1[i,2])^2 + (NMDSCoordinatesNormalizedCutoff1[i+1,1]-NMDSCoordinatesNormalizedCutoff1[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistancesNormalizedCutoff1[i]) #Print the distance value
}

sampleDistancesNormalizedCutoff1 #Confirm vector creation. This vector has all consecutive pairs, including those that are not paris of technical replicates. I need to retain just the odd numbered rows.
technicalReplicatesNormalizedCuttof1 <- rownames(NMDSCoordinatesNormalizedCutoff1) #Save rownames as a new vector
technicalReplicatesNormalizedCuttof1 #Confirm vector creation
technicalReplicateDistancesNormalizedCuttof1 <- data.frame(Sample = technicalReplicatesNormalizedCuttof1[seq(from = 1, to = nSamples, by = 2)], 
                                                              Distance = sampleDistancesNormalizedCutoff1[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistancesNormalizedCuttof1) #Confirm dataframe creation
tail(technicalReplicateDistancesNormalizedCuttof1) #Confirm dataframe creation

#Plot distance between technical replicate ordinations and save
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Ordination-Distances-Normalized-Cutoff1.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistancesNormalizedCuttof1$Sample, y = technicalReplicateDistancesNormalizedCuttof1$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()

#### CUTOFF = 0.8, NONNORMALIZED DATA ####

#### CUTOFF = 0.8, NORMALIZED DATA ####