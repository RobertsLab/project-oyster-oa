#In this script, I'll use three different adjusted R-squared cutoffs to screen transition data and reassess my technical replication. I'll try these cutoffs both with and without normalizing my data.

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

sampleDistancesNonNormalizedCutoff1 #Confirm vector creation. This vector has all consecutive pairs, including those that are not pairs of technical replicates. I need to retain just the odd numbered rows.
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

#### CUTOFF = 0.7, NONNORMALIZED DATA ####

#Remove transitions with R-squared values below 0.7 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 4, 5, 6, 18, 19, 21, 22, 26, 28, 30, 43, 47, 55, 67, 68, 69, 70, 76, 77, 79, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 103, 106, 107, 108, 109, 110, 111

SRMDataNMDSNonNormalizedPivotedCutoff2 <- SRMDataNMDSNonNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNonNormalizedPivotedCutoff2 <- SRMDataNMDSNonNormalizedPivotedCutoff2[-c(1, 2, 3, 4, 5, 6, 18, 19, 21, 22, 26, 28, 30, 43, 47, 55, 67, 68, 69, 70, 76, 77, 79, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 103, 106, 107, 108, 109, 110, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNonNormalizedPivotedCutoff2) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining
SRMDataNMDSNonNormalizedPivotedCutoff2 <- SRMDataNMDSNonNormalizedPivotedCutoff2[-c(1:3, 12, 15:17, 48:50, 56, 57:65),] #Remove peptides and proteins that don't fit the above criteria. I cut more than I expected.
rownames(SRMDataNMDSNonNormalizedPivotedCutoff2) #Confirm changes

#Load the source file for the biostats package
#source("biostats.R") #Either load the source R script or copy paste
#install.packages("vegan") #Install vegan package
#library(vegan)

#Format data for NMDS
SRMDataNMDSNonNormalizedPivotedCorrectedCutoff2 <- SRMDataNMDSNonNormalizedPivotedCutoff2 #Duplicate dataframe
SRMDataNMDSNonNormalizedPivotedCorrectedCutoff2[is.na(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff2)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff2) #Confirm there are no NAs

area.t3 <- t(SRMDataNMDSNonNormalizedPivotedCorrectedCutoff2) #Transpose the file so that rows and columns are switched
head(area.t3) #Confirm transposition
area.tra3 <- (area.t3+1) #Add 1 to all values before transforming
area.tra3 <- data.trans(area.tra3, method = 'log', plot = FALSE) #log(x+1) transformation

#Rerun NMDS
proc.nmds.nonnorm.cutoff2.euclidean <- metaMDS(area.t3, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.nonnorm.cutoff2.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.cutoff2.euclidean) #Plot basic NMDS
#vec.proc.nmds.nonnorm.cutoff2.euclidean <- envfit(proc.nmds.nonnorm.cutoff2.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.nonnorm.cutoff2.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.nonnorm.cutoff2.euclidean, p.max=.01, col='blue') #Plot eigenvectors

jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-NonNormalized-Cutoff2.jpeg", width = 1000, height = 1000)
ordiplot(proc.nmds.nonnorm.cutoff2.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
dev.off()

#Calculate distance between ordinations
NMDSCoordinatesNonNormalizedCutoff2 <- proc.nmds.nonnorm.cutoff2.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinatesNonNormalizedCutoff2) #Confirm dataframe creation
nSamples <- length(NMDSCoordinatesNonNormalizedCutoff2)/2 #Calculate the number of samples
sampleDistancesNonNormalizedCutoff2 <- vector(length = nSamples) #Create an empty vector to store distance values

for(i in 1:nSamples) { #For rows in NMDSCoordinatesNonNormalizedCutoff1
  sampleDistancesNonNormalizedCutoff2[i] <- sqrt((NMDSCoordinatesNonNormalizedCutoff2[i,1]-NMDSCoordinatesNonNormalizedCutoff2[i,2])^2 + (NMDSCoordinatesNonNormalizedCutoff2[i+1,1]-NMDSCoordinatesNonNormalizedCutoff2[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistancesNonNormalizedCutoff2[i]) #Print the distance value
}

sampleDistancesNonNormalizedCutoff2 #Confirm vector creation. This vector has all consecutive pairs, including those that are not pairs of technical replicates. I need to retain just the odd numbered rows.
technicalReplicatesNonNormalizedCuttof2 <- rownames(NMDSCoordinatesNonNormalizedCutoff2) #Save rownames as a new vector
technicalReplicatesNonNormalizedCuttof2 #Confirm vector creation
technicalReplicateDistancesNonNormalizedCuttof2 <- data.frame(Sample = technicalReplicatesNonNormalizedCuttof2[seq(from = 1, to = nSamples, by = 2)], 
                                                              Distance = sampleDistancesNonNormalizedCutoff2[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistancesNonNormalizedCuttof2) #Confirm dataframe creation
tail(technicalReplicateDistancesNonNormalizedCuttof2) #Confirm dataframe creation

#Plot distance between technical replicate ordinations and save
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Ordination-Distances-NonNormalized-Cutoff2.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistancesNonNormalizedCuttof2$Sample, y = technicalReplicateDistancesNonNormalizedCuttof2$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()

#### CUTOFF = 0.7, NORMALIZED DATA ####

#Remove transitions with R-squared values below 0.7 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 4, 5, 6, 18, 19, 21, 22, 26, 28, 30, 43, 47, 55, 67, 68, 69, 70, 76, 77, 79, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 103, 106, 107, 108, 109, 110, 111

SRMDataNMDSNormalizedPivotedCutoff2 <- SRMDataNMDSNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCutoff2 <- SRMDataNMDSNormalizedPivotedCutoff2[-c(1, 2, 3, 4, 5, 6, 18, 19, 21, 22, 26, 28, 30, 43, 47, 55, 67, 68, 69, 70, 76, 77, 79, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 103, 106, 107, 108, 109, 110, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNormalizedPivotedCutoff2) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining
SRMDataNMDSNormalizedPivotedCutoff2 <- SRMDataNMDSNormalizedPivotedCutoff2[-c(1:3, 12, 15:17, 48:50, 56, 57:65),] #Remove peptides and proteins that don't fit the above criteria. I cut more than I expected.
rownames(SRMDataNMDSNormalizedPivotedCutoff2) #Confirm changes

#Load the source file for the biostats package
#source("biostats.R") #Either load the source R script or copy paste
#install.packages("vegan") #Install vegan package
#library(vegan)

#Format data for NMDS
SRMDataNMDSNormalizedPivotedCorrectedCutoff2 <- SRMDataNMDSNormalizedPivotedCutoff2 #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCorrectedCutoff2[is.na(SRMDataNMDSNormalizedPivotedCorrectedCutoff2)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNormalizedPivotedCorrectedCutoff2) #Confirm there are no NAs

area.t4 <- t(SRMDataNMDSNormalizedPivotedCorrectedCutoff2) #Transpose the file so that rows and columns are switched
head(area.t4) #Confirm transposition
area.tra4 <- (area.t4+1) #Add 1 to all values before transforming
area.tra4 <- data.trans(area.tra4, method = 'log', plot = FALSE) #log(x+1) transformation

#Rerun NMDS
proc.nmds.norm.cutoff2.euclidean <- metaMDS(area.t4, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.norm.cutoff2.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff2.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.cutoff2.euclidean <- envfit(proc.nmds.norm.cutoff2.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.cutoff2.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.cutoff1.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#Save plot
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Normalized-Cutoff2.jpeg", width = 1000, height = 1000)
#ordiplot(proc.nmds.norm.cutoff2.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#Calculate distance between ordinations
NMDSCoordinatesNormalizedCutoff2 <- proc.nmds.norm.cutoff2.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinatesNormalizedCutoff2) #Confirm dataframe creation
nSamples <- length(NMDSCoordinatesNormalizedCutoff2)/2 #Calculate the number of samples
sampleDistancesNormalizedCutoff2 <- vector(length = nSamples) #Create an empty vector to store distance values

for(i in 1:nSamples) { #For rows in NMDSCoordinatesNormalizedCutoff2
  sampleDistancesNormalizedCutoff2[i] <- sqrt((NMDSCoordinatesNormalizedCutoff2[i,1]-NMDSCoordinatesNormalizedCutoff2[i,2])^2 + (NMDSCoordinatesNormalizedCutoff2[i+1,1]-NMDSCoordinatesNormalizedCutoff2[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistancesNormalizedCutoff2[i]) #Print the distance value
}

sampleDistancesNormalizedCutoff2 #Confirm vector creation. This vector has all consecutive pairs, including those that are not paris of technical replicates. I need to retain just the odd numbered rows.
technicalReplicatesNormalizedCuttof2 <- rownames(NMDSCoordinatesNormalizedCutoff2) #Save rownames as a new vector
technicalReplicatesNormalizedCuttof2 #Confirm vector creation
technicalReplicateDistancesNormalizedCuttof2 <- data.frame(Sample = technicalReplicatesNormalizedCuttof2[seq(from = 1, to = nSamples, by = 2)], 
                                                           Distance = sampleDistancesNormalizedCutoff2[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistancesNormalizedCuttof2) #Confirm dataframe creation
tail(technicalReplicateDistancesNormalizedCuttof2) #Confirm dataframe creation

#Plot distance between technical replicate ordinations and save
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Ordination-Distances-Normalized-Cutoff2.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistancesNormalizedCuttof2$Sample, y = technicalReplicateDistancesNormalizedCuttof2$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()

#### CUTOFF = 0.8, NONNORMALIZED DATA ####
#Since my past two attempts showed that normalized data is cleaner than nonnormalized data, I will not be writing code for this section.

#### CUTOFF = 0.8, NORMALIZED DATA ####

#Remove transitions with R-squared values below 0.8 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 43, 44, 45, 46, 47, 48, 50, 52, 55, 56, 57, 64, 67, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 103, 104, 106, 107, 108, 109, 110, 111

SRMDataNMDSNormalizedPivotedCutoff3 <- SRMDataNMDSNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCutoff3 <- SRMDataNMDSNormalizedPivotedCutoff3[-c(1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 43, 44, 45, 46, 47, 48, 50, 52, 55, 56, 57, 64, 67, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 103, 104, 106, 107, 108, 109, 110, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNormalizedPivotedCutoff3) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining
SRMDataNMDSNormalizedPivotedCutoff3 <- SRMDataNMDSNormalizedPivotedCutoff3[-c(1:3, 7, 16:17, 24:31),] #Remove peptides and proteins that don't fit the above criteria. There are maybe two proteins left
rownames(SRMDataNMDSNormalizedPivotedCutoff3) #Confirm changes

#Format data for NMDS
SRMDataNMDSNormalizedPivotedCorrectedCutoff3 <- SRMDataNMDSNormalizedPivotedCutoff3 #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCorrectedCutoff3[is.na(SRMDataNMDSNormalizedPivotedCorrectedCutoff3)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNormalizedPivotedCorrectedCutoff3) #Confirm there are no NAs

area.t5 <- t(SRMDataNMDSNormalizedPivotedCorrectedCutoff3) #Transpose the file so that rows and columns are switched
head(area.t5) #Confirm transposition
area.tra5 <- (area.t5+1) #Add 1 to all values before transforming
area.tra5 <- data.trans(area.tra5, method = 'log', plot = FALSE) #log(x+1) transformation

#Rerun NMDS
proc.nmds.norm.cutoff3.euclidean <- metaMDS(area.t5, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.norm.cutoff3.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff3.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.cutoff3.euclidean <- envfit(proc.nmds.norm.cutoff3.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.cutoff3.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.cutoff3.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#Save plot
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Normalized-Cutoff3.jpeg", width = 1000, height = 1000)
#ordiplot(proc.nmds.norm.cutoff3.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#Calculate distance between ordinations
NMDSCoordinatesNormalizedCutoff3 <- proc.nmds.norm.cutoff3.euclidean$points #Save NMDS coordinates of each point in a new dataframe
head(NMDSCoordinatesNormalizedCutoff3) #Confirm dataframe creation
nSamples <- length(NMDSCoordinatesNormalizedCutoff3)/2 #Calculate the number of samples
sampleDistancesNormalizedCutoff3 <- vector(length = nSamples) #Create an empty vector to store distance values

for(i in 1:nSamples) { #For rows in NMDSCoordinatesNormalizedCutoff2
  sampleDistancesNormalizedCutoff3[i] <- sqrt((NMDSCoordinatesNormalizedCutoff3[i,1]-NMDSCoordinatesNormalizedCutoff3[i,2])^2 + (NMDSCoordinatesNormalizedCutoff3[i+1,1]-NMDSCoordinatesNormalizedCutoff3[i+1,2])^2) #Calculate distance between ordinations
  print(sampleDistancesNormalizedCutoff3[i]) #Print the distance value
}

sampleDistancesNormalizedCutoff3 #Confirm vector creation. This vector has all consecutive pairs, including those that are not paris of technical replicates. I need to retain just the odd numbered rows.
technicalReplicatesNormalizedCuttof3 <- rownames(NMDSCoordinatesNormalizedCutoff3) #Save rownames as a new vector
technicalReplicatesNormalizedCuttof3 #Confirm vector creation
technicalReplicateDistancesNormalizedCuttof3 <- data.frame(Sample = technicalReplicatesNormalizedCuttof3[seq(from = 1, to = nSamples, by = 2)], 
                                                           Distance = sampleDistancesNormalizedCutoff3[seq(from = 1, to = nSamples, by = 2)]) #Create a new dataframe with just odd numbered row distances (technical replicate pairs)
head(technicalReplicateDistancesNormalizedCuttof3) #Confirm dataframe creation
tail(technicalReplicateDistancesNormalizedCuttof3) #Confirm dataframe creation

#Plot distance between technical replicate ordinations and save
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-13-NMDS-TechnicalReplication-Ordination-Distances-Normalized-Cutoff3.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistancesNormalizedCuttof3$Sample, y = technicalReplicateDistancesNormalizedCuttof3$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()