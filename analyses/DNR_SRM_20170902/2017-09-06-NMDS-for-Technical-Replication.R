#In this script, I'll use an NMDS plot to see if my technical replicates are similar.

#### IMPORT DATA ####

SRMAreas <- read.csv("2017-09-11-Gigas-SRM-ReplicatesOnly-PostDilutionCurve-NoPivot-Report.csv", na.strings = "#N/A") #Specify Skyline's special way of designating N/A values
head(SRMAreas) #Confirm import
tail(SRMAreas) #Confirm import

#### CREATE A MASTER DATAFRAME ####

#I want to merge my Skyline data with sample names, sites, and eelgrass condition to create a master dataframe will all possible information

sequenceFile <- read.csv("2017-07-28-SRM-Samples-Sequence-File.csv", na.strings = "N/A") # Import sequence file
head(sequenceFile) #Confirm import
sequenceFile <- sequenceFile[,c(2,3,8)] #Keep the Replicate.Name, Comment and TIC columns
names(sequenceFile) <- c("Replicate.Name", "Sample.Number", "TIC")
head(sequenceFile) #Confirm change
masterSRMData <- merge(x = SRMAreas, y = sequenceFile, by = "Replicate.Name") #Merge the sample names and replicate names to use for analysis.Note: rows lost because 16's technical replicate is no longer included in the analysis.
head(masterSRMData) #Confirm merge
tail(masterSRMData) #Confirm merge

biologicalReplicates <- read.csv("2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
tail(biologicalReplicates) #Confirm import
masterSRMDataBiologicalReplicates <- merge(x = masterSRMData, y = biologicalReplicates, by = "Sample.Number") #Add biological replicate information to master list.
head(masterSRMDataBiologicalReplicates)
masterSRMDataBiologicalReplicates <- masterSRMDataBiologicalReplicates[,-8] #Remove TIC Area column since it is empty
head(masterSRMDataBiologicalReplicates) #Confirm change
#write.csv(x = masterSRMDataBiologicalReplicates, file = "2017-09-07-Master-SRM-Data-BiologicalReplicates-NoBlanks-NoPivot.csv") #Write out master dataframe

#### SUBSET DATA FOR NMDS PLOT ####

#For the NMDS, I want only the protein/peptide/transition information and peak area

SRMDataNMDS <- masterSRMDataBiologicalReplicates #Duplicate master list into a new dataframe
head(SRMDataNMDS) #Confirm copy
tail(SRMDataNMDS) #Confirm copy
SRMDataNMDS <- SRMDataNMDS[,-c(2, 5, 7, 10, 11)] #Remove extraneous columns: Replicate.Name, Transition, Peptide.Retention.Time, Site, Eelgrass
head(SRMDataNMDS) #Confirm column removal
SRMDataNMDS <- SRMDataNMDS[! SRMDataNMDS$Protein.Name %in% "PRTC peptides", ] #Remove PRTC peptide data
head(SRMDataNMDS) #Confirm removal
transform(SRMDataNMDS, Area = as.numeric(Area)) #Make sure Area is recognized as a numeric variable
is.numeric(SRMDataNMDS$Area) #Confirm change
transform(SRMDataNMDS, TIC = as.numeric(TIC)) #Make sure TIC is recognized as a numeric variable
is.numeric(SRMDataNMDS$TIC) #Confirm change

#### MAKE NMDS WITHOUT NORMALIZING ####

#The goal is to have the row names of my new dataframe be Protein/Peptides/Transitions, with the column names as the sample number

SRMDataNonNormalizedNMDS <- SRMDataNMDS #Create a duplicate dataframe
SRMDataNonNormalizedNMDS <- SRMDataNMDS[,-6] #Remove TIC column
head(SRMDataNonNormalizedNMDS) #Confirm creation

#My first step is to change my dataframe from long to wide (i.e. cast it)
library(reshape2) #Instal package to pivot table
SRMDataNMDSNonNormalizedPivoted <- dcast(SRMDataNonNormalizedNMDS, Protein.Name + Peptide.Sequence + Fragment.Ion ~ Sample.Number) #Cast table! Protein/Peptides/Transitions remain as columns with Sample Number as column headers. Normalized.Area used as value column by default.
head(SRMDataNMDSNonNormalizedPivoted) #Confirm cast.
SRMDataNMDSNonNormalizedPivoted$RowNames <- paste(SRMDataNMDSNonNormalizedPivoted$Protein.Name, SRMDataNMDSNonNormalizedPivoted$Peptide.Sequence, SRMDataNMDSNonNormalizedPivoted$Fragment.Ion) #Merge Protein, Peptide and Transition information into one column
head(SRMDataNMDSNonNormalizedPivoted) #Confirm column merge
SRMDataNMDSNonNormalizedPivoted <- SRMDataNMDSNonNormalizedPivoted[,-c(1:3)] #Remove unmerged columns
head(SRMDataNMDSNonNormalizedPivoted) #Confirm column removal
write.csv(SRMDataNMDSNonNormalizedPivoted, file = "2017-09-07-SRM-Data-NMDS-Pivoted.csv") #Wrote out as .csv to make future analyses easier.

#Now I can make an NMDS plot

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste
install.packages("vegan") #Install vegan package
library(vegan)

SRMDataNMDSNonNormalizedPivotedCorrected <- SRMDataNMDSNonNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNonNormalizedPivotedCorrected[is.na(SRMDataNMDSNonNormalizedPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNonNormalizedPivotedCorrected) #Confirm there are no NAs

area.protID <- SRMDataNMDSNonNormalizedPivotedCorrected[-99] #Save all area data as a new dataframe
rownames(area.protID) <- SRMDataNMDSNonNormalizedPivotedCorrected[,99] #Make sure last column of protein names is recognized as row names instead of values
head(area.protID) #Confirm changes

area.t <- t(area.protID) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.nonnorm.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.nonnorm.euclidean) #Make Shepard plot
ordiplot(proc.nmds.nonnorm.euclidean) #Plot basic NMDS
vec.proc.nmds.nonnorm.euclidean <- envfit(proc.nmds.nonnorm.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.nonnorm.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
plot(vec.proc.nmds.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#proc.nmds.nonnorm.euclidean.log <- metaMDS(area.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance
#stressplot(proc.nmds.nonnorm.euclidean.log) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.euclidean.log) #Plot basic NMDS
#ordiplot(proc.nmds.nonnorm.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names

#proc.nmds.nonnorm.euclidean.autotransform <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = TRUE) #Make MDS dissimilarity matrix using euclidean distance and autotransformation
#stressplot(proc.nmds.nonnorm.euclidean.autotransform) #Make Shepard plot
#ordiplot(proc.nmds.nonnorm.euclidean.autotransform) #Plot basic NMDS
#ordiplot(proc.nmds.nonnorm.euclidean.autotransform, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names

#jpeg(filename = "2017-09-11-NMDS-TechnicalReplication-NonNormalized.jpeg", width = 1000, height = 1000)
#ordiplot(proc.nmds.nonnorm.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#Only the euclidean non-transformed version gives me an NMDS plot, but there looks like there's variation that should be accounted for by normalizing.


#### NORMALIZE BY TIC VALUES ####

SRMNormalizedDataNMDS <- SRMDataNMDS #Duplicate dataframe
SRMNormalizedDataNMDS$Normalized.Area <- SRMNormalizedDataNMDS$Area/SRMDataNMDS$TIC #Divide areas by corresponding TIC values
head(SRMNormalizedDataNMDS) #Confirm division
SRMNormalizedDataNMDS <- SRMNormalizedDataNMDS[,-c(5,6)] #Remove nonnormalized area and TIC columns
head(SRMNormalizedDataNMDS) #Confirm column removal

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

area.protID2 <- SRMDataNMDSPivotedCorrected[-99] #Save all area data as a new dataframe
rownames(area.protID2) <- SRMDataNMDSPivotedCorrected[,99] #Make sure last column of protein names is recognized as row names instead of values
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