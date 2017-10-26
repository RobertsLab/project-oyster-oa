#In this script, I'll use an R-squared cutoff of 0.6 to eliminate transitions from my dataset. Then I will average my technical replicates and proceed with NMDS and ANOSIM analysis.\

#### SET WORKING DIRECTORY ####

setwd("../../..") #Set working directory as the main SRM directory, and not the subdirectory with the script.
getwd()

#### IMPORT DATA ####
SRMDataNMDSNormalizedPivoted <- read.csv("2017-09-11-SRM-Data-Normalized-NMDS-Pivoted.csv", header = TRUE) #Import pivoted and normalized data from first technical replication script
head(SRMDataNMDSNormalizedPivoted) #Confirm import
rownames(SRMDataNMDSNormalizedPivoted) <- SRMDataNMDSNormalizedPivoted$RowNames #Assign rownames
SRMDataNMDSNormalizedPivoted <- SRMDataNMDSNormalizedPivoted[,-1] #Remove extraneous column
SRMDataNMDSNormalizedPivoted <- SRMDataNMDSNormalizedPivoted[,-93] #Remove RowNames column
head(SRMDataNMDSNormalizedPivoted) #Confirm changes

#### ELIMINATE TRANSITIONS ####

#Remove transitions with R-squared values below 0.6 cutoff. Based on analyses in 2017-10-10-Correlating-Transitions-in-Technical-Replicates.R, the following rows (individual transitions) should be removed: 1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111

SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivoted #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1[-c(1, 2, 3, 18, 21, 22, 28, 43, 55, 76, 85, 86, 87, 88, 89, 90, 91, 92, 93, 103, 106, 109, 111),] #Remove rows that don't make the cutoff
rownames(SRMDataNMDSNormalizedPivotedCutoff1) #See if any peptide has less than 2 transitions, and if any protein has less than 2 peptides remaining

#I found CHOYP_PSA.1.1|m.27259 had one peptide with only one transition, so I removed that as well.
SRMDataNMDSNormalizedPivotedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1[-88,]
rownames(SRMDataNMDSNormalizedPivotedCutoff1) #Confirm changes

#### LOAD DEPENDENCIES ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste
install.packages("vegan") #Install vegan package
library(vegan)

#### REFORMAT DATA FOR NMDS ####
SRMDataNMDSNormalizedPivotedCorrectedCutoff1 <- SRMDataNMDSNormalizedPivotedCutoff1 #Duplicate dataframe
SRMDataNMDSNormalizedPivotedCorrectedCutoff1[is.na(SRMDataNMDSNormalizedPivotedCorrectedCutoff1)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Confirm there are no NAs
#write.csv(SRMDataNMDSNormalizedPivotedCorrectedCutoff1, "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-25-Cutoff-0.6/2017-10-25-SRM-Data-NMDS-Normalized-Pivoted-Corrected-Cutoff0.6-Filtered.csv") #Save file for future use.

area.t <- t(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

#### MADE NMDS FOR TECHNICAL REPLICATION ####

proc.nmds.norm.cutoff1.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance.
stressplot(proc.nmds.norm.cutoff1.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.cutoff1.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.cutoff1.euclidean <- envfit(proc.nmds.norm.cutoff1.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.cutoff1.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.cutoff1.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### AVERAGE TECHNICAL REPLICATES ####

head(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Dataset I'll use to average technical replicates
SRMDataNMDSNormalizedCutoff1Averaged <- data.frame(x = rep(x = 0, times = length(SRMDataNMDSNormalizedPivotedCorrectedCutoff1$O01.1)),
                                               y = rep(x = 0, times = length(SRMDataNMDSNormalizedPivotedCorrectedCutoff1$O01.1))) #Create an empty dataframe to store averaged values
row.names(SRMDataNMDSNormalizedCutoff1Averaged) <- row.names(SRMDataNMDSNormalizedPivotedCorrectedCutoff1) #Add row names
head(SRMDataNMDSNormalizedCutoff1Averaged) #Confirm changes

for(i in 1:(length(SRMDataNMDSNormalizedPivotedCorrectedCutoff1)-1)) { #Average normalized area values for consecutive columns
  SRMDataNMDSNormalizedCutoff1Averaged[,i] <- (SRMDataNMDSNormalizedPivotedCorrectedCutoff1[,i]+SRMDataNMDSNormalizedPivotedCorrectedCutoff1[,i+1])/2
}
head(SRMDataNMDSNormalizedCutoff1Averaged) #Confirm averaging

SRMDataNMDSNormalizedAveraged <- SRMDataNMDSNormalizedCutoff1Averaged[seq(from = 1, to = (length(SRMDataNMDSNormalizedPivotedCorrectedCutoff1)-1), by = 2)] #Remove even-numbered columns, since those consecutive columns are not technical replicates
head(SRMDataNMDSNormalizedAveraged) #Confirm column removal

sampleIDs <- c("O01", "O04", "O06", "O08", "O10", "O100", "O101", "O102", "O103", "O106", "O118", "O121", "O122", "O124", "O128", "O131", "O137", "O14", "O140", "O145", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O49", "O51", "O52", "O56", "O60", "O64", "O66", "O71", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector
colnames(SRMDataNMDSNormalizedAveraged) <- sampleIDs #Add column names
colnames(SRMDataNMDSNormalizedAveraged) #Confirm column naming

SRMDataNMDSNormalizedAveragedCorrected <- SRMDataNMDSNormalizedAveraged #Duplicate dataframe
SRMDataNMDSNormalizedAveragedCorrected[is.na(SRMDataNMDSNormalizedAveragedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSNormalizedAveragedCorrected) #Confirm there are no NAs

#### NMDS FOR SITE AND EELGRASS CLUSTERING ####

area.prot2ID <- SRMDataNMDSNormalizedAveragedCorrected #Save all area data as a new dataframe
head(area.prot2ID) #Confirm changes

area2.t <- t(area.prot2ID) #Transpose the file so that rows and columns are switched
head(area2.t) #Confirm transposition
area2.tra <- (area2.t+1) #Add 1 to all values before transforming
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.averaged.euclidean <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.norm.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.norm.averaged.euclidean <- envfit(proc.nmds.norm.averaged.euclidean$points, area4.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.norm.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### ASSIGN COLORS AND SHAPES ####
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

#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-25-Cutoff-0.6/2017-10-25-NMDS-Norm-Analysis-Averaged-Cutoff0.6.jpeg", width = 1000, height = 1000)
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

dissimArea.t <- vegdist(area2.t, "euclidean") #Calculate dissimilarity matrix
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
simper(proc.nmds.norm.averaged.euclidean, ANOSIMReplicates$site)

eelgrassNormANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassNormANOSIM)
plot(eelgrassNormANOSIM)
simper(proc.nmds.norm.averaged.euclidean, ANOSIMReplicates$eelgrassCondition)

#### CALCULATE DISTANCES BETWEEN ORDINATIONS ####
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

#### IDENTIFY SAMPLES WITH LARGE DISTANCES BETWEEN ORDINATIONS ####
#To identify samples with ordination distances that are outliers, I'm going to define an upper fence and remove samples above it.

histogram(technicalReplicateDistancesNormalizedCuttof1$Distance) #Make a histogram of how many distance fall in bins. Looks like a conservative upper fence would be 0.15, which would allow me to remove samples with the highest ordination distances.
removeThese1 <- technicalReplicateDistancesNormalizedCuttof1$Sample[technicalReplicateDistancesNormalizedCuttof1$Distance >= 0.15] #Identify samples that need to be removed.
removeThese2 <- gsub("\\.1\\>", "\\.2\\", removeThese1) #Make list of replicates that should also be removed
removeThese2 #Confirm duplication
removeThese1 <- gsub("\\.2\\>", "\\.1\\", removeThese2) #Take residual levels out of list
removeThese1 #Confirm changes
removeThese <- c(removeThese1, removeThese2)
removeThese #Confirm changes

#Could continue to remove samples with large ordination distances, but I stopped here because Steven is not confident with the value of an R-squared cutoff. I will focus my energy on CV filtering instead.