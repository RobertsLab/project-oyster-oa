#In this script, I'll use an NMDS plot to see how my technical replication fared after Steven filtered out transitions with a coefficient of variation greater than 15.

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
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-NMDS-TechnicalReplication-Normalized-after-CV15-Filtering.jpeg", width = 1000, height = 1000) #Save plot
ordiplot(proc.nmds.norm.euclidean, choices = c(1,2), type = "text", display = "sites", cex = .8) #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.euclidean, p.max=.01, col='blue') #Plot eigenvectors
#dev.off() #Turn off plotting mechanism

#proc.nmds.norm.euclidean.log <- metaMDS(area.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance
#stressplot(proc.nmds.norm.euclidean.log) #Make Shepard plot
#ordiplot(proc.nmds.norm.euclidean.log, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names

#### CALCULATE DISTANCES BETWEEN TECHNICAL REPLICATE ORDINATIONS ####

NMDSCoordinates <- proc.nmds.norm.euclidean$points #Save NMDS coordinates of each point in a new dataframe
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

#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-NMDS-TechnicalReplication-Ordination-Distances-after-CV15-Filtering.jpeg", width = 1000, height = 1000)
plot(x = technicalReplicateDistances$Sample, y = technicalReplicateDistances$Distance, type = "line", xlab = "Sample", ylab = "Distance between Ordinations")
#dev.off()

#### IDENTIFY SAMPLES WITH LARGE DISTANCES BETWEEN ORDINATIONS ####
#To identify samples with ordination distances that are outliers, I'm going to define an upper fence and remove samples above it.

histogram(technicalReplicateDistances$Distance) #Make a histogram of how many distance fall in bins. Looks like a conservative upper fence would be 0.2, which would allow me to remove samples with the highest ordination distances.
removeThese <- technicalReplicateDistances$Sample[technicalReplicateDistances$Distance >= 0.2] #Identify samples that need to be removed.
removeThese <- gsub(" -1", "", removeThese) #Remove " -1" from the end of each sample ID
removeThese #Confirm changes

#### REMOVE SAMPLES WITH LARGE DISTANCES BETWEEN ORDINATIONS ####

SRMModifiedAreasAdjusted <- SRMModifiedAreas[! SRMModifiedAreas$Sample %in% removeThese, ] #Duplicate original dataframe, but remove samples with large ordination distances saved in the vector removeThese

#### REMAKE NMDS FOR TECHNICAL REPLICATION ####

#Reformat data
SRMModifiedAreasAdjustedReplicate1 <- SRMModifiedAreasAdjusted[, 1:3] #Duplicate dataframe, keeping only Replicate 1 information
SRMModifiedAreasAdjustedReplicate1$Sample.Number <- paste(SRMModifiedAreasAdjustedReplicate1$Sample, "-1") #Create a new column with sample and replicate number
SRMModifiedAreasAdjustedReplicate1 <- SRMModifiedAreasAdjustedReplicate1[, -2] #Remove Sample column
colnames(SRMModifiedAreasAdjustedReplicate1) <- c("Protein.Name", "Area", "Sample") #Change column names
head(SRMModifiedAreasAdjustedReplicate1) #Confirm changes

SRMModifiedAreasAdjustedReplicate2 <- SRMModifiedAreasAdjusted[, c(1, 2, 4)] #Duplicate dataframe, keeping only Replicate 2 information
SRMModifiedAreasAdjustedReplicate2$Sample.Number <- paste(SRMModifiedAreasAdjustedReplicate2$Sample, "-2") #Create a new column with sample and replicate number
SRMModifiedAreasAdjustedReplicate2 <- SRMModifiedAreasAdjustedReplicate2[, -2] #Remove Sample column
colnames(SRMModifiedAreasAdjustedReplicate2) <- c("Protein.Name", "Area", "Sample") #Change column names
head(SRMModifiedAreasAdjustedReplicate2) #Confirm changes

SRMModifiedAreasAdjustedLong <- rbind(SRMModifiedAreasAdjustedReplicate1, SRMModifiedAreasAdjustedReplicate2) #Paste dataframes together
head(SRMModifiedAreasAdjustedLong) #Confirm changes
tail(SRMModifiedAreasAdjustedLong) #Confirm changes
transform(SRMModifiedAreasAdjustedLong, Area = as.numeric(Area)) #Make sure Area is recognized as a numeric variable
is.numeric(SRMModifiedAreasAdjustedLong$Area) #Confirm change

SRMModifiedAreasAdjustedPivoted <- dcast(SRMModifiedAreasAdjustedLong, Protein.Name ~ Sample, value.var = "Area") #Cast table! Protein.Name remains as a column with Sample Number as column headers. Area column defined as values.
head(SRMModifiedAreasAdjustedPivoted) #Confirm cast
rownames(SRMModifiedAreasAdjustedPivoted) <- SRMModifiedAreasAdjustedPivoted[, 1] #Save Protein.Name column as rownames
SRMModifiedAreasAdjustedPivoted <- SRMModifiedAreasAdjustedPivoted[, -1] #Remove Protein.Name column
head(SRMModifiedAreasAdjustedPivoted) #Confirm changes

SRMModifiedAreasAdjustedPivotedCorrected <- SRMModifiedAreasAdjustedPivoted #Duplicate dataframe
SRMModifiedAreasAdjustedPivotedCorrected[is.na(SRMModifiedAreasAdjustedPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMModifiedAreasAdjustedPivotedCorrected) #Confirm there are no NAs

#Make NMDS plot

area.prot2ID <- SRMModifiedAreasAdjustedPivotedCorrected #Save all area data as a new dataframe
head(area.prot2ID) #Confirm changes

area2.t <- t(area.prot2ID) #Transpose the file so that rows and columns are switched
head(area2.t) #Confirm transposition
area2.tra <- (area2.t+1) #Add 1 to all values before transforming
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.adj.euclidean <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.norm.adj.euclidean) #Make Shepard plot
#ordiplot(proc.nmds.norm.adj.euclidean) #Plot basic NMDS
#vec.proc.nmds.norm.adj.euclidean <- envfit(proc.nmds.norm.adj.euclidean$points, area.t, perm = 1000) #Calculate loadings
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-NMDS-TechnicalReplication-Normalized-after-CV15-and-Distance-Filtering.jpeg", width = 1000, height = 1000) #Save plot
ordiplot(proc.nmds.norm.adj.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#plot(vec.proc.nmds.norm.adj.euclidean, p.max=.01, col='blue') #Plot eigenvectors
#dev.off() #Turn off plotting mechanism

#### REMAKE NMDS FOR ANALYSES ####

#Average technical replicates
SRMAveragedAreas <- SRMModifiedAreasAdjusted #Dupliate dataframe
SRMAveragedAreas$Average.Area <- ((SRMAveragedAreas$Replicate1 + SRMAveragedAreas$Replicate2)/2) #Average peak areas and save as a new column
SRMAveragedAreas <- SRMAveragedAreas[, -c(3:4)] #Remove replicate area columns
head(SRMAveragedAreas) #Confirm changes

#Reformat data
SRMAveragedAreasPivoted <- dcast(SRMAveragedAreas, Protein.Name ~ Sample) #Cast table!
head(SRMAveragedAreasPivoted) #Confirm table cast!

SRMAveragedAreasPivotedCorrected <- SRMAveragedAreasPivoted #Duplicate dataframe
SRMAveragedAreasPivotedCorrected[is.na(SRMAveragedAreasPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMAveragedAreasPivotedCorrected) #Confirm there are no NAs

rownames(SRMAveragedAreasPivotedCorrected) <- SRMAveragedAreasPivotedCorrected$Protein.Name #Save protein name column as rownames
SRMAveragedAreasPivotedCorrected <- SRMAveragedAreasPivotedCorrected[, -1] #Remove Protein.Name column
head(SRMAveragedAreasPivotedCorrected) #Confirm changes
write.csv(SRMAveragedAreasPivotedCorrected, "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-Averaged-Pivoted-Corrected-SRM-Data-after-CV15-and-Distance-Filtering.csv") #Wrote out to make subsequent analyses easier.

#Make base NMDS plot
area.prot3ID <- SRMAveragedAreasPivotedCorrected #Save all area data as a new dataframe
head(area.prot3ID) #Confirm change

area3.t <- t(area.prot3ID) #Transpose the file so that rows and columns are switched
head(area3.t) #Confirm transposition
area3.tra <- (area3.t+1) #Add 1 to all values before transforming
area3.tra <- data.trans(area3.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.averaged.adjusted.euclidean <- metaMDS(area3.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.norm.averaged.adjusted.euclidean) #Make Shepard plot
#vec.proc.nmds.norm.averaged.adjusted.euclidean <- envfit(proc.nmds.norm.averaged.adjusted.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.averaged.adjusted.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.norm.averaged.adjusted.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#Assign colors and shapes
sampleIDs <- c("O01", "O04", "O06", "O08", "O10", "O100", "O101", "O102", "O103", "O106", "O118", "O121", "O122", "O124", "O128", "O131", "O137", "O14", "O140", "O145", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O49", "O51", "O52", "O56", "O60", "O64", "O66", "O71", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector
sampleIDs <- sampleIDs[! sampleIDs %in% removeThese] #Remove eliminated samples from ID vector

biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
biologicalReplicates <- biologicalReplicates[biologicalReplicates$sampleIDs %in% sampleIDs, ] #Match sample IDs in vector

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

#Refine NMDS
#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-NMDS-Norm-Analysis-Averaged-Adjusted-after-CV15-and-Distance-Filtering.jpeg", width = 1000, height = 1000)
fig.nmds <- ordiplot(proc.nmds.norm.averaged.adjusted.euclidean, choices = c(1,2), type = "none", display = "sites", xlab= "Axis 1", ylab= "Axis 2", cex = 0.5) #Save NMDS as a new object

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

#### REPEAT ANOSIM ####
dissimArea.t <- vegdist(area3.t, "euclidean") #Calculate dissimilarity matrix
ANOSIMReplicates <- biologicalReplicates #Duplicate dataframe
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
head(ANOSIMReplicates) #Confirm changes

ANOSIMReplicates$eelgrassCondition <- factor(ANOSIMReplicates$eelgrassCondition) #Make sure residual factors are no longer present
ANOSIMReplicates$site <- factor(ANOSIMReplicates$site) #Make sure residual factors are no longer present
str(ANOSIMReplicates) #Confirm new factors

siteNormAdjANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site presence
summary(siteNormAdjANOSIM)
plot(siteNormAdjANOSIM)
simper(proc.nmds.norm.averaged.adjusted.euclidean, ANOSIMReplicates$site)

eelgrassNormAdjANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassNormAdjANOSIM)
plot(eelgrassNormAdjANOSIM)
simper(proc.nmds.norm.averaged.adjusted.euclidean, ANOSIMReplicates$eelgrassCondition)