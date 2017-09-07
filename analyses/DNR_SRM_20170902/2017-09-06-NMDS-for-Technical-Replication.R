#In this script, I'll use an NMDS plot to see if my technical replicates are similar.

#### IMPORT DATA ####

SRMAreas <- read.csv("2017-09-06-Gigas-SRM-ReplicatesOnly-PostDilutionCurve-NoPivot-Report.csv", na.strings = "#N/A") #Specify Skyline's special way of designating N/A values
head(SRMAreas) #Confirm import

#### CREATE A MASTER DATAFRAME ####

#I want to merge my Skyline data with sample names, sites, and eelgrass condition to create a master dataframe will all possible information

sequenceFile <- read.csv("2017-07-28-SRM-Samples-Sequence-File.csv", na.strings = "N/A") # Import sequence file
head(sequenceFile) #Confirm import
sequenceFile <- sequenceFile[,c(2,3)] #Keep only the Replicate.Name and Comment columns
names(sequenceFile) <- c("Replicate.Name", "Sample.Number")
head(sequenceFile) #Confirm change
masterSRMData <- merge(x = SRMAreas, y = sequenceFile, by = "Replicate.Name") #Merge the sample names and replicate names. 
head(masterSRMData) #Confirm merge
tail(masterSRMData) #Confirm merge

biologicalReplicates <- read.csv("2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
tail(biologicalReplicates) #Confirm import
masterSRMDataBiologicalReplicates <- merge(x = masterSRMData, y = biologicalReplicates, by = "Sample.Number") #Add biological replicate information to master list. OBLNK2-1 and OBLNK2-2 not included.
head(masterSRMDataBiologicalReplicates)
masterSRMDataBiologicalReplicates <- masterSRMDataBiologicalReplicates[,-8] #Remove TIC Area column since it is empty
head(masterSRMDataBiologicalReplicates) #Confirm change
write.csv(x = masterSRMDataBiologicalReplicates, file = "2017-09-07-Master-SRM-Data-BiologicalReplicates-NoBlanks-NoPivot.csv")

#### MODIFY DATAFRAME ####

masterSRMDataCorrected <- masterSRMData
masterSRMDataCorrected[is.na(masterSRMData)] <- 0 #Replace NAs with 0s and save as a new dataframe
head(masterSRMDataCorrected) #Confirm changes
transform(masterSRMDataCorrected, Peptide.Retention.Time = as.numeric(Peptide.Retention.Time)) #Make sure retention times are recognized as numeric variables
transform(masterSRMDataCorrected, Area = as.numeric(Area)) #Make sure area is recognized as a numeric variable
sapply(masterSRMDataCorrected, class) #Confirm that changes were made
is.numeric(masterSRMDataCorrected$Area) #Confirm changes were made
is.numeric(masterSRMDataCorrected$Peptide.Retention.Time) #Confirm changes were made

SRMDataTargetsOnly <- masterSRMDataCorrected[! masterSRMDataCorrected$Protein.Name %in% "PRTC peptides", ] #Create a dataset without PRTC peptides
head(SRMDataTargetsOnly) #Confirm changes
SRMDataTargetsOnly <- SRMDataTargetsOnly[, -c(5, 7, 9, 10)] #Remove everything that wasn't directly from Skyline and the transition ion
head(SRMDataTargetsOnly) #Confirm changes

#### NON-NORMALIZED NMDS PLOT ####

#Before normalizing, I want to see how my technical replicates cluster

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- SRMDataTargetsOnly[-7]
rownames(area.protID2) <- SRMDataTargetsOnly[,7]

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
area2.t <- t(area.protID2[,1:10])
area2.tra <- (area2.t+1)
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area2.tra, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5))) #Looks like all of the points are on top of eachother.......
legend(0,0.06, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black')) 

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area2.tra dataset.

#Create heatmap
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "fullHeatmap.png")
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)
dev.off()