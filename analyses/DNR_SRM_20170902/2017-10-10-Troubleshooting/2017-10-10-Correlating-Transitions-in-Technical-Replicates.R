#In this script, I'll assess if there are certain transitions contributing to the poor technical replication of my samples.

#### IMPORT DATA ####

SRMAreas <- read.csv("2017-09-12-Gigas-SRM-ReplicatesOnly-PostDilutionCurve-NoPivot-RevisedSettings-Report.csv", na.strings = "#N/A") #Specify Skyline's special way of designating N/A values
head(SRMAreas) #Confirm import
tail(SRMAreas) #Confirm import

#### CREATE A MASTER DATAFRAME ####

#I want to merge my Skyline data with sample names, sites, and eelgrass condition to create a master dataframe will all possible information

sequenceFile <- read.csv("2017-07-28-SRM-Samples-Sequence-File.csv", na.strings = "N/A") # Import sequence file
head(sequenceFile) #Confirm import
sequenceFile <- sequenceFile[,c(2,3,8)] #Keep the Replicate.Name, Comment and TIC columns
names(sequenceFile) <- c("Replicate.Name", "Sample.Number", "TIC")
head(sequenceFile) #Confirm change
masterSRMData <- merge(x = SRMAreas, y = sequenceFile, by = "Replicate.Name") #Merge the sample names and replicate names to use for analysis.
head(masterSRMData) #Confirm merge
tail(masterSRMData) #Confirm merge

biologicalReplicates <- read.csv("2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
tail(biologicalReplicates) #Confirm import
masterSRMDataBiologicalReplicates <- merge(x = masterSRMData, y = biologicalReplicates, by = "Sample.Number") #Add biological replicate information to master list.
head(masterSRMDataBiologicalReplicates) #Confirm change
#write.csv(x = masterSRMDataBiologicalReplicates, file = "2017-09-07-Master-SRM-Data-BiologicalReplicates-NoBlanks-NoPivot.csv") #Write out master dataframe

#### SUBSET TARGET PROTEIN DATA ####

#I want only the protein/peptide/transition information and peak area

SRMDataTargetsOnly <- masterSRMDataBiologicalReplicates #Duplicate master list into a new dataframe
head(SRMDataTargetsOnly) #Confirm copy
tail(SRMDataTargetsOnly) #Confirm copy
SRMDataTargetsOnly <- SRMDataTargetsOnly[,-c(2, 5, 7, 10, 11)] #Remove extraneous columns: Replicate.Name, Transition, Peptide.Retention.Time, Site, Eelgrass
head(SRMDataTargetsOnly) #Confirm column removal
SRMDataPRTCOnly <- SRMDataTargetsOnly[SRMDataTargetsOnly$Protein.Name %in% "PRTC peptides", ] #Save PRTC peptide data as a new dataframe
SRMDataTargetsOnly <- SRMDataTargetsOnly[! SRMDataTargetsOnly$Protein.Name %in% "PRTC peptides", ] #Remove PRTC peptide data from target protein dataframe
head(SRMDataTargetsOnly) #Confirm removal
transform(SRMDataTargetsOnly, Area = as.numeric(Area)) #Make sure Area is recognized as a numeric variable
is.numeric(SRMDataTargetsOnly$Area) #Confirm change
transform(SRMDataTargetsOnly, TIC = as.numeric(TIC)) #Make sure TIC is recognized as a numeric variable
is.numeric(SRMDataTargetsOnly$TIC) #Confirm change

#### REFORMAT TARGET ONLY DATAFRAME ####

#First I will make a new dataframe be Protein/Peptides/Transitions, with the column names as the sample number

SRMDataTargetsOnly <- SRMDataTargetsOnly[,-6] #Remove TIC column
head(SRMDataTargetsOnly) #Confirm creation

#My first step is to change my dataframe from long to wide (i.e. cast it)
library(reshape2) #Instal package to pivot table
SRMDataTargetsOnlyPivoted <- dcast(SRMDataTargetsOnly, Protein.Name + Peptide.Sequence + Fragment.Ion ~ Sample.Number) #Cast table! Protein/Peptides/Transitions remain as columns with Sample Number as column headers. Normalized.Area used as value column by default.
head(SRMDataTargetsOnlyPivoted) #Confirm cast.
SRMDataTargetsOnlyPivoted$RowNames <- paste(SRMDataTargetsOnlyPivoted$Protein.Name, SRMDataTargetsOnlyPivoted$Peptide.Sequence, SRMDataTargetsOnlyPivoted$Fragment.Ion) #Merge Protein, Peptide and Transition information into one column
head(SRMDataTargetsOnlyPivoted) #Confirm column merge
SRMDataTargetsOnlyPivoted <- SRMDataTargetsOnlyPivoted[,-c(1:3)] #Remove unmerged columns
head(SRMDataTargetsOnlyPivoted) #Confirm column removal

SRMDataTargetsOnlyPivotedCorrected <- SRMDataTargetsOnlyPivoted #Duplicate dataframe
SRMDataTargetsOnlyPivotedCorrected[is.na(SRMDataTargetsOnlyPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataTargetsOnlyPivotedCorrected) #Confirm there are no NAs
rownames(SRMDataTargetsOnlyPivotedCorrected) <- SRMDataTargetsOnlyPivotedCorrected$RowNames #Set RowNames column as dataframe rownames
SRMDataTargetsOnlyPivotedCorrected <- subset(SRMDataTargetsOnlyPivotedCorrected, select = -c(RowNames)) #Remove RowNames column
head(SRMDataTargetsOnlyPivotedCorrected) #Confirm changes

#Now I'll take my dataframe and split it into two: one for each batch of technical replicates.

SRMDataTargetsReplicateOne <- SRMDataTargetsOnlyPivotedCorrected[, c(seq(from = 1, to = (length(SRMDataTargetsOnlyPivotedCorrected) - 1), by = 2))] #Subset all odd columns (first replicate)
colnames(SRMDataTargetsReplicateOne) #Confirm subset
SRMDataTargetsReplicateTwo <- SRMDataTargetsOnlyPivotedCorrected[, c(seq(from = 2, to = length(SRMDataTargetsOnlyPivotedCorrected), by = 2))] #Subset all even columns (first replicate)
colnames(SRMDataTargetsReplicateTwo) #Confirm subset

#Finally, I'll transpose each dataframe. The resulting dataframes will have transitions in the columns and samples in the rows.

SRMDataTransposedReplicateOne <- t(SRMDataTargetsReplicateOne) #Transpose Replicate 1 dataframe
head(SRMDataTransposedReplicateOne) #Confirm transposition
SRMDataTransposedReplicateTwo <- t(SRMDataTargetsReplicateTwo) #Transpose Replicate 2 dataframe
head(SRMDataTransposedReplicateTwo) #Confirm transposition

#### CREATE PLOTS FOR EACH TRANSITION ####

plot(x= SRMDataTargetsOnlyPivotedCorrected$`O137-1`, y = SRMDataTargetsOnlyPivotedCorrected$`O137-2`) #x = first column of first dataframe, y = first column of second dataframe

#### REFORMAT PRTC ONLY DATAFRAME ####

transform(SRMDataPRTCOnly, Area = as.numeric(Area)) #Make sure Area is recognized as a numeric variable
is.numeric(SRMDataPRTCOnly$Area) #Confirm change
transform(SRMDataPRTCOnly, TIC = as.numeric(TIC)) #Make sure TIC is recognized as a numeric variable
is.numeric(SRMDataPRTCOnly$TIC) #Confirm change