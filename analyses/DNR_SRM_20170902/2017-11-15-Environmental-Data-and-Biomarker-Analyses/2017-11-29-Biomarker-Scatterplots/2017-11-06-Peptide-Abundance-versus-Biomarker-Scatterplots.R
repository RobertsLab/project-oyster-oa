#In this script, I'll regress peptide abundance against different biomarkers to see if the biomarker variables explain variation in abundance.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

SRMDataNMDSAveragedCorrected <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv", header = TRUE) #Import modified dataset. This dataset has rownames as the first column, column names as sample IDs. Area data is averaged and corrected (no NAs)
rownames(SRMDataNMDSAveragedCorrected) <- SRMDataNMDSAveragedCorrected[,1]
SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveragedCorrected[,-1] #Remove first column of rownames
head(SRMDataNMDSAveragedCorrected) #Confirm import.

#### REFORMAT DATA ####

#Transpose dataframe
SRMDataNMDSAveragedCorrectedTransposed <- data.frame(t(SRMDataNMDSAveragedCorrected)) #Transpose the data
SRMDataNMDSAveragedCorrectedTransposed$Sample.Number <- rownames(SRMDataNMDSAveragedCorrectedTransposed) #Save rownames as a new column
head(SRMDataNMDSAveragedCorrectedTransposed) #Confirm changes

#### IMPORT BIOMARKER DATA ####

biomarkerData <- read.csv("../../data/DNR/2017-11-21-Alex-Data-Yaamini-Samples-Only.csv") #Import dataset
biomarkerData <- biomarkerData[,-c(1:3, 5, 8:11, 13:14)] #Remove empty columns and columns I don't need (Round, Spp, Bag, Rep.x, FAvial, topValve, bothValves)
head(biomarkerData) #Confirm changes
colnames(biomarkerData) <- c("Sample.Number", "Site", "Habitat", "finalHeight", "tissueMass", "shellThickness", "peakLoad", "Strength", "delC", "percentC", "delN", "percentN", "CNRatio")
head(biomarkerData)

#### MERGE DATAFRAMES ####

peptideBiomarkerData <- merge(x = SRMDataNMDSAveragedCorrectedTransposed, y = biomarkerData, by = "Sample.Number") #Merge dataframes by sample number
head(peptideBiomarkerData)
#write.csv(peptideBiomarkerData, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Biomarker-Scatterplots/2017-11-29-Peptide-and-Biomarker-Data.csv") #Write out dataframe
