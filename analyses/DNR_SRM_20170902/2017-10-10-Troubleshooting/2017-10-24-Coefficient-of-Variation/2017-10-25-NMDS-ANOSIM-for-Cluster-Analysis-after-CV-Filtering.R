#Steven went through my technical replicates and calculated a coefficient of variation for each transition. He removed any transitions that had CVs greater than 20. Using this dataset, I'll average my technical replicates and proceed with an NMDS and ANOSIM.

#### SET WORKING DIRECTORY ####

setwd("../..")
getwd()

#### IMPORT DATA ####

SRMModifiedAreas <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-24-Norm-modSR.csv", header = TRUE) #Import Steven's modified dataset
head(SRMModifiedAreas) #Confirm import. There are eight columns. The first column and the last two columns are null, and the coefficient of variation column is not needed.
SRMModifiedAreas <- SRMModifiedAreas[,-c(1, 6:8)] #Remove unnecessary columns
head(SRMModifiedAreas) #Confirm changes. Now I only have Protein.Name, Sample, Replicate1 and Replicate2 (peak areas from Skyline, which are a proxy for protein abundance)

#### AVERAGE TECHNICAL REPLICATES ####

SRMAveragedAreas <- SRMModifiedAreas #Dupliate dataframe
SRMAveragedAreas$Average.Area <- ((SRMAveragedAreas$Replicate1 + SRMAveragedAreas$Replicate2)/2) #Average peak areas and save as a new column
head(SRMAveragedAreas) #Confirm changes
