#In this script, I'll depict normalized protein area across samples as bar charts after integrating my transitions. These boxplots will be at the peptide level and only differentiated by sites. This is because my NMDS analyses did not show any significant differences in clustering based on habitat type.

#### SET WORKING DIRECTORY ####
setwd("../../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

SRMDataNMDSAveragedCorrected <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv", header = TRUE) #Import modified dataset. This dataset has rownames as the first column, column names as sample IDs. Area data is averaged and corrected (no NAs)
rownames(SRMDataNMDSAveragedCorrected) <- SRMDataNMDSAveragedCorrected[,1]
SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveragedCorrected[,-1] #Remove first column of rownames
head(SRMDataNMDSAveragedCorrected) #Confirm import.

#### REFORMAT DATA ####

SRMDataNMDSAveragedCorrectedTransposed <- data.frame(t(SRMDataNMDSAveragedCorrected)) #Transpose the data
SRMDataNMDSAveragedCorrectedTransposed$Sample.Number <- rownames(SRMDataNMDSAveragedCorrectedTransposed) #Save rownames as a new column
head(SRMDataNMDSAveragedCorrectedTransposed) #Confirm changes

biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
colnames(biologicalReplicates) <- c("Sample.Number", "Site", "Eelgrass.Condition")
rownames(biologicalReplicates) <- biologicalReplicates$Sample.Number #Set sample number as row names
head(biologicalReplicates) #Confirm changes
biologicalReplicates$Site <- factor(biologicalReplicates$Site) #Remove 0 as a factor
biologicalReplicates$Eelgrass.Condition <- factor(biologicalReplicates$Eelgrass.Condition) #Remove 0 as a factor
str(biologicalReplicates) #Confirm factor reset

boxplotData <- merge(x = biologicalReplicates, y = SRMDataNMDSAveragedCorrectedTransposed, by = "Sample.Number") #Merge together
head(boxplotData) #Confirm merge
rownames(boxplotData) <- boxplotData$Sample.Number #Set sample number as row names
boxplotData <- boxplotData[-1] #Remove Sample.Number column
head(boxplotData) #Confirm changes

attach(boxplotData)
boxplotData <- boxplotData[order(Site),] #Reorder so sites are sorted alphabetically
detach(boxplotData)
boxplotData <- boxplotData[-2] #Remove habitat information
head(boxplotData) #Confirm changes

#### BREAKUP INDIVIDUAL DATASETS ####

caseInletData <- subset(x = boxplotData, subset = boxplotData$Site == "CI") #Subset Case Inlet data
fidalgoBayData <- subset(x = boxplotData, subset = boxplotData$Site == "FB") #Subset Fidalgo Bay data
portGambleData <- subset(x = boxplotData, subset = boxplotData$Site == "PG") #Subset Port Gamble Bay data
skokomishRiverData <- subset(x = boxplotData, subset = boxplotData$Site == "SK") #Subset Skokomish River Delta data
willapaBayData <- subset(x = boxplotData, subset = boxplotData$Site == "WB") #Subset Willapa Bay data

caseInletData <- caseInletData[-1] #Remove site classification
fidalgoBayData <- fidalgoBayData[-1] #Remove site classification
portGambleData <- portGambleData[-1] #Remove site classification
skokomishRiverData <- skokomishRiverData[-1] #Remove site classification
willapaBayData <- willapaBayData[-1] #Remove site classification

#### AVERAGE ACROSS SAMPLES ####

peptideNames <- colnames(caseInletData) #Isolate peptide names

#Case Inlet
caseInletAverages <- rep(0, times = length(peptideNames)) #Create an empty vector to store data
for(i in 1:length(peptideNames)) {
  caseInletAverages[i] <- mean(caseInletData[,i])
} #Average each column and save it in the caseInletAverages vector

caseInletAverages <- data.frame("peptide" = peptideNames,
                                "site" = rep("CI", times = length(peptideNames)),
                                "averageNormalizedAbundance" = caseInletAverages) #Create a new dataframe with peptide names, site, and average normalized protein abundances
head(caseInletAverages) #Confirm dataframe creation

#Fidalgo Bay
fidalgoBayAverages <- rep(0, times = length(peptideNames)) #Create an empty vector to store data
for(i in 1:length(peptideNames)) {
  fidalgoBayAverages[i] <- mean(fidalgoBayData[,i])
} #Average each column and save it in the fidalgoBayAverages vector

fidalgoBayAverages <- data.frame("peptide" = peptideNames,
                                "site" = rep("FB", times = length(peptideNames)),
                                "averageNormalizedAbundance" = fidalgoBayAverages) #Create a new dataframe with peptide names, site, and average normalized protein abundances
head(fidalgoBayAverages) #Confirm dataframe creation

#Port Gamble Bay
portGambleAverages <- rep(0, times = length(peptideNames)) #Create an empty vector to store data
for(i in 1:length(peptideNames)) {
  portGambleAverages[i] <- mean(portGambleData[,i])
} #Average each column and save it in the portGambleAverages vector

portGambleAverages <- data.frame("peptide" = peptideNames,
                                "site" = rep("PG", times = length(peptideNames)),
                                "averageNormalizedAbundance" = portGambleAverages) #Create a new dataframe with peptide names, site, and average normalized protein abundances
head(portGambleAverages) #Confirm dataframe creation

#Skokomish River Delta
skokomishRiverAverages <- rep(0, times = length(peptideNames)) #Create an empty vector to store data
for(i in 1:length(peptideNames)) {
  skokomishRiverAverages[i] <- mean(skokomishRiverData[,i])
} #Average each column and save it in the skokomishRiverAverages vector

skokomishRiverAverages <- data.frame("peptide" = peptideNames,
                                "site" = rep("SK", times = length(peptideNames)),
                                "averageNormalizedAbundance" = skokomishRiverAverages) #Create a new dataframe with peptide names, site, and average normalized protein abundances
head(skokomishRiverAverages) #Confirm dataframe creation

#Willapa Bay
willapaBayAverages <- rep(0, times = length(peptideNames)) #Create an empty vector to store data
for(i in 1:length(peptideNames)) {
  willapaBayAverages[i] <- mean(willapaBayData[,i])
} #Average each column and save it in the willapaBayAverages vector

willapaBayAverages <- data.frame("peptide" = peptideNames,
                                "site" = rep("WB", times = length(peptideNames)),
                                "averageNormalizedAbundance" = willapaBayAverages) #Create a new dataframe with peptide names, site, and average normalized protein abundances
head(willapaBayAverages) #Confirm dataframe creation

averagePeptideData <- rbind(caseInletAverages, fidalgoBayAverages, portGambleAverages, skokomishRiverAverages, willapaBayAverages) #Merge all averaged peptide data into a single dataframe
#write.csv(averagePeptideData, "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-DNR-Paper-Figure/2018-02-16-Average-Peptide-Data-by-Site.csv") #Wrote out table for future analyses
