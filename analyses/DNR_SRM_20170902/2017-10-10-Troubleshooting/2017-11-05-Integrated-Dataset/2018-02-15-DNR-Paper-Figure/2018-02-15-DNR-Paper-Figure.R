#In this script, I'll make figure options for the DNR paper.

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

averagePeptideData <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-DNR-Paper-Figure/2018-02-16-Average-Peptide-Data-by-Site.csv", header = TRUE) #Import average peptide dataset
head(averagePeptideData) #Confirm import
averagePeptideData <- averagePeptideData[,-1] #Remove first column
head(averagePeptideData) #Confirm column deletion

#### SPECIFY COLORS AND SHAPES ####
#I want each protein to have the same color, with the corresponding peptides as different shapes. I also think that proteins with similar functions should have similar colors.
averagePeptideData$shapes <- c(16, 17, 15, 16, 17, 16, 17, 15, 16, 17, 16, 17, 15, 16, 17, 15, 16, 16, 17, 16, 17, 16, 17, 16, 17, 15, 16, 17, 16, 17, 15, 16, 17, 15, 16, 17, 15) #Specify shapes. 16 = circle, 17 = triangle, 15 = square
head(averagePeptideData) #Confirm addition
averagePeptideData$colors <- c(rep("goldenrod1", times = 3), rep("turquoise1", times = 2), rep("deepskyblue1", times = 3), rep("seagreen1", times = 2), rep("purple1", times = 3), rep("dodgerblue1", times = 3), rep("dodgerblue3", times = 1), rep("sandybrown", times = 2), rep("orangered", times = 2), rep("palegreen3", times = 2), rep("chocolate4", times = 3), rep("violet", times = 2), rep("olivedrab2", times = 3), rep("lawngreen", times = 3), rep("springgreen4", times = 3))

#### CREATE PLOT WITH ALL DATA ####

#jpeg("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-DNR-Paper-Figure/2018-02-16-All-Peptide-Abundances-Across-Sites.jpeg", height = 1000, width = 1000) #Save file
dotchart(x = averagePeptideData$averageNormalizedAbundance, groups = averagePeptideData$site, pch = averagePeptideData$shapes, color = averagePeptideData$colors, pt.cex = 2.5, main = "Peptide Abundances Across Sites", xlab = "Normalized Peptide Abundance", ylab = "Sites", cex.main = 3, cex.lab = 1.5) #Create plot with all peptide abundance data
#dev.off() #Turn off plotting device

#### ISOLATE ONLY DIFFERENTIALLY EXPRESSED PEPTIDES ####

diffExpPeptides <- peptideNames[c(4, 6:10, 14:15, 18:19, 22, 24, 26)] #Copy differentially expressed peptides into a new vector
diffExpPeptides <- data.frame("peptide" = diffExpPeptides,
                              "temp" = rep(0, times = length(diffExpPeptides))) #Add a temporary column and make diffExpPepties a dataframe
head(diffExpPeptides) #Confirm column addition
diffExpPeptides <- merge(x = diffExpPeptides, y = averagePeptideData, by = "peptide") #Merge by peptide name
diffExpPeptides <- diffExpPeptides[,-2] #Remove "temp" column
head(diffExpPeptides) #Confirm merge and column removal

attach(diffExpPeptides)
diffExpPeptides <- diffExpPeptides[order(site),] #Reorder so sites are sorted alphabetically
detach(diffExpPeptides)
head(diffExpPeptides) #Confirm sort

#### CREATE PLOT WITH DIFFERENTIALLY EXPRESSED PEPTIDES ####
#jpeg("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-DNR-Paper-Figure/2018-02-18-Differentially-Expressed-Peptides-Across-Sites.jpeg", height = 1000, width = 1000) #Save file
dotchart(x = diffExpPeptides$averageNormalizedAbundance, groups = diffExpPeptides$site, pch = diffExpPeptides$shapes, color = diffExpPeptides$colors, pt.cex = 2.5, main = "Differentially Expressed Peptides Across Sites", xlab = "Normalized Peptide Abundance", ylab = "Sites", cex.main = 3, cex.lab = 1.5) #Create plot with differentially expressed peptide abundance data
#dev.off() #Turn off plotting device

