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
colnames(biomarkerData) <- c("Sample.Number", "Site", "Habitat", "Final Shell Height", "Tissue Mass", "Shell Thickness", "Peak Load", "Shell Strength", "delta C", "Percent C", "delta N", "Percent N", "C to N Ratio")
head(biomarkerData)

#### MERGE DATAFRAMES ####

peptideBiomarkerData <- merge(x = SRMDataNMDSAveragedCorrectedTransposed, y = biomarkerData, by = "Sample.Number") #Merge dataframes by sample number
head(peptideBiomarkerData) #Confirm changes. Peptides are columns 2-38, biomarkers 41-50
#write.csv(peptideBiomarkerData, "2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Biomarker-Scatterplots/2017-11-29-Peptide-and-Biomarker-Data.csv") #Write out dataframe

#### ASSIGN COLORS BY SITE ####

attach(peptideBiomarkerData) #Attach dataframe
peptideBiomarkerData <- peptideBiomarkerData[order(Site),] #Reorder so sites are sorted alphabetically
head(peptideBiomarkerData) #Confirm sorting
tail(peptideBiomarkerData) #Confirm sorting
detach(peptideBiomarkerData)
peptideBiomarkerData$Colors <- c(rep(x = "red", times = sum(peptideBiomarkerData$Site == "CI")),
                                 rep(x = "blue", times = sum(peptideBiomarkerData$Site == "FB")),
                                 rep(x = "magenta", times = sum(peptideBiomarkerData$Site == "PG")),
                                 rep(x = "green", times = sum(peptideBiomarkerData$Site == "SK")),
                                 rep(x = "black", times = sum(peptideBiomarkerData$Site == "WB"))) #Create a color vector
head(peptideBiomarkerData) #Confirm addition
tail(peptideBiomarkerData) #Confirm addition

#### CHANGE WORKING DIRECTORY ####

setwd("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Biomarker-Scatterplots/") #Change working directory
getwd() #Confirm changes

#### MAKE SCATTERPLOTS ####

nPeptides <- 38 #Columns 2-38 are peptide names
for(i in 2:nPeptides) { #For all peptides
  peptideBiomarkerModel <- lm(peptideBiomarkerData[,i] ~ peptideBiomarkerData[,41], na.action = na.omit)
  fileName <- paste(colnames(peptideBiomarkerData)[i], "vs.", colnames(peptideBiomarkerData)[41], ".jpeg")
  jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
  plot(x = peptideBiomarkerData[,41], y = peptideBiomarkerData[,i], xlab = "Final Shell Height", ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerData)[i], "vs.", colnames(peptideBiomarkerData)[41]), cex.main = 2) #Create plot, but do not plot points
  text(x = peptideBiomarkerData[,41], y = peptideBiomarkerData[,i], labels = peptideBiomarkerData$Sample.Number, cex = 0.8, col = peptideBiomarkerData$Colors, font = 2) #Plot sample ID instead of points
  abline(peptideBiomarkerModel) #Plot regression
  legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
  dev.off() #Turn off plotting device
}
