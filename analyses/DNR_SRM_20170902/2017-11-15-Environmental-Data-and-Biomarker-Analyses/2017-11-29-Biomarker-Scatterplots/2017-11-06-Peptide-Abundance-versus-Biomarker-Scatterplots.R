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

nPeptides <- 38 #Columns 2-38 are peptides
nBiomarkers <- 50 #Columns 41-50 are biomarkers
for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerData[,i] ~ peptideBiomarkerData[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerData)[i], "vs.", colnames(peptideBiomarkerData)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerData[,j], y = peptideBiomarkerData[,i], xlab = colnames(peptideBiomarkerData)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerData)[i], "vs.", colnames(peptideBiomarkerData)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerData[,j], y = peptideBiomarkerData[,i], labels = peptideBiomarkerData$Sample.Number, cex = 2, col = peptideBiomarkerData$Colors, font = 2) #Plot sample ID instead of points
    abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
}

#### PUT IMPORTANT INFORMATION INTO NEW DATAFRAME ####

comparisonStatistics <- data.frame("ComparisonName" = rep(x = 0, times = 370),
                                   "Peptide" = rep(x = 0, times = 370),
                                   "Biomarker" = rep(x = 0, times = 370),
                                   "R.Squared" = rep(x = 0, times = 370),
                                   "Slope" = rep(x = 0, times = 370)) #Create empty dataframe to hold all results
head(comparisonStatistics) #Confirm changes

nComparisons <- length(comparisonStatistics$ComparisonName) #Calculate number of comparisons
for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    testName <- paste(colnames(peptideBiomarkerData)[i], "vs.", colnames(peptideBiomarkerData)[j]) #Assign comparison name
    print(testName)
  }
}

peptideBiomarkerModel <- lm(peptideBiomarkerData[,2] ~ peptideBiomarkerData[,41], na.action = na.omit) #Make model
comparisonStatistics$ComparisonName[1] <- paste(colnames(peptideBiomarkerData)[2], "vs.", colnames(peptideBiomarkerData)[41]) #Assign comparison name
comparisonStatistics$Peptide[1] <-  colnames(peptideBiomarkerData)[2] #Specify peptide
comparisonStatistics$Biomarker[1] <- colnames(peptideBiomarkerData)[41] #Specify biomarker
comparisonStatistics$R.Squared[1] <- format(summary(peptideBiomarkerModel)$adj.r.squared, digits = 4) #Adjusted R-squared
comparisonStatistics$Slope[1] <- format(summary(peptideBiomarkerModel)$coeff[2], digits = 4) #Regression's slope

#### MAKE SITE SPECIFIC SCATTERPLOTS ####

#Case Inlet
peptideBiomarkerCaseInlet <- subset(peptideBiomarkerData, subset = peptideBiomarkerData$Site == "CI") #Create Case Inlet subset
head(peptideBiomarkerCaseInlet) #Confirm subset

setwd("2017-12-01-Case-Inlet-Scatterplots") #Change working directory
getwd() #Confirm changes

for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerCaseInlet[,i] ~ peptideBiomarkerCaseInlet[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerCaseInlet)[i], "vs.", colnames(peptideBiomarkerCaseInlet)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerCaseInlet[,j], y = peptideBiomarkerCaseInlet[,i], xlab = colnames(peptideBiomarkerCaseInlet)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerCaseInlet)[i], "vs.", colnames(peptideBiomarkerCaseInlet)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerCaseInlet[,j], y = peptideBiomarkerCaseInlet[,i], labels = peptideBiomarkerCaseInlet$Sample.Number, cex = 2, font = 2) #Plot sample ID instead of points
    abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
}

#Fidalgo Bay

peptideBiomarkerFidalgoBay <- subset(peptideBiomarkerData, subset = peptideBiomarkerData$Site == "FB") #Create Fidalgo Bay subset
head(peptideBiomarkerFidalgoBay) #Confirm subset

setwd("../2017-12-01-Fidalgo-Bay-Scatterplots/") #Change working directory
getwd() #Confirm changes

for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerFidalgoBay[,i] ~ peptideBiomarkerFidalgoBay[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerFidalgoBay)[i], "vs.", colnames(peptideBiomarkerFidalgoBay)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerFidalgoBay[,j], y = peptideBiomarkerFidalgoBay[,i], xlab = colnames(peptideBiomarkerFidalgoBay)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerFidalgoBay)[i], "vs.", colnames(peptideBiomarkerFidalgoBay)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerFidalgoBay[,j], y = peptideBiomarkerFidalgoBay[,i], labels = peptideBiomarkerFidalgoBay$Sample.Number, cex = 2, font = 2) #Plot sample ID instead of points
    #abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
} #Couldn't plot regression on all graphs, might have been an issue with a and b being finite for model

#Port Gamble Bay

peptideBiomarkerPortGamble <- subset(peptideBiomarkerData, subset = peptideBiomarkerData$Site == "PG") #Create Port Gamble subset
head(peptideBiomarkerPortGamble) #Confirm subset

setwd("../2017-12-01-Port-Gamble-Scatterplots/") #Change working directory
getwd() #Confirm changes

for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerPortGamble[,i] ~ peptideBiomarkerPortGamble[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerPortGamble)[i], "vs.", colnames(peptideBiomarkerPortGamble)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerPortGamble[,j], y = peptideBiomarkerPortGamble[,i], xlab = colnames(peptideBiomarkerPortGamble)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerPortGamble)[i], "vs.", colnames(peptideBiomarkerPortGamble)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerPortGamble[,j], y = peptideBiomarkerPortGamble[,i], labels = peptideBiomarkerPortGamble$Sample.Number, cex = 2, font = 2) #Plot sample ID instead of points
    abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
}

#Skokomish River Delta

peptideBiomarkerSkokomishRiver <- subset(peptideBiomarkerData, subset = peptideBiomarkerData$Site == "SK") #Create Skokomish River subset
head(peptideBiomarkerSkokomishRiver) #Confirm subset

setwd("../2017-12-01-Skokomish-River-Scatterplots/") #Change working directory
getwd() #Confirm changes

for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerSkokomishRiver[,i] ~ peptideBiomarkerSkokomishRiver[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerSkokomishRiver)[i], "vs.", colnames(peptideBiomarkerSkokomishRiver)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerSkokomishRiver[,j], y = peptideBiomarkerSkokomishRiver[,i], xlab = colnames(peptideBiomarkerSkokomishRiver)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerSkokomishRiver)[i], "vs.", colnames(peptideBiomarkerSkokomishRiver)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerSkokomishRiver[,j], y = peptideBiomarkerSkokomishRiver[,i], labels = peptideBiomarkerSkokomishRiver$Sample.Number, cex = 2, font = 2) #Plot sample ID instead of points
    abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
}

#Willapa Bay

peptideBiomarkerSkokomishRiver <- subset(peptideBiomarkerData, subset = peptideBiomarkerData$Site == "WB") #Create Willapa Bay subset
head(peptideBiomarkerSkokomishRiver) #Confirm subset

setwd("../2017-12-01-Skokomish-River-Scatterplots/") #Change working directory
getwd() #Confirm changes

for(i in 2:nPeptides) { #For all peptides
  for (j in 41:nBiomarkers) { #For all biomarkers
    peptideBiomarkerModel <- lm(peptideBiomarkerSkokomishRiver[,i] ~ peptideBiomarkerSkokomishRiver[,j], na.action = na.omit)
    fileName <- paste(colnames(peptideBiomarkerSkokomishRiver)[i], "vs.", colnames(peptideBiomarkerSkokomishRiver)[j], ".jpeg")
    jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
    plot(x = peptideBiomarkerSkokomishRiver[,j], y = peptideBiomarkerSkokomishRiver[,i], xlab = colnames(peptideBiomarkerSkokomishRiver)[j], ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideBiomarkerSkokomishRiver)[i], "vs.", colnames(peptideBiomarkerSkokomishRiver)[j]), cex.main = 1.75) #Create plot, but do not plot points
    text(x = peptideBiomarkerSkokomishRiver[,j], y = peptideBiomarkerSkokomishRiver[,i], labels = peptideBiomarkerSkokomishRiver$Sample.Number, cex = 2, font = 2) #Plot sample ID instead of points
    abline(peptideBiomarkerModel) #Plot regression
    legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideBiomarkerModel)$adj.r.squared, digits=4))) #Plot R-squared value
    dev.off() #Turn off plotting device
  }
}
