#In this script, I'll test for bare vs. eelgrass protein expression differences at the individual site level. Brent and Emma suggested I do this to double check that there is no eelgrass effect, since I did not see a global eelgrass effect.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
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

#### BREAKUP INDIVIDUAL DATASETS ####

caseInletData <- subset(x = boxplotData, subset = boxplotData$Site == "CI") #Subset Case Inlet data
fidalgoBayData <- subset(x = boxplotData, subset = boxplotData$Site == "FB") #Subset Fidalgo Bay data
portGambleData <- subset(x = boxplotData, subset = boxplotData$Site == "PG") #Subset Port Gamble Bay data
skokomishRiverData <- subset(x = boxplotData, subset = boxplotData$Site == "SK") #Subset Skokomish River Delta data
willapaBayData <- subset(x = boxplotData, subset = boxplotData$Site == "WB") #Subset Willapa Bay data

#### ASSIGN FILENAMES ####

caseInletFilenames <- data.frame(protein = colnames(caseInletData),
                               modifier = rep("BvECaseInlet.jpeg", length(caseInletData))) #Make filename sheet
caseInletFilenames$siteFilenames <- paste(caseInletFilenames$protein, caseInletFilenames$modifier) #Make a new column for the actual filenames
head(caseInletFilenames) #Confirm changes

fidalgoBayFilenames <- data.frame(protein = colnames(fidalgoBayData),
                                 modifier = rep("BvEFidalgoBay.jpeg", length(fidalgoBayData))) #Make filename sheet
fidalgoBayFilenames$siteFilenames <- paste(fidalgoBayFilenames$protein, fidalgoBayFilenames$modifier) #Make a new column for the actual filenames
head(fidalgoBayFilenames) #Confirm changes

portGambleFilenames <- data.frame(protein = colnames(portGambleData),
                                 modifier = rep("BvEPortGamble.jpeg", length(portGambleData))) #Make filename sheet
portGambleFilenames$siteFilenames <- paste(portGambleFilenames$protein, portGambleFilenames$modifier) #Make a new column for the actual filenames
head(portGambleFilenames) #Confirm changes

skokomishRiverFilenames <- data.frame(protein = colnames(skokomishRiverData),
                                 modifier = rep("BvESkokomishRiver.jpeg", length(skokomishRiverData))) #Make filename sheet
skokomishRiverFilenames$siteFilenames <- paste(skokomishRiverFilenames$protein, skokomishRiverFilenames$modifier) #Make a new column for the actual filenames
head(skokomishRiverFilenames) #Confirm changes

willapaBayFilenames <- data.frame(protein = colnames(willapaBayData),
                                 modifier = rep("BvEWillapaBay.jpeg", length(willapaBayData))) #Make filename sheet
willapaBayFilenames$siteFilenames <- paste(willapaBayFilenames$protein, willapaBayFilenames$modifier) #Make a new column for the actual filenames
head(willapaBayFilenames) #Confirm changes

#### CHANGE WORKING DIRECTORY ####

setwd("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-Individual-Site-Eelgrass-Differences/")
getwd()

#### BOXPLOTS ####

#Case Inlet
nPeptides <- (length(caseInletData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- caseInletFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(caseInletData[,i] ~ caseInletData$Eelgrass.Condition, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(caseInletData[,i] ~ caseInletData$Eelgrass.Condition, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(caseInletData[,i] ~ caseInletData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(caseInletFilenames$protein[i])
  dev.off() #Close file
}

#Fidgalo Bay
nPeptides <- (length(fidalgoBayData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- fidalgoBayFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(fidalgoBayData[,i] ~ fidalgoBayData$Eelgrass.Condition, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(fidalgoBayData[,i] ~ fidalgoBayData$Eelgrass.Condition, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(fidalgoBayData[,i] ~ fidalgoBayData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(fidalgoBayFilenames$protein[i])
  dev.off() #Close file
}

#Port Gamble Bay
nPeptides <- (length(portGambleData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- portGambleFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(portGambleData[,i] ~ portGambleData$Eelgrass.Condition, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(portGambleData[,i] ~ portGambleData$Eelgrass.Condition, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(portGambleData[,i] ~ portGambleData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(portGambleFilenames$protein[i])
  dev.off() #Close file
}

#Skokomish River Delta
nPeptides <- (length(skokomishRiverData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- skokomishRiverFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(skokomishRiverData[,i] ~ skokomishRiverData$Eelgrass.Condition, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(skokomishRiverData[,i] ~ skokomishRiverData$Eelgrass.Condition, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(skokomishRiverData[,i] ~ skokomishRiverData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(skokomishRiverFilenames$protein[i])
  dev.off() #Close file
}

#Willapa Bay
nPeptides <- (length(willapaBayData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- willapaBayFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(willapaBayData[,i] ~ willapaBayData$Eelgrass.Condition, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(willapaBayData[,i] ~ willapaBayData$Eelgrass.Condition, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(willapaBayData[,i] ~ willapaBayData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(willapaBayFilenames$protein[i])
  dev.off() #Close file
}

#### CREATE NEW DATA TABLE ####

#Case Inlet
caseInletANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(caseInletData),
                                    "ANOVA.Fstatistic" = rep(x = 0, times = length(caseInletData)),
                                    "ANOVA.pvalue" = rep(x = 0, times = length(caseInletData))) #Create a dataframe to hold all results
caseInletANOVATukeyResults <- caseInletANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(caseInletANOVATukeyResults) #Confirm changes

#Fidalgo Bay
fidalgoBayANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(fidalgoBayData),
                                         "ANOVA.Fstatistic" = rep(x = 0, times = length(fidalgoBayData)),
                                         "ANOVA.pvalue" = rep(x = 0, times = length(fidalgoBayData))) #Create a dataframe to hold all results
fidalgoBayANOVATukeyResults <- fidalgoBayANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(fidalgoBayANOVATukeyResults) #Confirm changes

#Port Gamble
portGambleANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(portGambleData),
                                         "ANOVA.Fstatistic" = rep(x = 0, times = length(portGambleData)),
                                         "ANOVA.pvalue" = rep(x = 0, times = length(portGambleData))) #Create a dataframe to hold all results
portGambleANOVATukeyResults <- portGambleANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(portGambleANOVATukeyResults) #Confirm changes

#Skokomish River Delta
skokomishRiverANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(skokomishRiverData),
                                         "ANOVA.Fstatistic" = rep(x = 0, times = length(skokomishRiverData)),
                                         "ANOVA.pvalue" = rep(x = 0, times = length(skokomishRiverData))) #Create a dataframe to hold all results
skokomishRiverANOVATukeyResults <- skokomishRiverANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(skokomishRiverANOVATukeyResults) #Confirm changes

#Willapa Bay
willapaBayANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(willapaBayData),
                                         "ANOVA.Fstatistic" = rep(x = 0, times = length(willapaBayData)),
                                         "ANOVA.pvalue" = rep(x = 0, times = length(willapaBayData))) #Create a dataframe to hold all results
willapaBayANOVATukeyResults <- willapaBayANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(willapaBayANOVATukeyResults) #Confirm changes

#### ADJUST P-VALUES FOR MULTIPLE COMPARISONS ####
#FDR = 0.10

#Case Inlet
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  caseInletANOVA <- aov(caseInletData[,i] ~ caseInletData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between eelgrass conditions
  caseInletANOVATukeyResults[(i-2), 2] <- summary(caseInletANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  caseInletANOVATukeyResults[(i-2), 3] <- summary(caseInletANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
} #Add ANOVA F statistics and p-values to the table
head(caseInletANOVATukeyResults) #Confirm that tests were completed
caseInletANOVATukeyResults$ANOVA.adjusted.BH.pvalue <- p.adjust(p = caseInletANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values using B-H method and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(caseInletANOVATukeyResults) #Confirm addition
#write.csv(caseInletANOVATukeyResults, "2018-02-15-OneWayANOVA-TukeyHSD-by-Habitat-CaseInlet-pValues.csv") #Wrote out table for future analyses

#Fidalgo Bay
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fidalgoBayANOVA <- aov(fidalgoBayData[,i] ~ fidalgoBayData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between eelgrass conditions
  fidalgoBayANOVATukeyResults[(i-2), 2] <- summary(fidalgoBayANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  fidalgoBayANOVATukeyResults[(i-2), 3] <- summary(fidalgoBayANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
} #Add ANOVA F statistics and p-values to the table
head(fidalgoBayANOVATukeyResults) #Confirm that tests were completed
fidalgoBayANOVATukeyResults$ANOVA.adjusted.BH.pvalue <- p.adjust(p = fidalgoBayANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values using B-H method and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(fidalgoBayANOVATukeyResults) #Confirm addition
#write.csv(fidalgoBayANOVATukeyResults, "2018-02-15-OneWayANOVA-TukeyHSD-by-Habitat-FidalgoBay-pValues.csv") #Wrote out table for future analyses

#Port Gamble Bay
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  portGambleANOVA <- aov(portGambleData[,i] ~ portGambleData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between eelgrass conditions
  portGambleANOVATukeyResults[(i-2), 2] <- summary(portGambleANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  portGambleANOVATukeyResults[(i-2), 3] <- summary(portGambleANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
} #Add ANOVA F statistics and p-values to the table
head(portGambleANOVATukeyResults) #Confirm that tests were completed
portGambleANOVATukeyResults$ANOVA.adjusted.BH.pvalue <- p.adjust(p = portGambleANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values using B-H method and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(portGambleANOVATukeyResults) #Confirm addition
#write.csv(portGambleANOVATukeyResults, "2018-02-15-OneWayANOVA-TukeyHSD-by-Habitat-PortGamble-pValues.csv") #Wrote out table for future analyses

#Skokomish River Delta
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  skokomishRiverANOVA <- aov(skokomishRiverData[,i] ~ skokomishRiverData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between eelgrass conditions
  skokomishRiverANOVATukeyResults[(i-2), 2] <- summary(skokomishRiverANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  skokomishRiverANOVATukeyResults[(i-2), 3] <- summary(skokomishRiverANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
} #Add ANOVA F statistics and p-values to the table
head(skokomishRiverANOVATukeyResults) #Confirm that tests were completed
skokomishRiverANOVATukeyResults$ANOVA.adjusted.BH.pvalue <- p.adjust(p = skokomishRiverANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values using B-H method and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(skokomishRiverANOVATukeyResults) #Confirm addition
#write.csv(skokomishRiverANOVATukeyResults, "2018-02-15-OneWayANOVA-TukeyHSD-by-Habitat-SkokomishRiver-pValues.csv") #Wrote out table for future analyses

#Willapa Bay
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  willapaBayANOVA <- aov(willapaBayData[,i] ~ willapaBayData$Eelgrass.Condition) #Perform an ANOVA to test for significant differences between eelgrass conditions
  willapaBayANOVATukeyResults[(i-2), 2] <- summary(willapaBayANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  willapaBayANOVATukeyResults[(i-2), 3] <- summary(willapaBayANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
} #Add ANOVA F statistics and p-values to the table
head(willapaBayANOVATukeyResults) #Confirm that tests were completed
willapaBayANOVATukeyResults$ANOVA.adjusted.BH.pvalue <- p.adjust(p = willapaBayANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values using B-H method and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(willapaBayANOVATukeyResults) #Confirm addition
#write.csv(willapaBayANOVATukeyResults, "2018-02-15-OneWayANOVA-TukeyHSD-by-Habitat-WillapaBay-pValues.csv") #Wrote out table for future analyses