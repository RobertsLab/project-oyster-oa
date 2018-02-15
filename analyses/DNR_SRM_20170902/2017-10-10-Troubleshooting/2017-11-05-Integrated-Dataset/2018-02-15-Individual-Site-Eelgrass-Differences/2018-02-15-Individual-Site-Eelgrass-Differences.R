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

#### PERFORM TUKEY HSD POST-HOC TEST ####
#This test can be used to understand where significant ANOVA results come from

siteANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(boxplotData),
                                    "ANOVA.Fstatistic" = rep(x = 0, times = length(boxplotData)),
                                    "ANOVA.pvalue" = rep(x = 0, times = length(boxplotData)),
                                    "FB-CI" = rep(x = 0, times = length(boxplotData)),
                                    "PG-CI" = rep(x = 0, times = length(boxplotData)),
                                    "SK-CI" = rep(x = 0, times = length(boxplotData)),
                                    "WB-CI" = rep(x = 0, times = length(boxplotData)),
                                    "PG-FB" = rep(x = 0, times = length(boxplotData)),
                                    "SK-FB" = rep(x = 0, times = length(boxplotData)),
                                    "WB-FB" = rep(x = 0, times = length(boxplotData)),
                                    "SK-PG" = rep(x = 0, times = length(boxplotData)),
                                    "WB-PG" = rep(x = 0, times = length(boxplotData)),
                                    "WB-SK" = rep(x = 0, times = length(boxplotData))) #Create a dataframe to hold all results
siteANOVATukeyResults <- siteANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(siteANOVATukeyResults) #Confirm changes

#Perform Tukey HSD
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  siteANOVA <- aov(boxplotData[,i] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
  siteANOVATukeyResults[(i-2), 2] <- summary(siteANOVA)[[1]][["F value"]][[1]] #Paste ANOVA F-statistic in table
  siteANOVATukeyResults[(i-2), 3] <- summary(siteANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
  siteTukeyHSD <- TukeyHSD(siteANOVA) #Perform Tukey Honest Significant Difference post-hoc test to determine where ANOVA significance is coming from
  siteANOVATukeyResults[(i-2),4:13] <- siteTukeyHSD$`boxplotData$Site`[,4] #Paste Tukey results into table
} #Add all ANOVA and Tukey HSD p-values to the table
head(siteANOVATukeyResults) #Confirm that tests were completed

#Adjust p-values for multiple comparisons (2-15-2018). After talking to Brent, he mentioned that I should use the Benjamini correction for multiple comparisons by controlling the FDR. I will use a FDR of 10% (0.1).

siteANOVATukeyResults$ANOVA.adjusted.pvalue <- p.adjust(p = siteANOVATukeyResults$ANOVA.pvalue, method = "BH") #Adjust p-values and add a column to the table. I can then compare these p-values to my FDR of 0.1.
head(siteANOVATukeyResults) #Confirm addition
#write.csv(siteANOVATukeyResults, "2017-11-06-OneWayANOVA-TukeyHSD-by-Site-pValues.csv") #Wrote out table for future analyses