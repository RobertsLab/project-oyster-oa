#In this script, I'll regress peptide abundance against different biomarkers to see if the biomarker variables explain variation in abundance.

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT GROWTH DATA ####

growthData <- read.csv("../../data/DNR/2017-12-05-OysterGrowth.csv", header = TRUE) #Import growth data
head(growthData) #Confirm import
growthData <- growthData[, -c(1, 3, 4, 5, 8:12)] #Remove extra columns: Round, Habitat, Exclosure, Rep, Ishell1-5
head(growthData) #Confirm column deletion

#### CALCUALTE MEAN SIZE AT BEGINNING OF OUTPLANT ####

mean(growthData$AvgIshell, na.rm = TRUE) #Mean = 27.19098

#### CALCULATE PERCENT GROWTH ####

growthData$diffLength <- growthData$Fshell - growthData$AvgIshell #Find the difference in shell length
growthData$percentGrowth <- ((growthData$diffLength)/(growthData$AvgIshell))*100 #Calculate percent growth
head(growthData) #Confirm calculations

#### INITIAL BOXPLOT ####

boxplot(growthData$percentGrowth ~ growthData$Site, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot using all growth information, not just those subsetted for proteomic analyses
title(ylab = "Percent Growth", line = 2.3, cex.lab = 2) #Add the y-axis label
stripchart(growthData$percentGrowth ~ growthData$Site, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
siteANOVA <- aov(growthData$percentGrowth ~ growthData$Site) #Perform an ANOVA to test for significant differences between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA

TukeyHSD(siteANOVA) #Signficant differences between FB-CI (0.0012929), PG-FB (0.0388860), SK-FB (0.0000009), SK-PG (0.0223310), WB-SK (0.0017054)

#### SUBSET SAMPLES USED FOR PROTEOMIC ANALYSES ####

biologicalReplicates <- read.csv("../../data/DNR/2017-11-21-Alex-Data-Yaamini-Samples-Only.csv") #Import biomarker dataset with sample numbers
head(biologicalReplicates) #Confirm import
biologicalReplicates <- biologicalReplicates[, c(4, 6)] #Want only PRVial and Site information
head(biologicalReplicates) #Confirm column deletions

growthSamplesOnly <- merge(x = growthData, y = biologicalReplicates, by = "PRVial") #Merge datasets by PRVial information
colnames(growthSamplesOnly) #Get column names
growthSamplesOnly <- growthSamplesOnly[, -7] #Remove Site.x column
head(growthSamplesOnly) #Confirm changes

#### SAMPLES ONLY BOXPLOT ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-19-Growth-Data-Analyses/2017-12-19-Growth-Differences-Between-Sites.jpeg", height = 1000, width = 2000)
boxplot(growthSamplesOnly$percentGrowth ~ growthSamplesOnly$Site, main = "Growth Differences between Sites", cex.main = 3, cex.axis = 1.5) #Create the boxplot using all growth information, not just those subsetted for proteomic analyses
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Percent Growth", cex.lab = 2.5, line = 2.2) #Add y-axis label
stripchart(growthSamplesOnly$percentGrowth ~ growthSamplesOnly$Site, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
siteANOVA <- aov(growthSamplesOnly$percentGrowth ~ growthSamplesOnly$Site) #Perform an ANOVA to test for significant differences between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
#dev.off()

TukeyHSD(siteANOVA) #Significant pairwise differences: FB-CI (0.0313908) and SK-FB (0.0102010).

#### IMPORT AND FORMAT PEPTIDE ABUNDANCE DATA ####

SRMDataNMDSAveragedCorrected <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv", header = TRUE) #Import modified dataset. This dataset has rownames as the first column, column names as sample IDs. Area data is averaged and corrected (no NAs)
rownames(SRMDataNMDSAveragedCorrected) <- SRMDataNMDSAveragedCorrected[,1]
SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveragedCorrected[,-1] #Remove first column of rownames
head(SRMDataNMDSAveragedCorrected) #Confirm import.

SRMDataNMDSAveragedCorrectedTransposed <- data.frame(t(SRMDataNMDSAveragedCorrected)) #Transpose the data
SRMDataNMDSAveragedCorrectedTransposed$PRVial <- rownames(SRMDataNMDSAveragedCorrectedTransposed) #Save rownames as a new column
head(SRMDataNMDSAveragedCorrectedTransposed) #Confirm changes

#### MERGE DATASETS ####

peptideGrowthData <- merge(x = SRMDataNMDSAveragedCorrectedTransposed, y = growthSamplesOnly, by = "PRVial") #Merge dataframes by PRVial
head(peptideGrowthData) #Confirm changes. Peptides are columns 2-38, percentGrowth is 43.

#### ASSIGN COLORS BY SITE ####

attach(peptideGrowthData) #Attach dataframe
peptideGrowthData <- peptideGrowthData[order(Site),] #Reorder so sites are sorted alphabetically
head(peptideGrowthData) #Confirm sorting
tail(peptideGrowthData) #Confirm sorting
detach(peptideGrowthData)
peptideGrowthData$Colors <- c(rep(x = "red", times = sum(peptideGrowthData$Site == "CI")),
                                 rep(x = "blue", times = sum(peptideGrowthData$Site == "FB")),
                                 rep(x = "magenta", times = sum(peptideGrowthData$Site == "PG")),
                                 rep(x = "green", times = sum(peptideGrowthData$Site == "SK")),
                                 rep(x = "black", times = sum(peptideGrowthData$Site == "WB"))) #Create a color vector
head(peptideGrowthData) #Confirm addition
tail(peptideGrowthData) #Confirm addition

#### CHANGE WORKING DIRECTORY ####

setwd("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-19-Growth-Data-Analyses/2017-12-19-Peptide-Growth-Scatterplots/") #Change working directory
getwd() #Confirm changes

#### CREATE SCATTERPLOTS ####

nPeptides <- 38 #Columns 2-38 are peptides
for(i in 2:nPeptides) { #For all peptides
  peptideGrowthModel <- lm(peptideGrowthData[,i] ~ peptideGrowthData$percentGrowth, na.action = na.omit)
  fileName <- paste(colnames(peptideGrowthData)[i], "vs.", "Percent Growth", ".jpeg")
  jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
  plot(x = peptideGrowthData$percentGrowth, y = peptideGrowthData[,i], xlab = "Percent Growth", ylab = "Abundance", type = "n", cex.lab = 1.5, cex.axis = 1.5, main = paste(colnames(peptideGrowthData)[i], "vs.", "Percent Growth"), cex.main = 1.75) #Create plot, but do not plot points
  text(x = peptideGrowthData$percentGrowth, y = peptideGrowthData[,i], labels = peptideGrowthData$PRVial, cex = 2, col = peptideGrowthData$Colors, font = 2) #Plot sample ID instead of points
  abline(peptideGrowthModel) #Plot regression
  legend("topleft", bty = "n", legend = paste("R2 =", format(summary(peptideGrowthModel)$adj.r.squared, digits=4))) #Plot R-squared value
  dev.off() #Turn off plotting device
}

#Because none of the R-squared values were much different from zero, I will not be saving any of the information in a separate table.