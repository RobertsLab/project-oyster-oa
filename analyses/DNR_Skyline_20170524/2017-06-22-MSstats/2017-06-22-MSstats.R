##### BARE VS. EELGRASS #####

#### INSTALL MSSTATS ####
source("http://bioconductor.org/biocLite.R")
biocLite("MSstats")
library(MSstats)

#### IMPORT AND PROCESS DATA ####

rawPeakAreas <- read.csv("2017-06-22-skyline-to-msstats-peak-areas.csv", na.strings = "#N/A") #Import peak area data. This does not include oyster 2 (the first O107 file before I had to make a new sample and redo it).
names(rawPeakAreas) <- c("ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge", "PrecursorMz", "FragmentIon", "ProductCharge", "ProductMz", "IsotopeLabelType", "Condition", "BioReplicate", "FileName", "Area", "StandardType", "Truncated", "DetectionQValue") #Rename columns
head(rawPeakAreas)
peakAreas <- SkylinetoMSstatsFormat(rawPeakAreas) #Convert raw Skyline format to MSstats format. Went from 6964020 obs of 16 variables to 189100 obs of 10 variables.
head(peakAreas)
QuantData <- dataProcess(peakAreas)
head(QuantData$ProcessedData) #View processed data

#### CREATE A CONTRAST MATRIX ####

levels(QuantData$ProcessedData$GROUP_ORIGINAL) #Levels: "Bare" and "Eelgrass"
comparison <- matrix(c(-1, 1), nrow = 1)
row.names(comparison) <- "Eelgrass-Bare"

#### GROUP COMPARISON ####

testResultOneComparison <- groupComparison(contrast.matrix = comparison, data = QuantData) #1: In optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp),: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q. 2: In optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp),: convergence code 3 from bobyqa: bobyqa -- a trust region step failed to reduce q
testResultOneComparison$ComparisonResult #View results
proteinComparisonResults <- testResultOneComparison$ComparisonResult #Save as new dataframe
write.csv(proteinComparisonResults, "2017-06-23-MSstats-BarevEelgrass-Differential-Expression.csv", col.names = c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue"), row.names = F) #Write out data frame

#### GROUP COMPARISON PLOTS ####

#Can't do a heatmap, because at least two comparisons are needed. There were no significantly different proteins expressed between bare and eelgrass treatments.

groupComparisonPlots(data = proteinComparisonResults, type = "VolcanoPlot") #Volcano plot, alpha = 0.05.
groupComparisonPlots(data = proteinComparisonResults, type = "ComparisonPlot") #Comparison plot, alpha = 0.05

#I will repeat this process by site, to see if that will give me any differentially expressed proteins.

##### BETWEEN SITE VARIATION #####

#### IMPORT AND PROCESS DATA ####

rawPeakAreas <- read.csv("2017-06-30-skyline-to-msstats-peak-areas-sites-only.csv", na.strings = "#N/A") #Import peak area data.
names(rawPeakAreas) <- c("ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge", "PrecursorMz", "FragmentIon", "ProductCharge", "ProductMz", "IsotopeLabelType", "Condition", "BioReplicate", "FileName", "Area", "StandardType", "Truncated", "DetectionQValue") #Rename columns
head(rawPeakAreas)
peakAreas <- SkylinetoMSstatsFormat(rawPeakAreas) #Convert raw Skyline format to MSstats format. Went from 6964020 obs of 16 variables to 189100 obs of 10 variables.
head(peakAreas)
QuantData <- dataProcess(peakAreas)
head(QuantData$ProcessedData) #View processed data

#### CREATE A CONTRAST MATRIX ####

levels(QuantData$ProcessedData$GROUP_ORIGINAL) #Levels: "CI" "FB" "PG" "SK" "WB"
comparison1 <- matrix(c(-1, 1, 0, 0, 0), nrow = 1) #CI-FB
comparison2 <- matrix(c(-1, 0, 1, 0, 0), nrow = 1) #CI-PG
comparison3 <- matrix(c(-1, 0, 0, 1, 0), nrow = 1) #CI-SK
comparison4 <- matrix(c(-1, 0, 0, 0, 1), nrow = 1) #CI-WB
comparison5 <- matrix(c(0, -1, 1, 0, 0), nrow = 1) #FB-PG
comparison6 <- matrix(c(0, -1, 0, 1, 0), nrow = 1) #FB-SK
comparison7 <- matrix(c(0, -1, 0, 0, 1), nrow = 1) #FB-WB
comparison8 <- matrix(c(0, 0, -1, 1, 0), nrow = 1) #PG-SK
comparison9 <- matrix(c(0, 0, -1, 0, 1), nrow = 1) #PG-WB
comparison10 <- matrix(c(0, 0, 0, -1, 1), nrow = 1) #SK-WB
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6, comparison7, comparison8, comparison9, comparison10) #Merge all pairwise comparisons together
row.names(comparison) <- c("CI-FB", "CI-PG", "CI-SK", "CI-WB", "FB-PG", "FB-SK", "FB-WB", "PG-SK", "PG-WB", "SK-WB")

#### GROUP COMPARISON ####

testResultSiteComparison <- groupComparison(contrast.matrix = comparison, data = QuantData) #Compare pairs
testResultSiteComparison$ComparisonResult #View results
proteinComparisonSiteResults <- testResultSiteComparison$ComparisonResult #Save as new dataframe
write.csv(proteinComparisonSiteResults, "2017-06-30-MSstats-Sites-Differential-Expression.csv", col.names = c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue"), row.names = F) #Write out data frame

#### GROUP COMPARISON PLOTS ####

groupComparisonPlots(data = proteinComparisonSiteResults, type = "VolcanoPlot") #Volcano plot, alpha = 0.05.
groupComparisonPlots(data = proteinComparisonSiteResults, type = "ComparisonPlot") #Comparison plot, alpha = 0.05
groupComparisonPlots(data = proteinComparisonSiteResults, type = "Heatmap") #Heatmap

#The last step is to compare across both sites and eelgrass conditions.

##### SITES AND EELGRASS CONDITIONS #####

#### IMPORT AND PROCESS DATA ####
rawPeakAreas <- read.csv("2017-06-30-skyline-to-msstats-peak-areas-sites-eelgrass.csv", na.strings = "#N/A") #Import peak area data.
names(rawPeakAreas) <- c("ProteinName", "PeptideSequence", "PeptideModifiedSequence", "PrecursorCharge", "PrecursorMz", "FragmentIon", "ProductCharge", "ProductMz", "IsotopeLabelType", "Condition", "BioReplicate", "FileName", "Area", "StandardType", "Truncated", "DetectionQValue") #Rename columns
head(rawPeakAreas)
peakAreas <- SkylinetoMSstatsFormat(rawPeakAreas) #Convert raw Skyline format to MSstats format. Went from 6964020 obs of 16 variables to 189100 obs of 10 variables.
head(peakAreas)
QuantData <- dataProcess(peakAreas)
head(QuantData$ProcessedData) #View processed data

#### CREATE A CONTRAST MATRIX ####

levels(QuantData$ProcessedData$GROUP_ORIGINAL) #Levels: "BareCI" "BareFB" "BarePG" "BareSK" "BareWB" "EelgrassCI" "EelgrassFB" "EelgrassPG" "EelgrassSK" "EelgrassWB"
comparison1 <- matrix(c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1) #BareCI-BareFB
comparison2 <- matrix(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 1) #BareCI-BarePG
comparison3 <- matrix(c(-1, 0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 1) #BareCI-BareSK
comparison4 <- matrix(c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0), nrow = 1) #BareCI-BareWB
comparison5 <- matrix(c(-1, 0, 0, 0, 0, 1, 0, 0, 0, 0), nrow = 1) #BareCI-EelgrassCI
comparison6 <- matrix(c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 1) #BareFB-BarePG
comparison7 <- matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 1) #BareFB-BareSK
comparison8 <- matrix(c(0, -1, 0, 0, 1, 0, 0, 0, 0, 0), nrow = 1) #BareFB-BareWB
comparison9 <- matrix(c(0, -1, 0, 0, 0, 0, 1, 0, 0, 0), nrow = 1) #BareFB-EelgrassFB
comparison10 <- matrix(c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0), nrow = 1) #BarePG-BareSK
comparison11 <- matrix(c(0, 0, -1, 0, 1, 0, 0, 0, 0, 0), nrow = 1) #BarePG-BareWB
comparison12 <- matrix(c(0, 0, -1, 0, 0, 0, 0, 1, 0, 0), nrow = 1) #BarePG-EelgrassPG
comparison13 <- matrix(c(0, 0, 0, -1, 1, 0, 0, 0, 0, 0), nrow = 1) #BareSK-BareWB
comparison14 <- matrix(c(0, 0, 0, -1, 0, 0, 0, 0, 1, 0), nrow = 1) #BareSK-EelgrassSK
comparison15 <- matrix(c(0, 0, 0, 0, -1, 0, 0, 0, 0, 1), nrow = 1) #BareWB-EelgrassWB
comparison16 <- matrix(c(0, 0, 0, 0, 0, -1, 1, 0, 0, 0), nrow = 1) #EelgrassCI-EelgrassFB
comparison17 <- matrix(c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0), nrow = 1) #EelgrassCI-EelgrassPG
comparison18 <- matrix(c(0, 0, 0, 0, 0, -1, 0, 0, 1, 0), nrow = 1) #EelgrassCI-EelgrassSK
comparison19 <- matrix(c(0, 0, 0, 0, 0, -1, 0, 0, 0, 1), nrow = 1) #EelgrassCI-EelgrassWB
comparison20 <- matrix(c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0), nrow = 1) #EelgrassFB-EelgrassPG
comparison21 <- matrix(c(0, 0, 0, 0, 0, 0, -1, 0, 1, 0), nrow = 1) #EelgrassFB-EelgrassSK
comparison22 <- matrix(c(0, 0, 0, 0, 0, 0, -1, 0, 0, 1), nrow = 1) #EelgrassFB-EelgrassWB
comparison23 <- matrix(c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0), nrow = 1) #EelgrassPG-EelgrassSK
comparison24 <- matrix(c(0, 0, 0, 0, 0, 0, 0, -1, 0, 1), nrow = 1) #EelgrassPG-EelgrassWB
comparison25 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1), nrow = 1) #EelgrassSK-EelgrassWB
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6, comparison7, comparison8, comparison9, comparison10, comparison11, comparison12, comparison13, comparison14, comparison15, comparison16, comparison17, comparison18, comparison19, comparison20, comparison21, comparison22, comparison23, comparison24, comparison25) #Merge all pairwise comparisons together
row.names(comparison) <- c("BareCI-BareFB", "BareCI-BarePG", "BareCI-BareSK", "BareCI-BareWB", "BareCI-EelgrassCI", "BareFB-BarePG", "BareFB-BareSK", "BareFB-BareWB", "BareFB-EelgrassFB", "BarePG-BareSK", "BarePG-BareWB", "BarePG-EelgrassPG", "BareSK-BareWB", "BareSK-EelgrassSK", "BareWB-EelgrassWB", "EelgrassCI-EelgrassFB", "EelgrassCI-EelgrassPG", "EelgrassCI-EelgrassSK", "EelgrassCI-EelgrassWB", "EelgrassFB-EelgrassPG", "EelgrassFB-EelgrassSK", "EelgrassFB-EelgrassWB", "EelgrassPG-EelgrassSK", "EelgrassPG-EelgrassWB", "EelgrassSK-EelgrassWB")

#### GROUP COMPARISON ####

testResultSiteEelgrassComparison <- groupComparison(contrast.matrix = comparison, data = QuantData) #Compare pairs
testResultSiteEelgrassComparison$ComparisonResult #View results
proteinComparisonSiteEelgrassResults <- testResultSiteEelgrassComparison$ComparisonResult #Save as new dataframe
write.csv(proteinComparisonSiteEelgrassResults, "2017-06-30-MSstats-Sites-Eelgrass-Differential-Expression.csv", col.names = c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue"), row.names = F) #Write out data frame

#### GROUP COMPARISON PLOTS ####

groupComparisonPlots(data = proteinComparisonSiteEelgrassResults, type = "VolcanoPlot") #Volcano plot, alpha = 0.05.
groupComparisonPlots(data = proteinComparisonSiteEelgrassResults, type = "ComparisonPlot") #Comparison plot, alpha = 0.05
groupComparisonPlots(data = proteinComparisonSiteEelgrassResults, type = "Heatmap") #Heatmap