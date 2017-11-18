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

#### ASSIGN FILENAMES ####

boxplotFilenames <- data.frame(protein = colnames(boxplotData),
                               modifier = rep(".jpeg", length(boxplotData))) #Make filename sheet
boxplotFilenames$siteFilenames <- paste(boxplotFilenames$protein, boxplotFilenames$modifier) #Make a new column for the site only filenames
head(boxplotFilenames) #Confirm changes

#### CHANGE WORKING DIRECTORY ####

setwd("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-06-Boxplots/")
getwd()

#### MAKE BOXPLOTS JUST BASED ON SITES ####

nPeptides <- (length(boxplotData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- boxplotFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(boxplotData[,i] ~ boxplotData$Site, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(boxplotData[,i] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(boxplotFilenames$protein[i])
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
#write.csv(siteANOVATukeyResults, "2017-11-06-OneWayANOVA-TukeyHSD-by-Site-pValues.csv") #Wrote out table for future analyses

#### POWER ANALYSIS ####

#Install dependencies
install.packages("pwr") #Install the power calculation package
library(pwr) #Load package

#Determine what kind of power I have for a small, medium or large effect sizes. I have k = 5 groups, roughly n = 7 observations per group, and a significance level of 0.05. A small effect size is denoted by f = 0.1, medium is f = 0.25, and large is f = 0.4. These values are suggested by the creator of the package.
pwr.anova.test(k = 5, n = 7, f = 0.1, sig.level = 0.05, power = NULL)[5] #Power = 0.06537487
pwr.anova.test(k = 5, n = 7, f = 0.25, sig.level = 0.05, power = NULL)[5] #Power = 0.163053
pwr.anova.test(k = 5, n = 7, f = 0.4, sig.level = 0.05, power = NULL)[5] #Power = 0.381159

#Determine what kind of effect size I can detect for my experimental design at different power levels. A small effect size is denoted by f = 0.1, medium is f = 0.25, and large is f = 0.4. These values are suggested by the creator of the package.
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 1)[3] #f = 1e+07
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.95)[3] #f = 0.7887625
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.90)[3] #f = 0.7180425
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.85)[3] #f = 0.6699974
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.80)[3] #f = 0.6316528