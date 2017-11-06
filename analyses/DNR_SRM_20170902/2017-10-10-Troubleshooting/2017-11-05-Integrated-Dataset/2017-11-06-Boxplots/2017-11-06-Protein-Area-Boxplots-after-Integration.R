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
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "Abundance") #Create the boxplot
  stripchart(boxplotData[,i] ~ boxplotData$Site, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(boxplotData[,i] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("ANOVA p-value =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits=4))) #Plot p-value from ANOVA
  title(fileName)
  dev.off() #Close file
}

#### PERFORM TUKEY HSD POST-HOC TEST ####
#This test can be used to understand where significant ANOVA results come from

temp <- aov(boxplotData[,3] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
summary(temp)[[1]][["Pr(>F)"]][[1]] #Plot p-value from ANOVA
temp2 <- TukeyHSD(temp) #Perform Tukey Honest Significant Difference post-hoc test to determine where ANOVA significance is coming from
temp2$`boxplotData$Site`[,4]
