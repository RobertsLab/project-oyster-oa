#In this script, I'll depict normalized protein area across samples as bar charts after filtering transitions with a coefficient of variation > 15 and technical replicates with ordination distances > 0.20.

#### SET WORKING DIRECTORY ####
setwd("../../../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

SRMDataNMDSAveragedCorrected <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-Averaged-Pivoted-Corrected-SRM-Data-after-CV15-and-Distance-Filtering.csv", header = TRUE) #Import modified dataset. This dataset has rownames as the first column, column names as sample IDs. Area data is averaged and corrected (no NAs)
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
                               modifier = rep("Eelgrass", length(boxplotData))) #Make filename sheet
boxplotFilenames$siteFilenames <- paste(boxplotFilenames$protein, ".jpeg") #Make a new column for the site only filenames
boxplotFilenames$eelgrassFilenames <-paste(boxplotFilenames$protein, boxplotFilenames$modifier, ".jpeg") #Mae a new column for the site and eelgrass filenames
head(boxplotFilenames) #Confirm changes

#### CHANGE WORKING DIRECTORY ####

setwd("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-26-CV-15/2017-10-26-CV-15-Boxplots/")
getwd()

#### MAKE BOXPLOTS JUST BASED ON SITES ####

nTransitions <- (length(boxplotData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nTransitions) { #For all of my columns with transition IDs
  fileName <- boxplotFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "Abundance") #Create the boxplot
  title(fileName)
  dev.off() #Close file
}

#### MAKE BOXPLOTS BASED ON SITES AND EELGRASS CONDITIONS ####

for(i in 3:nTransitions) { #For all of my columns with transition IDs
  fileName <- boxplotFilenames$eelgrassFilenames[i] #Set the file name choices as the third column
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site + boxplotData$Eelgrass.Condition, xlab = "Sites", ylab = "Abundance", names = c("CI.Bare", "FB.Bare", "PG.Bare", "SK.Bare", "WB.Bare", "CI.Eelgrass", "FB.Eelgrass", "PG.Eelgrass", "SK.Eelgrass", "WB.Eelgrass")) #Create the boxplot
  title(fileName)
  dev.off() #Close file
}
