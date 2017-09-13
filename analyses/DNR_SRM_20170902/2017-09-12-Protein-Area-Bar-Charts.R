#Before going through this script, I went through "2017-09-10-NMDS-ANOSIM-for-Cluster-Analysis." In this script, I'll depict normalized protein area across samples as bar charts

#### DATA MANIPULATION ####
SRMDataNMDSAveragedCorrected #From 2017-09-06-NMDS-for-Technical-Replication. Average normalized area data.
boxplotData <- data.frame(t(SRMDataNMDSAveragedCorrected)) #Transpose the data
boxplotData$Sample.Number <- rownames(boxplotData) #Save rownames as a new column
head(boxplotData) #Confirm changes

biologicalReplicates <- read.csv("2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
rownames(biologicalReplicates) <- biologicalReplicates$Sample.Number #Set sample number as row names
head(biologicalReplicates) #Confirm changes
biologicalReplicates <- biologicalReplicates[-c(50,100),] #Remove blanks
biologicalReplicates$Site <- factor(biologicalReplicates$Site) #Remove 0 as a factor
biologicalReplicates$Eelgrass.Condition <- factor(biologicalReplicates$Eelgrass.Condition) #Remove 0 as a factor
str(biologicalReplicates) #Confirm factor reset

boxplotData <- merge(x = biologicalReplicates, y = boxplotData, by = "Sample.Number") #Merge together
head(boxplotData) #Confirm merge
rownames(boxplotData) <- boxplotData$Sample.Number #Set sample number as row names
boxplotData <- boxplotData[-1] #Remove Sample.Number column
head(boxplotData) #Confirm changes

#### MAKE BOXPLOTS JUST BASED ON SITES ####

nTransitions <- (length(boxplotData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
boxplotFilenames <- data.frame(filenames = colnames(boxplotData),
                               modifier = rep("Eelgrass", 113)) #Make filename sheet
for(i in 3:nTransitions) { #For all of my columns with transition IDs
  fileName <- boxplotFilenames$filenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "Abundance") #Create the boxplot
  dev.off() #Close file
}

#### MAKE BOXPLOTS BASED ON SITES AND EELGRASS CONDITIONS ####

boxplotFilenames$Eelgrass <- paste(boxplotFilenames$filenames, boxplotFilenames$modifier) #Create unique names for eelgrass plots
for(i in 3:nTransitions) { #For all of my columns with transition IDs
  fileName <- boxplotFilenames$Eelgrass[i] #Set the file name choices as the third column
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site + boxplotData$Eelgrass.Condition, xlab = "Sites", ylab = "Abundance", names = c("CI.Bare", "FB.Bare", "PG.Bare", "SK.Bare", "WB.Bare", "CI.Eelgrass", "FB.Eelgrass", "PG.Eelgrass", "SK.Eelgrass", "WB.Eelgrass")) #Create the boxplot
  dev.off() #Close file
}

 #Create the boxplot
?boxplot










