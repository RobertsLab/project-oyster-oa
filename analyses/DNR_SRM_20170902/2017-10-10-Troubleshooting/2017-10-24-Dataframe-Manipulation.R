#In this script, I'll manipulate my dataframe for Steven's use.

#### NONNORMALIZED DATA ####

#### IMPORT DATA ####

nonnormalizedData <- read.csv("2017-09-07-SRM-Data-NMDS-Pivoted.csv", header = TRUE) #Import pivoted nonnormalized data
nonnormalizedData <- nonnormalizedData[,-1] #Remove irrelevant first column
head(nonnormalizedData) #Confirm changes

#### CREATE SAMPLE ID VECTOR ####
sampleIDs <- c("O01", "O04", "O06", "O08", "010", "O100", "O101", "O102", "O103", "O106", "O118", "O121", "O122", "O124", "O128", "O131", "O137", "O14", "O140", "O145", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O49", "O51", "O52", "O56", "O60", "O64", "O66", "O71", "O78", "O90", "O91", "O96", "O99")

#### REFORMAT DATAFRAME ####
technicalReplicatesSampleIDsNonNormalized <- data.frame("Protein.Name" = rep(x = 0, times = 10), "Sample" = rep(x = 0, times = 10), "Replicate1" = rep(x = 0, times = 10), "Replicate2" = rep(x = 0, times = 10))
nSamples <- length(sampleIDs) #Count the number of sample IDs
for(i in 1:nSamples) {
  temp <- data.frame("Protein.Name" = nonnormalizedData[,93], #First column is protein name
                     "Sample" = rep(x = sampleIDs[i], times = length(nonnormalizedData$RowNames)), #Second column is the sample ID
                     "Replicate1" = nonnormalizedData[,((2*i)-1)], #Third column is the area data from the first technical replicate
                     "Replicate2" = nonnormalizedData[,2*i]) #Fourth column is the area data from the second technical replicate
  technicalReplicatesSampleIDsNonNormalized <- rbind(technicalReplicatesSampleIDsNonNormalized, temp)
}
technicalReplicatesSampleIDsNonNormalized <- technicalReplicatesSampleIDsNonNormalized[-c(1:10),] #Remove empty first ten rows
head(technicalReplicatesSampleIDsNonNormalized) #Confirm changes
#write.csv(x = technicalReplicatesSampleIDsNonNormalized, file = "2017-10-24-NonNormalized-Adjacent-Technical-Replicates")

#### NORMALIZED DATA ####

#### IMPORT DATA ####

normalizedData <- read.csv("2017-09-11-SRM-Data-Normalized-NMDS-Pivoted.csv", header = TRUE) #Import pivoted nonnormalized data
normalizedData <- normalizedData[,-1] #Remove irrelevant first column
head(normalizedData) #Confirm changes

#### REFORMAT DATAFRAME ####

technicalReplicatesSampleIDsNormalized <- data.frame("Protein.Name" = rep(x = 0, times = 10), "Sample" = rep(x = 0, times = 10), "Replicate1" = rep(x = 0, times = 10), "Replicate2" = rep(x = 0, times = 10))
nSamples <- length(sampleIDs) #Count the number of sample IDs
for(i in 1:nSamples) {
  temp <- data.frame("Protein.Name" = normalizedData[,93], #First column is protein name
                     "Sample" = rep(x = sampleIDs[i], times = length(normalizedData$RowNames)), #Second column is the sample ID
                     "Replicate1" = normalizedData[,((2*i)-1)], #Third column is the area data from the first technical replicate
                     "Replicate2" = normalizedData[,2*i]) #Fourth column is the area data from the second technical replicate
  technicalReplicatesSampleIDsNormalized <- rbind(technicalReplicatesSampleIDsNormalized, temp)
}
technicalReplicatesSampleIDsNormalized <- technicalReplicatesSampleIDsNormalized[-c(1:10),] #Remove empty first ten rows
head(technicalReplicatesSampleIDsNormalized) #Confirm changes
#write.csv(x = technicalReplicatesSampleIDsNormalized, file = "2017-10-24-Normalized-Adjacent-Technical-Replicates")