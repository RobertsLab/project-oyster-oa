#In this script, I'll plot my larval abundance data with error bars. Then, I'll see if analysis is worth it

#### SET WORKING DIRECTORY ####

getwd()
setwd("../../") #Set working directory as repository

#### IMPORT DATA ####
larvalCounts <- read.csv("data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2018-03-29-Larval-Counts.csv", header = TRUE) #Import larval count data
head(larvalCounts) #Confirm import. The headers are not dates, so I need to fix that.

#### FORMAT DATA ####
bucketNumbers <- larvalCounts$Bucket #Save bucket numbers as a new vector
larvalCounts <- larvalCounts[,-1] #Remove bucket number column. The bucket number is the same as the row number, so I don't need to add it back in.
head(larvalCounts) #Confirm changes

headerDates <- as.Date(c("2017-07-30", "2017-08-01", "2017-08-03", "2017-08-05", "2017-08-07", "2017-08-09", "2017-08-11", "2017-08-13", "2017-08-15", "2017-08-17")) #Create string of dates for column names
colnames(larvalCounts) <- headerDates #Rename columns
head(larvalCounts) #Confirm changes.

larvalCounts.trans <- data.frame(t(larvalCounts)) #Transpose dataframe and save it as a new dataframe
colnames(larvalCounts.trans) <- bucketNumbers #Use bucket numbers as column names
larvalCounts.trans$Date <- rownames(larvalCounts.trans) #Save rownames as a new column
head(larvalCounts.trans) #Confirm changes



#### PLOT DATA ####

