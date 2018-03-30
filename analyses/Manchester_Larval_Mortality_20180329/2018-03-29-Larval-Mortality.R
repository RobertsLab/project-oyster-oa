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

bucketNames <- paste("Bucket", bucketNumbers, sep = "") #Create a new vector of bucket names
nBuckets <- length(larvalCounts.trans) - 1 #Save number of buckets as a new value

larvalCountsBuckets <- NULL #Create an empty dataframe to store information

for(i in 1:nBuckets){
  tempData <- data.frame("Date" = larvalCounts.trans$Date,
                               "Count" = larvalCounts.trans[,i],
                               "Bucket" = rep(bucketNames[i], times = length(larvalCounts.trans$Date))) #Create an individual data frame for counts from a specific bucket
  larvalCountsBuckets <- rbind(larvalCountsBuckets, tempData) #Save that in the new combined dataframe
} #Create and populate a new dataframe
head(larvalCountsBuckets) #Confirm changes

attach(larvalCountsBuckets)
larvalCountsBuckets <- larvalCountsBuckets[order(Bucket),] #Reorder so buckets are sorted alphabetically
head(larvalCountsBuckets) #Confirm sorting
tail(larvalCountsBuckets) #Confirm sorting
detach(larvalCountsBuckets)

larvalCountsBuckets$Colors <- c(rep("red", times = 50),
                                rep("pink", times = 50),
                                rep("green", times = 50),
                                rep("blue", times = 50),
                                rep("orange", times = 40)) #Create a column to assign colors
#Buckets 1-5 are FLxML: red
#Buckets 6-10 are FLxMA: pink
#Buckets 11-15 are FAxML: green
#Buckets 16-20 are FAxMA: blue
#Buckets 21-24 are HSxHS: orange
head(larvalCountsBuckets) #Confirm column addition

#### PLOT DATA ####

plot(x = larvalCountsBuckets$Date, y = larvalCountsBuckets$Count, type = "l")


?plot
