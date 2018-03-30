#In this script, I'll plot my larval abundance data. Then, I'll see if analysis is worth it.

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

headerDates <- as.Date(c("07-30-2017", "08-01-2017", "08-03-2017", "08-05-2017", "08-07-2017", "08-09-2017", "08-11-2017", "08-13-2017", "08-15-2017", "08-17-2017"), format = "%m-%d-%Y") #Create string of dates for column names
colnames(larvalCounts) <- headerDates #Rename columns
head(larvalCounts) #Confirm changes.

larvalCounts.trans <- data.frame(t(larvalCounts)) #Transpose dataframe and save it as a new dataframe
bucketNames <- paste("Bucket", bucketNumbers, sep = "") #Create a new vector of bucket names
colnames(larvalCounts.trans) <- bucketNames #Use bucket numbers as column names
larvalCounts.trans$Date <- rownames(larvalCounts.trans) #Save rownames as a new column
head(larvalCounts.trans) #Confirm changes

#### PLOT DATA ####
#I will first plot larval count data from the four pH treatment families. I will not include error bars on these estimates since my goal is to first see how much the data overlap to begin with.

#Buckets 1-5 are FLxML: deeppink
#Buckets 6-10 are FLxMA: purple
#Buckets 11-15 are FAxML: green2
#Buckets 16-20 are FAxMA: royalblue
#Buckets 21-24 are HSxHS: not plotting

countRange <- c(0, 80000) #Define range of plot

#jpeg("analyses/Manchester_Larval_Mortality_20180329/2018-03-30-Larval-Counts-Over-Time.jpeg", height = 1000, width = 1500)
plot(larvalCounts.trans$Bucket1, ylab = "Number Live Larvae", ylim = countRange, xaxt = "n", xlab = "", cex.lab = 3, col = "deeppink", type = "p") #Plot just bucket 1 counts
axis(side = 1, at = seq(from = 1, to = 10, by = 1), lab = larvalCounts.trans$Date[seq(from = 1, to = length(larvalCounts.trans$Date), by = 1)]) #Add x axis
lines(larvalCounts.trans$Bucket2, type = "p", col = "deeppink") #Add Bucket 2
lines(larvalCounts.trans$Bucket3, type = "p", col = "deeppink") #Add Bucket 3
lines(larvalCounts.trans$Bucket4, type = "p", col = "deeppink") #Add Bucket 4
lines(larvalCounts.trans$Bucket5, type = "p", col = "deeppink") #Add Bucket 5

lines(larvalCounts.trans$Bucket6, type = "p", col = "purple") #Add Bucket 6
lines(larvalCounts.trans$Bucket7, type = "p", col = "purple") #Add Bucket 7
lines(larvalCounts.trans$Bucket8, type = "p", col = "purple") #Add Bucket 8
lines(larvalCounts.trans$Bucket9, type = "p", col = "purple") #Add Bucket 9
lines(larvalCounts.trans$Bucket10, type = "p", col = "purple") #Add Bucket 10

lines(larvalCounts.trans$Bucket16, type = "p", col = "royalblue") #Add Bucket 16
lines(larvalCounts.trans$Bucket17, type = "p", col = "royalblue") #Add Bucket 17
lines(larvalCounts.trans$Bucket18, type = "p", col = "royalblue") #Add Bucket 18
lines(larvalCounts.trans$Bucket19, type = "p", col = "royalblue") #Add Bucket 19
lines(larvalCounts.trans$Bucket20, type = "p", col = "royalblue") #Add Bucket 20

legend("topright", cex = 1, pch = rep(1, times = 4), legend = c("FL x ML", "FL x MA", "FA x ML", "FA x MA"), col = c("deeppink", "purple", "green2", "royalblue")) #Add legend

#dev.off() #Turn off plotting device

#Counts from each treatment overlap a bit. I don't think anything has a significantly different count. Rates could be different between certain days, but without slopes, it's hard to really tell. Now I will plot counts from each parental treatment family separately with connecting lines to understand the trends better. I will keep the same color scheme as above.

#FL x ML
#jpeg("analyses/Manchester_Larval_Mortality_20180329/2018-03-30-FL-ML-Larval-Counts-Over-Time.jpeg", height = 1000, width = 1500)
plot(larvalCounts.trans$Bucket1, col = "deeppink", type = "b", ylab = "Number Live Larvae", ylim = countRange, xaxt = "n", xlab = "", cex.lab = 3, main = "Low pH Female x Low pH Male Larvae", cex.main = 5) #Plot just bucket 1 counts
axis(side = 1, at = seq(from = 1, to = 10, by = 1), lab = larvalCounts.trans$Date[seq(from = 1, to = length(larvalCounts.trans$Date), by = 1)]) #Add x axis
lines(larvalCounts.trans$Bucket2, type = "b", col = "deeppink") #Add Bucket 2
lines(larvalCounts.trans$Bucket3, type = "b", col = "deeppink") #Add Bucket 3
lines(larvalCounts.trans$Bucket4, type = "b", col = "deeppink") #Add Bucket 4
lines(larvalCounts.trans$Bucket5, type = "b", col = "deeppink") #Add Bucket 5
#dev.off()

#FL x MA
#jpeg("analyses/Manchester_Larval_Mortality_20180329/2018-03-30-FL-MA-Larval-Counts-Over-Time.jpeg", height = 1000, width = 1500)
plot(larvalCounts.trans$Bucket6, col = "purple", type = "b", ylab = "Number Live Larvae", ylim = countRange, xaxt = "n", xlab = "", cex.lab = 3, main = "Low pH Female x Ambient pH Male Larvae", cex.main = 5) #Plot just bucket 6 counts
axis(side = 1, at = seq(from = 1, to = 10, by = 1), lab = larvalCounts.trans$Date[seq(from = 1, to = length(larvalCounts.trans$Date), by = 1)]) #Add x axis
lines(larvalCounts.trans$Bucket7, type = "b", col = "purple") #Add Bucket 7
lines(larvalCounts.trans$Bucket8, type = "b", col = "purple") #Add Bucket 8
lines(larvalCounts.trans$Bucket9, type = "b", col = "purple") #Add Bucket 9
lines(larvalCounts.trans$Bucket10, type = "b", col = "purple") #Add Bucket 10
#dev.off()

#FA x ML
#jpeg("analyses/Manchester_Larval_Mortality_20180329/2018-03-30-FA-ML-Larval-Counts-Over-Time.jpeg", height = 1000, width = 1500)
plot(larvalCounts.trans$Bucket11, col = "green2", type = "b", ylab = "Number Live Larvae", ylim = countRange, xaxt = "n", xlab = "", cex.lab = 3, main = "Ambient pH Female x Low pH Male Larvae", cex.main = 5) #Plot just bucket 11 counts
axis(side = 1, at = seq(from = 1, to = 10, by = 1), lab = larvalCounts.trans$Date[seq(from = 1, to = length(larvalCounts.trans$Date), by = 1)]) #Add x axis
lines(larvalCounts.trans$Bucket12, type = "b", col = "green2") #Add Bucket 12
lines(larvalCounts.trans$Bucket13, type = "b", col = "green2") #Add Bucket 13
lines(larvalCounts.trans$Bucket14, type = "b", col = "green2") #Add Bucket 14
lines(larvalCounts.trans$Bucket15, type = "b", col = "green2") #Add Bucket 15
#dev.off()

#FA x MA
jpeg("analyses/Manchester_Larval_Mortality_20180329/2018-03-30-FA-MA-Larval-Counts-Over-Time.jpeg", height = 1000, width = 1500)
plot(larvalCounts.trans$Bucket16, col = "royalblue", type = "b", ylab = "Number Live Larvae", ylim = countRange, xaxt = "n", xlab = "", cex.lab = 3, main = "Ambient pH Female x Ambient pH Male Larvae", cex.main = 5) #Plot just bucket 16 counts
axis(side = 1, at = seq(from = 1, to = 10, by = 1), lab = larvalCounts.trans$Date[seq(from = 1, to = length(larvalCounts.trans$Date), by = 1)]) #Add x axis
lines(larvalCounts.trans$Bucket17, type = "b", col = "royalblue") #Add Bucket 17
lines(larvalCounts.trans$Bucket18, type = "b", col = "royalblue") #Add Bucket 18
lines(larvalCounts.trans$Bucket19, type = "b", col = "royalblue") #Add Bucket 19
lines(larvalCounts.trans$Bucket20, type = "b", col = "royalblue") #Add Bucket 20
#dev.off()