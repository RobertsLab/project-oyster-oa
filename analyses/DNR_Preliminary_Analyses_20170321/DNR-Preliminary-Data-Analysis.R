#In this script, I'll visualize my preliminary data to see what differences there are between eelgrass and bare patches across sites.

#Step 1: Import data
averageAreaAdjustedMerged <- read.csv(file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Skyline_20170314/Oyster-AverageAdjustedMergedArea.csv", header = TRUE, na.strings = "NA")
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Removing the extra column "X"
averageAreaAdjustedMerged[is.na(averageAreaAdjustedMerged)] <- 0 #Replace NAs with 0s
averageAreaAdjustedMerged #Confirm changes

#Step 2: Nonmetric multidimensional scaling plot (NMDS)
#Load the source file for the biostats package
source("biostats.R")
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- averageAreaAdjustedMerged[-1]
rownames(area.protID2) <- averageAreaAdjustedMerged[,1]

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
area.t <- t(area.protID2)
area.tra <- (area.t+1)
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area.tra, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
points(fig.nmds, 'sites', col=c(1), pch=(19))
legend(x=, y=, pch=(19), legend=c('legend text 1', 'legend text 2', ...), col=c(1))