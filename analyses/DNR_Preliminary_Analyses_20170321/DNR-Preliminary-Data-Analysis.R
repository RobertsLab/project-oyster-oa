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
points(fig.nmds, 'sites', col=c(1:5), pch=c(15, 15, 15, 15, 15, 20, 20, 20, 20, 20))
legend(x= , y=, pch=c(15, 15, 15, 15, 15, 20, 20, 20, 20, 20), legend=c('bareCaseInlet', 'bareFidalgoBay', 'bareWillapaBay', 'bareSkokomishRiver', 'barePortGamble', 'eelgrassCaseInlet', 'eelgrassFidalgoBay', 'eelgrassWillapaBay', 'eelgrassSkokomishRiver', 'eelgrassPortGamble'), col=c(1:5))

#Step 3: Heatmap

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area.tra dataset.

#Create heatmap
pheatmap(area.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "preliminaryHeatmap.png")
pheatmap(area.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)
dev.off()

#Modify heatmap parameters

#Step 4: Bubble Plots
#ggplot bubble plots
library(ggplot2)
library(ggthemes)

#I had to use geom_jitter because my points were too close together, you may want to use geom_point
#width and height are specific to jitter
#alpha makes the bubbles transparent
#geom_hline puts horizontal lines between y-axis values
ggplot(data) + geom_jitter(aes(x=column name with x values, y=column name with y values, colour=column to dictate point color, size=column to dictate bubble sizes), alpha=0.5, width=0.2, height=0.2) + labs(y='y-axis label') + geom_hline(yintercept=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5), color='grey80') + theme_tufte()