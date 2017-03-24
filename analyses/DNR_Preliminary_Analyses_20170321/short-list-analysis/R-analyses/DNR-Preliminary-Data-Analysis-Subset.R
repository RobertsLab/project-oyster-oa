#In this script, I'll visualize my stress subset data to see what differences there are between eelgrass and bare patches across sites.

#Step 1: Import data
stressProteins <- read.table("/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/short-list-analysis/nmds-heatmap-short-list.txt", header = TRUE, na.strings = "NA")
stressProteins[is.na(stressProteins)] <- 0 #Replace NAs with 0s
stressProteins #Confirm changes

#Step 2: Nonmetric multidimensional scaling plot (NMDS)
#Load the source file for the biostats package
source("/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/biostats.R")
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- stressProteins[-1]
rownames(area.protID2) <- stressProteins[,1]

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
area.t <- t(area.protID2[,1:10])
area.tra <- (area.t+1)
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area.tra, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5)), cex = 1.5)
legend(x = -0.10, y = 0.05, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), bty = "n")

###################################################

#Step 3: Heatmap

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area.tra dataset.

#Create heatmap
pheatmap(area.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = T)

#Export preliminary heatmap as a .png
png(filename = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/short-list-analysis/R-analyses/subsetHeatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(area.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = T)
dev.off()