#In this script, I'll visualize my preliminary data to see what differences there are between eelgrass and bare patches across sites.

#Step 1: Import data
averageAreaAdjustedMerged <- read.csv(file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Skyline_20170314/Oyster-AverageAdjustedMergedArea.csv", header = TRUE)
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Removing the extra column "X"

#Step 2: Nonmetric multidimensional scaling plot (NMDS)
#Load the source file for the biostats package
source("biostats.R")
require(vegan) #Require the vegan package

#make sure first column of protein names is recognized as row names instead of values
nsaf.protID2<-nsaf.dat[-1]
rownames(nsaf.protID2)<-nsaf.dat[,1]

#transpose the file so that rows and columns are switched and normalized by log(x+1)
nsaf.t<-t(nsaf.protID2)
nsaf.tra<-(nsaf.t+1)
nsaf.tra<-data.trans(nsaf.tra, method='log', plot=F)

#make MDS dissimilarity matrix
proc.nmds<-metaMDS(nsaf.tra, distance='bray', k=2, trymax=100, autotransform=F)

#make figure
fig.nmds<-ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
points(fig.nmds, 'sites', col=c(whatever colors you want), pch=(whatever points you want))
legend(x=, y=, pch=(points you chose above), legend=c('legend text 1', 'legend text 2', ...), col=c(colors you chose above))