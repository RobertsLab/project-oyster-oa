#normalized spectral abundance factor (NSAF) calculation
#Your input file should have a column of protein names, a column of protein lengths, and columns of spectral counts (spc) for each sample
#read in your file
#divide all individual spc by protein length (L)
spcL<-filename[,1st column of:last column of spc]/filename[,column with protein lengths]

#sum spc/L for each replicate
sum.spcL<-apply(spcL, 2, sum)

#divide each spc/L by the sum of all spc/L for that replicate
spcL.w.sum<-rbind(spcL, sum.spcL)
spcL.sum.mat<-as.matrix(spcL.w.sum)
nsaf<-t(t(spcL.sum.mat)/spcL.sum.mat[row containing sum of spc/L,])

#rejoin with protein IDs
nsaf.dat<-cbind(filename[,column of protein names], nsaf[1:last row of data,])

################################################################################################################
#Nonmetric multidimensional scaling plot (NMDS)
library(vegan)
#load sourcefile library for biostats

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

###########################################################################################################
#heat map
library(pheatmap)
#read in your file
#your data should be log(x) or log(x+1) transformed
pheatmap(data name, cluster_rows=T, cluster_cols=T, clustering_distance_rows='euclidean', clustering_distance_cols='euclidean', clustering_method='average', show_rownames=F)
#check out the command and play around with the parameters to get the heatmap you want

#ggplot bubble plots
library(ggplot2)
library(ggthemes)

#I had to use geom_jitter because my points were too close together, you may want to use geom_point
#width and height are specific to jitter
#alpha makes the bubbles transparent
#geom_hline puts horizontal lines between y-axis values
ggplot(data) + geom_jitter(aes(x=column name with x values, y=column name with y values, colour=column to dictate point color, size=column to dictate bubble sizes), alpha=0.5, width=0.2, height=0.2) + labs(y='y-axis label') + geom_hline(yintercept=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5), color='grey80') + theme_tufte()