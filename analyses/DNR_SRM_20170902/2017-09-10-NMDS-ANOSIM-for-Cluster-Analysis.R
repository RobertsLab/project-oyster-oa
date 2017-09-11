#Before going through this script, I went through "2017-09-06-NMDS-for-Technical-Replication"

#### AVERAGE TECHNICAL REPLICATES ####

#After examining how my technical replicates are clustering together, I will average and proceed with an ANOSIM and NMDS

head(SRMDataNMDSPivotedCorrected) #Dataset I'll use to average technical replicates, from my first R script (NMDS for Technical Replication)
SRMDataNMDSAveraged <- data.frame(x = rep(x = 0, times = 111),
                                  y = rep(x = 0, times = 111)) #Create an empty dataframe to store averaged values
row.names(SRMDataNMDSAveraged) <- SRMDataNMDSPivotedCorrected$RowNames #Add row names
head(SRMDataNMDSAveraged) #Confirm changes
for(i in 1:89) { #Average normalized area values for consecutive columns
  SRMDataNMDSAveraged[,i] <- (SRMDataNMDSPivotedCorrected[,i]+SRMDataNMDSPivotedCorrected[,i+1])/2
}
head(SRMDataNMDSAveraged) #Confirm averaging
SRMDataNMDSAveraged <- SRMDataNMDSAveraged[seq(from = 1, to = 89, by = 2)] #Remove even-numbered columns, since those consecutive columns are not technical replicates
head(SRMDataNMDSAveraged) #Confirm column removal
colnames(SRMDataNMDSAveraged) <- technicalReplicates[seq(from = 1, to = 89, by = 2)] #Add column names
head(SRMDataNMDSAveraged) #Confirm column naming

#### NMDS FOR SITE AND EELGRASS CLUSTERING ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS. It can be found the project-oyster-oa repo >> analyses >> DNR_Preliminary_Analyses_20170321. It is also at the bottom of this script.
install.packages("vegan") #Install vegan package
library(vegan)

SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveraged #Duplicate dataframe
SRMDataNMDSAveragedCorrected[is.na(SRMDataNMDSAveragedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSAveragedCorrected) #Confirm there are no NAs

area.protID2 <- SRMDataNMDSAveragedCorrected[,-45] #Save all area data as a new dataframe except for OBLNK2
head(area.protID2) #Confirm changes

area2.t <- t(area.protID2) #Transpose the file so that rows and columns are switched
head(area2.t) #Confirm transposition
area2.tra <- (area2.t+1) #Add 1 to all values before transforming
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.euclidean <- metaMDS(area2.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.euclidean) #Make Shepard plot
vec.proc.nmds.euclidean <- envfit(proc.nmds.euclidean$points, area2.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
plot(vec.proc.nmds.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### ASSIGN COLORS AND SHAPES ####

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(Sample.Number = technicalReplicates,
                            y = rep(x = 0, times = length(technicalReplicates))) #Create a temporary dataframe with technical replicate names used in NMDS
head(temporaryData) #Confirm dataframe creation
NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "Sample.Number") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[-c(89:90),-2] #Remove OBLNK2 and empty column
tail(NMDSColorShapeCustomization) #Confirm removal
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[seq(from = 1, to = 87, by = 2),] #Keep only every other row
head(NMDSColorShapeCustomization)
tail(NMDSColorShapeCustomization)

#Create a color and shape palette
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Site),] #Reorder so sites are sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Colors <- c(rep(x = "red", times = sum(NMDSColorShapeCustomization$Site == "CI")),
            rep(x = "blue", times = sum(NMDSColorShapeCustomization$Site == "FB")),
            rep(x = "black", times = sum(NMDSColorShapeCustomization$Site == "PG")),
            rep(x = "green", times = sum(NMDSColorShapeCustomization$Site == "SK")),
            rep(x = "magenta", times = sum(NMDSColorShapeCustomization$Site == "WB"))) #Create a color vector
NMDSColorShapeCustomization[,4] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Eelgrass.Condition),] #Reorder so eelgrass condition is sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Shapes <- c(rep(x = 16, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Bare")),
                 rep(x = 17, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Eelgrass"))) #Make a shape vector
NMDSColorShapeCustomization[,5] <- NMDS.Shapes #Add the shape vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Sample.Number),] #Resort by sample number
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
colnames(NMDSColorShapeCustomization) <- c("Sample.Number", "Site", "Eelgrass.Condition", "Color", "Shape") #Change column names
head(NMDSColorShapeCustomization) #Confirm change

#### NMDS REFINEMENT ####

fig.nmds <- ordiplot(proc.nmds.euclidean, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("bottomleft", pch = c(rep(16, 17), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

#jpeg(filename = "2017-09-08-NMDS-TechnicalReplication-Normalized.jpeg", width = 1000, height = 1000)
#ordiplot(proc.nmds.euclidean, choices = c(1,2), type = "text", display = "sites") #Plot refined NMDS displaying only samples with their names
#dev.off()

#### NOTES FROM JULIAN ####

#shy away from cluster --> testing among groups instead (are the groups sig diff from eachother)
#ANOSIM: ranked dissimilarities (+NMDS) (can do two-way analysis)
#permanova + PCoA = PCA with euclidean
#different sensitivity to sample size