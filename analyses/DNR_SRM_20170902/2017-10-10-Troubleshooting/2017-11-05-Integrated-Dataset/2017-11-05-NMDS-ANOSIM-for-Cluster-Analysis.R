#In this script, I'll see if there are any clutering patterns in my data between site and eelgrass habitats.

#### SET WORKING DIRECTORY ####

setwd("../..") #Set working directory to the main SRM data file
getwd()

#### SOURCE PREVIOUS SCRIPT ####

source("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-for-Technical-Replication.R") #Information from previous script

#### AVERAGE TECHNICAL REPLICATES ####

#After examining how my technical replicates are clustering together, I will average and proceed with an ANOSIM and NMDS

head(SRMDataNMDSPivotedCorrected) #Dataset I'll use to average technical replicates, from my first R script (NMDS for Technical Replication)
SRMDataNMDSAveraged <- data.frame(x = rep(x = 0, times = length(SRMDataNMDSPivotedCorrected$`O01-1`)),
                                  y = rep(x = 0, times = length(SRMDataNMDSPivotedCorrected$`O01-1`))) #Create an empty dataframe to store averaged values
row.names(SRMDataNMDSAveraged) <- row.names(SRMDataNMDSPivotedCorrected) #Add row names
head(SRMDataNMDSAveraged) #Confirm changes
for(i in 1:(length(SRMDataNMDSPivotedCorrected))) { #Average normalized area values for consecutive columns
  SRMDataNMDSAveraged[,i] <- (SRMDataNMDSPivotedCorrected[,i]+SRMDataNMDSPivotedCorrected[,i+1])/2
}
head(SRMDataNMDSAveraged) #Confirm averaging
SRMDataNMDSAveraged <- SRMDataNMDSAveraged[seq(from = 1, to = (length(SRMDataNMDSPivotedCorrected)-1), by = 2)] #Remove even-numbered columns, since those consecutive columns are not technical replicates
head(SRMDataNMDSAveraged) #Confirm column removal
colnames(SRMDataNMDSAveraged) <- technicalReplicates[seq(from = 1, to = (length(SRMDataNMDSPivotedCorrected)-1), by = 2)] #Add column names
colnames(SRMDataNMDSAveraged) #Confirm column naming
head(SRMDataNMDSAveraged) #Confirm column naming

#### NMDS FOR SITE AND EELGRASS CLUSTERING ####

#Load the source file for the biostats package
#source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS. It can be found the project-oyster-oa repo >> analyses >> DNR_Preliminary_Analyses_20170321. It is also at the bottom of this script.
#install.packages("vegan") #Install vegan package
#library(vegan)

SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveraged #Duplicate dataframe
SRMDataNMDSAveragedCorrected[is.na(SRMDataNMDSAveragedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSAveragedCorrected) #Confirm there are no NAs

area.protID4 <- SRMDataNMDSAveragedCorrected #Save all area data as a new dataframe
head(area.protID4) #Confirm changes

area4.t <- t(area.protID4) #Transpose the file so that rows and columns are switched
head(area4.t) #Confirm transposition
area4.tra <- (area4.t+1) #Add 1 to all values before transforming
area4.tra <- data.trans(area4.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.averaged.euclidean <- metaMDS(area4.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
stressplot(proc.nmds.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.averaged.euclidean <- envfit(proc.nmds.averaged.euclidean$points, area4.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### ASSIGN COLORS AND SHAPES ####

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(Sample.Number = technicalReplicates,
                            y = rep(x = 0, times = length(technicalReplicates))) #Create a temporary dataframe with technical replicate names used in NMDS
head(temporaryData) #Confirm dataframe creation
NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "Sample.Number") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[,-2] #Remove empty column
head(NMDSColorShapeCustomization) #Confirm removal
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[seq(from = 1, to = 91, by = 2),] #Keep only every other row
head(NMDSColorShapeCustomization) #Confirm changes
NMDSColorShapeCustomization$Sample.Number #Confirm changes

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

fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("bottomleft", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

#jpeg(filename = "2017-09-11-NMDS-Analysis-Averaged", width = 1000, height = 1000)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

#points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
#legend("bottomleft", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))
#dev.off()

#### ANOSIM ####

dissimArea4.t <- vegdist(area4.t, "euclidean") #Calculate dissimilarity matrix
ANOSIMReplicates <- biologicalReplicates[c(1:49),] #Subset sample numbers used as IDs in ANOSIM
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
head(ANOSIMReplicates) #Confirm changes

str(ANOSIMReplicates) #Examine structure
ANOSIMReplicates$Site <- factor(ANOSIMReplicates$Site) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Eelgrass.Condition <- factor(ANOSIMReplicates$Eelgrass.Condition) #Make sure only preesnt factors are recognized
str(ANOSIMReplicates) #Confirm structure

siteANOSIM <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIM)
plot(siteANOSIM)
simper(proc.nmds.averaged.euclidean, ANOSIMReplicates$Site)

eelgrassANOSIM <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIM)
plot(eelgrassANOSIM)
simper(proc.nmds.averaged.euclidean, ANOSIMReplicates$Eelgrass.Condition)

#Two-way ANOSIM by Site and Eelgrass presence