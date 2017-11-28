#In this script, I'll see if there are any clutering patterns in my data between site and eelgrass habitats.

#### SET WORKING DIRECTORY ####

setwd("../..") #Set working directory to the main SRM data file
getwd()

#### IMPORT DATA ####

SRMDataNMDSPivotedCorrected <- read.csv("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Technical-Replicates-Pivoted.csv")
rownames(SRMDataNMDSPivotedCorrected) <- SRMDataNMDSPivotedCorrected$X #Set row names
SRMDataNMDSPivotedCorrected <- SRMDataNMDSPivotedCorrected[,-1] #Remove column of row names
head(SRMDataNMDSPivotedCorrected) #Confirm there are no NAs

sampleColumnNames <- c("O01", "O04", "O08", "O10", "O100", "O101", "O102", "O106", "O118", "O121", "O124", "O131", "O137", "O140", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O51", "O56", "O60", "O64", "O66", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector
biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
colnames(biologicalReplicates) <- c("Sample.Number", "Site", "Eelgrass.Condition") #Rename columns
head(biologicalReplicates) #Confirm import

#### AVERAGE TECHNICAL REPLICATES ####

#After examining how my technical replicates are clustering together, I will average and proceed with an ANOSIM and NMDS

SRMDataNMDSAveraged <- data.frame(x = rep(x = 0, times = length(SRMDataNMDSPivotedCorrected$O01.1)),
                                  y = rep(x = 0, times = length(SRMDataNMDSPivotedCorrected$O01.1))) #Create an empty dataframe to store averaged values
row.names(SRMDataNMDSAveraged) <- row.names(SRMDataNMDSPivotedCorrected) #Add row names
head(SRMDataNMDSAveraged) #Confirm changes
tail(SRMDataNMDSAveraged) #Confirm changes
for(i in 1:(length(SRMDataNMDSPivotedCorrected))) { #Average normalized area values for consecutive columns
  SRMDataNMDSAveraged[,i] <- (SRMDataNMDSPivotedCorrected[,i]+SRMDataNMDSPivotedCorrected[,i+1])/2
}
head(SRMDataNMDSAveraged) #Confirm averaging
SRMDataNMDSAveraged <- SRMDataNMDSAveraged[seq(from = 1, to = (length(SRMDataNMDSPivotedCorrected)-1), by = 2)] #Remove even-numbered columns, since those consecutive columns are not technical replicates
head(SRMDataNMDSAveraged) #Confirm column removal
colnames(SRMDataNMDSAveraged) <- sampleColumnNames #Add column names
colnames(SRMDataNMDSAveraged) #Confirm column naming
head(SRMDataNMDSAveraged) #Confirm column naming

#### NMDS FOR SITE AND EELGRASS CLUSTERING ####

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS. It can be found the project-oyster-oa repo >> analyses >> DNR_Preliminary_Analyses_20170321. It is also at the bottom of this script.
install.packages("vegan") #Install vegan package
library(vegan)

SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveraged #Duplicate dataframe
SRMDataNMDSAveragedCorrected[is.na(SRMDataNMDSAveragedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSAveragedCorrected) #Confirm there are no NAs
#write.csv(SRMDataNMDSAveragedCorrected, "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv") #Wrote out dataframe

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
temporaryData <- data.frame(Sample.Number = sampleColumnNames,
                            y = rep(x = 0, times = length(sampleColumnNames))) #Create a temporary dataframe with sample  names
head(temporaryData) #Confirm dataframe creation
NMDSColorShapeCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "Sample.Number") #Merge biological information with samples used
head(NMDSColorShapeCustomization) #Confirm merge
tail(NMDSColorShapeCustomization) #Confirm merge
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[,-2] #Remove empty column
head(NMDSColorShapeCustomization) #Confirm removal

#Add region information (Puget Sound vs. Willapa Bay)
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Site),] #Reorder so sites are sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization$Region <- c(rep("PS", times = (length(NMDSColorShapeCustomization$Site)-6)), rep("WB", times = 6)) #Add regional information
NMDSColorShapeCustomization$NMDS.Region.Shapes <- c(rep(20, times = (length(NMDSColorShapeCustomization$Site)-6)), rep(8, times = 6))
head(NMDSColorShapeCustomization) #Confirm changes
tail(NMDSColorShapeCustomization) #Confirm changes

#Create a color and shape palette
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Site),] #Reorder so sites are sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Colors <- c(rep(x = "red", times = sum(NMDSColorShapeCustomization$Site == "CI")),
            rep(x = "blue", times = sum(NMDSColorShapeCustomization$Site == "FB")),
            rep(x = "magenta", times = sum(NMDSColorShapeCustomization$Site == "PG")),
            rep(x = "green", times = sum(NMDSColorShapeCustomization$Site == "SK")),
            rep(x = "black", times = sum(NMDSColorShapeCustomization$Site == "WB"))) #Create a color vector
NMDSColorShapeCustomization[,6] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Eelgrass.Condition),] #Reorder so eelgrass condition is sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Shapes <- c(rep(x = 16, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Bare")),
                 rep(x = 17, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Eelgrass"))) #Make a shape vector
NMDSColorShapeCustomization[,7] <- NMDS.Shapes #Add the shape vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Sample.Number),] #Resort by sample number
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
colnames(NMDSColorShapeCustomization) <- c("Sample.Number", "Site", "Eelgrass.Condition", "Region", "Region.Shape", "Color", "Shape") #Change column names
head(NMDSColorShapeCustomization) #Confirm changes
tail(NMDSColorShapeCustomization) #Confirm changes

#### NMDS REFINEMENT ####

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-NMDS-Analysis-Averaged.jpeg", width = 1000, height = 750)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
legend("topright", pch = c(rep(x = 16, times = 6), 17), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'), cex = 0.5)
#dev.off()

#### ANOSIM ####

dissimArea4.t <- vegdist(area4.t, "euclidean") #Calculate dissimilarity matrix
ANOSIMReplicates <- biologicalReplicates #Subset sample numbers used as IDs in ANOSIM
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
ANOSIMReplicates$Site.Eelgrass <- paste(ANOSIMReplicates$Site, ANOSIMReplicates$Eelgrass.Condition)
head(ANOSIMReplicates) #Confirm changes

str(ANOSIMReplicates) #Examine structure
ANOSIMReplicates$Site <- factor(ANOSIMReplicates$Site) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Eelgrass.Condition <- factor(ANOSIMReplicates$Eelgrass.Condition) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Site.Eelgrass <- factor(ANOSIMReplicates$Site.Eelgrass) #Make sure only preesnt factors are recognized
str(ANOSIMReplicates) #Confirm structure

siteANOSIM <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIM)
plot(siteANOSIM)

eelgrassANOSIM <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIM)
plot(eelgrassANOSIM)

siteEelgrassANOSIM <- anosim(dat = dissimArea4.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIM)
plot(siteEelgrassANOSIM)

regionANOSIM <- anosim(dat = dissimArea4.t, grouping = NMDSColorShapeCustomization[,4]) #One-way ANOSIM by Region (Puget Sound vs. Willapa Bay)
summary(regionANOSIM)
plot(regionANOSIM)
regionSim <- simper(comm = area4.t, group = NMDSColorShapeCustomization$Region) #Calculate similarity percentages
summary(regionSim) #Show similarity percentages

#average = Average contribution to overall dissimilarity
#sd = Standard deviation of contribution
#ratio = Average to SD ratio
#ava, avb = Average abundances per group
#cumsum = Ordered cumulative contribution

#### NMDS REFINEMENT BY REGION ####
#I'm going to take my averaged normalized data and plot it by region (Puget Sound vs. Willapa Bay).

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-28-NMDS-Analysis-Averaged-by-Region.jpeg", width = 1000, height = 750)
fig.nmds.2 <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object
points(fig.nmds.2, "sites", pch = NMDSColorShapeCustomization$Region.Shape) #Add points
legend("topleft", bty = "n", legend = paste("R = 0.2368", "Significance = 0.031"), cex = 0.8) #Add R and p-value from ANOSIM
legend("topright", pch = c(20, 8), legend=c("Puget Sound", "Willapa Bay"), cex = 0.8)
vec.proc.nmds.averaged.euclidean <- envfit(ord = proc.nmds.averaged.euclidean$points, env = area4.t, perm = 1000, na.rm = TRUE) #Calculate loadings
plot(vec.proc.nmds.averaged.euclidean, p.max = .001, col= "blue", cex = 0.3, lty = 2) #Plot eigenvectors
#dev.off()