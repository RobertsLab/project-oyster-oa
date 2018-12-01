#In this script, I'll see if there are any clutering patterns in my data between site and eelgrass habitats.

#### LOAD DEPENDENCIES ####

source("analyses/DNR_SRM_20170902/biostats.R") #Load the source file for the biostats commands
#install.packages("vegan") #Install vegan package
require(vegan)

#### IMPORT DATA ####

SRMDataNMDSPivotedCorrected <- read.csv("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Technical-Replicates-Pivoted.csv")
rownames(SRMDataNMDSPivotedCorrected) <- SRMDataNMDSPivotedCorrected$X #Set row names
SRMDataNMDSPivotedCorrected <- SRMDataNMDSPivotedCorrected[,-1] #Remove column of row names
head(SRMDataNMDSPivotedCorrected) #Confirm there are no NAs

sampleColumnNames <- c("O01", "O04", "O08", "O10", "O100", "O101", "O102", "O106", "O118", "O121", "O124", "O131", "O137", "O140", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O51", "O56", "O60", "O64", "O66", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector
biologicalReplicates <- read.csv("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
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

#### TRANSFORM DATA ####

SRMDataNMDSAveragedCorrected <- SRMDataNMDSAveraged #Duplicate dataframe
#write.csv(SRMDataNMDSAveragedCorrected, "analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-11-27-Averaged-Areas-Pivoted-NAs.csv") #Wrote out dataframe with NAs
SRMDataNMDSAveragedCorrected[is.na(SRMDataNMDSAveragedCorrected)] <- 0 #Replace NAs with 0s
head(SRMDataNMDSAveragedCorrected) #Confirm there are no NAs
#write.csv(SRMDataNMDSAveragedCorrected, "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-05-Averaged-Areas-Pivoted-Corrected.csv") #Wrote out dataframe

area.protID4 <- SRMDataNMDSAveragedCorrected #Save all area data as a new dataframe
head(area.protID4) #Confirm changes
area4.t <- t(area.protID4) #Transpose the file so that rows and columns are switched
head(area4.t) #Confirm transposition

area4.tra <- data.trans(area4.t, method = 'hellingers', plot = FALSE) #Hellinger (asymmetric) transformation
head(area4.tra) #Confirm transformation

#### NMDS FOR SITE AND EELGRASS CLUSTERING ####

nmds.scree(area4.tra, distance = "euclidean", k = 10, autotransform = FALSE, trymax = 20) #Create a screeplot to compare the stress for solutions across different k values from 2 to 10. Use 20 different random start configurations. As the number of ordination axes increases, stress is minimized because the NMDS algorithm is trying to represent p dimensional data in k dimensions. Using 2 axes is appropriate.

proc.nmds.averaged.euclidean <- metaMDS(area4.tra, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix on hellinger transformed data using euclidean distance.
proc.nmds.averaged.euclidean$stress #Stress of NMDS is 0.07508988

nmds.monte(area4.tra, distance = "euclidean", k = 2, autotransform = FALSE, trymax = 20) #Perform a randomization test to determine if the solution for k dimensions is significant. The observed stress value, 0.07508988, is less than the expected stress value. P-value = 0.00990099
stressplot(proc.nmds.averaged.euclidean) #Make Shepard plot to visualize the relationship between original dissimilarities (distance matrix) and distnaces in ordination space. The non-metric R-squared value is 0.994 (redundant with observed stress value and p-value from the randomization test)

vec.proc.nmds.averaged.euclidean <- envfit(proc.nmds.averaged.euclidean$points, area4.tra, perm = 1000) #Calculate loadings by correlating NMDS scores with original variables
vec.proc.nmds.averaged.euclidean #Look at loadings

ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "text", display = "sites", xlab = "Axis 1", ylab = "Axis 2") #Plot basic NMDS
plot(vec.proc.nmds.averaged.euclidean, p.max = 0.001, col = 'blue') #Plot loadings that are significant at the 0.001 level

#### IMPORT AND FORMAT BIOLOGICAL DATA ####

biologicalReplicates <- read.csv("analyses/DNR_SRM_20170902/2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A", fileEncoding="UTF-8-BOM") #Import site and eelgrass condition information (i.e. biological replicate information), using specific file encoding information
head(biologicalReplicates) #Confirm import
biologicalReplicates$Sample.Number <- as.character(biologicalReplicates$Sample.Number) #Convert sample number to character string
biologicalReplicates$Sample.Number <- substr(biologicalReplicates$Sample.Number, 1, nchar(biologicalReplicates$Sample.Number)-2) #Remove -1 or -2 from end of sample number
biologicalReplicates <- biologicalReplicates[1:50,] #Keep only the first 50 rows, since everything repeats
head(biologicalReplicates) #Confirm changes
tail(biologicalReplicates) #Confirm changes

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
NMDS.Colors <- c(rep(x = "#00A9BD", times = sum(NMDSColorShapeCustomization$Site == "CI")),
            rep(x = "#38001C", times = sum(NMDSColorShapeCustomization$Site == "FB")),
            rep(x = "#440D82", times = sum(NMDSColorShapeCustomization$Site == "PG")),
            rep(x = "#017A74", times = sum(NMDSColorShapeCustomization$Site == "SK")),
            rep(x = "#EB8B0C", times = sum(NMDSColorShapeCustomization$Site == "WB"))) #Create a color vector
NMDSColorShapeCustomization[,6] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Eelgrass.Condition),] #Reorder so eelgrass condition is sorted alphabetically
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
NMDS.Shapes <- c(rep(x = 1, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Bare")),
                 rep(x = 16, times = sum(NMDSColorShapeCustomization$Eelgrass.Condition == "Eelgrass"))) #Make a shape vector
NMDSColorShapeCustomization[,7] <- NMDS.Shapes #Add the shape vector to the dataframe
head(NMDSColorShapeCustomization) #Confirm addition
attach(NMDSColorShapeCustomization)
NMDSColorShapeCustomization <- NMDSColorShapeCustomization[order(Sample.Number),] #Resort by sample number
head(NMDSColorShapeCustomization) #Confirm sorting
detach(NMDSColorShapeCustomization)
colnames(NMDSColorShapeCustomization) <- c("Sample.Number", "Site", "Eelgrass.Condition", "Region", "Region.Shape", "Color", "Shape") #Change column names
head(NMDSColorShapeCustomization) #Confirm changes
tail(NMDSColorShapeCustomization) #Confirm changes

#### NMDS BY SITE WITH CONFIDENCE ELLIPSE ####

#pdf("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-11-28-Protein-Abundance-NMDS-Ellipse.pdf", width = 11, height = 8.5)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object
text(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color) #Add oyster sample IDs to NMDS and color-code by site distinction

ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "CI", col = "#00A9BD") #Add confidence ellipse around the oyster samples from Case Inlet
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "FB", col = "#38001C") #Add confidence ellipse around the oyster samples from Fidalgo Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "PG", col = "#440D82") #Add confidence ellipse around the oyster samples from Port Gamble Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "SK", col = "#017A74") #Add confidence ellipse around the oyster samples from Skokomish River Delta
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "WB", col = "#EB8B0C") #Add confidence ellipse around the oyster samples from Willapa Bay

legend("topright", pch = rep(x = 16, times = 5), legend=c('Case Inlet', "Fidalgo Bay", "Port Gamble Bay", "Skokomish", "Willapa Bay"), col = c('#00A9BD', '#38001C', '#440D82', '#017A74', '#EB8B0C'), cex = 0.5, bty = "n")
#dev.off()

#### NMDS BY SITE WITH POLYGONS ####

#Legend for NMDS plot:
#Case Inlet = Red
#Fidalgo Bay = Blue
#Port Gamble Bay = Magenta
#Skokomish River Delta = Green
#Willapa Bay = Black

ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Create an empty plot
text(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color) #Add oyster sample IDs to NMDS and color-code by site distinction

ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "CI", col = "red") #Add confidence ellipse around the oyster samples from Case Inlet
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "FB", col = "blue") #Add confidence ellipse around the oyster samples from Fidalgo Bay
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "PG", col = "magenta") #Add confidence ellipse around the oyster samples from Port Gamble Bay
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "SK", col = "green") #Add confidence ellipse around the oyster samples from Skokomish River Delta
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "WB", col = "black") #Add confidence ellipse around the oyster samples from Willapa Bay

legend("topright", pch = rep(x = 16, times = 5), legend=c('Case Inlet', "Fidalgo Bay", "Port Gamble Bay", "Skokomish", "Willapa Bay"), col = c('red', 'blue', 'magenta', 'green', 'black'), cex = 0.5)

#### NMDS BY HABITAT WITH CONFIDENCE ELLIPSES ####

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-05-23-NMDS-Analysis-Averaged-HabitatOnly.jpeg", width = 1000, height = 750)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle

points(fig.nmds, "sites", col = "black", pch = NMDSColorShapeCustomization$Shape)
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Eelgrass.Condition, show.groups = "Eelgrass", col = "green") #Add confidence ellipse around the oyster samples from eelgrass habitats
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Eelgrass.Condition, show.groups = "Bare", col = "black") #Add confidence ellipse around the oyster samples from bare habitats
legend("topright", pch = c(16, 17), legend=c("Bare", "Eelgrass"), col=c("black", "green"), cex = 1)
#dev.off()

#### NMDS BY HABITAT WITH POLYGONS ####

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-05-23-NMDS-Analysis-Averaged-HabitatOnly.jpeg", width = 1000, height = 750)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Bare = circle
#Eelgrass = Triangle

points(fig.nmds, "sites", col = "black", pch = NMDSColorShapeCustomization$Shape)
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Eelgrass.Condition, show.groups = "Eelgrass", col = "green") #Add confidence ellipse around the oyster samples from eelgrass habitats
ordihull(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Eelgrass.Condition, show.groups = "Bare", col = "black") #Add confidence ellipse around the oyster samples from bare habitats
legend("topright", pch = c(16, 17), legend=c("Bare", "Eelgrass"), col=c("black", "green"), cex = 1)
#dev.off()

#### NMDS BY SITE AND HABITAT WITH CONFIDENCE ELLIPSE ####

#pdf("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-11-28-Protein-Abundance-Site-Habitat-NMDS-Ellipse.pdf", width = 11, height = 8.5)
fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object
points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape)
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 1, text = "NMDS1", line = 2)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 2, text = "NMDS2", line = 2)
box(col = "grey80")

ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "CI", col = "#00A9BD88") #Add confidence ellipse around the oyster samples from Case Inlet
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "FB", col = "#38001C88") #Add confidence ellipse around the oyster samples from Fidalgo Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "PG", col = "#440D8288") #Add confidence ellipse around the oyster samples from Port Gamble Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "SK", col = "#017A7488") #Add confidence ellipse around the oyster samples from Skokomish River Delta
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "WB", col = "#EB8B0C88") #Add confidence ellipse around the oyster samples from Willapa Bay

legend("topleft", pch = c(rep(x = 1, times = 6), 16), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('#00A9BD', '#38001C', '#440D82', '#017A74', '#EB8B0C', 'black', 'black'), cex = 0.5, bty = "n")

#dev.off()

#### PLOT JUST THE LOADINGS ####

sigLoadings <- envfit(proc.nmds.averaged.euclidean$points, area4.tra[,c(5, 12:13, 19, 22, 30, 33)], perm = 1000, na.rm = TRUE) #Only calculate loadings simper identified as driving differences between FB-WB and SK-WB. See ANOSIM section for simper results.
sigLoadings #View loadings

ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Create an empty plot
plot(sigLoadings, col = 'grey20', labels = c("1", "2", "3", "4", "5", "6", "7")) #Plot loadings that simper determined were significant
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 1, text = "NMDS1", line = 2)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 2, text = "NMDS2", line = 2)
box(col = "grey80")

#### NMDS BY REGION ####
#I'm going to take my averaged normalized data and plot it by region (Puget Sound vs. Willapa Bay).

#jpeg(filename = "2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-28-NMDS-Analysis-Averaged-by-Region.jpeg", width = 1000, height = 750)
fig.nmds.2 <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "Axis 1", ylab = "Axis 2", cex = 0.5) #Save NMDS as a new object
points(fig.nmds.2, "sites", pch = NMDSColorShapeCustomization$Region.Shape) #Add points
#legend("topleft", bty = "n", legend = paste("R = 0.2368", "Significance = 0.031"), cex = 0.8) #Add R and p-value from ANOSIM
legend("topright", pch = c(20, 8), legend=c("Puget Sound", "Willapa Bay"), cex = 0.8)
#vec.proc.nmds.averaged.euclidean <- envfit(ord = proc.nmds.averaged.euclidean$points, env = area4.t, perm = 1000, na.rm = TRUE) #Calculate loadings
#plot(vec.proc.nmds.averaged.euclidean, p.max = 0.001, col= "blue", cex = 0.3, lty = 2) #Plot eigenvectors
#dev.off()

#### GLOBAL ANOSIM ####

ANOSIMReplicates <- biologicalReplicates #Subset sample numbers used as IDs in ANOSIM
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
head(ANOSIMReplicates) #Confirm changes

ANOSIMReplicates <- ANOSIMReplicates[-50,] #Remove OBLNK row
tail(ANOSIMReplicates) #Confirm changes

ANOSIMReplicates$Site.Eelgrass <- paste(ANOSIMReplicates$Site, ANOSIMReplicates$Eelgrass.Condition) #Add a new column with site and eelgrass designation
head(ANOSIMReplicates) #Confirm changes

str(ANOSIMReplicates) #Examine structure of ANOSIMReplicates
ANOSIMReplicates$Site <- factor(ANOSIMReplicates$Site) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Eelgrass.Condition <- factor(ANOSIMReplicates$Eelgrass.Condition) #Make sure only preesnt factors are recognized
ANOSIMReplicates$Site.Eelgrass <- factor(ANOSIMReplicates$Site.Eelgrass) #Make sure only preesnt factors are recognized
str(ANOSIMReplicates) #Confirm structure

dissimArea4.t <- vegdist(area4.tra, "euclidean") #Calculate euclidean dissimilarity matrix

siteANOSIM <- anosim(dissimArea4.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site
summary(siteANOSIM)
siteANOSIM$statistic #R = 0.0645419
siteANOSIM$signif #p = 0.065
plot(siteANOSIM) #Obtain boxplots and permutation test histogram

eelgrassANOSIM <- anosim(dissimArea4.t, grouping = ANOSIMReplicates[,2]) #One-way ANOSIM by Eelgrass presence
summary(eelgrassANOSIM)
eelgrassANOSIM$statistic #R = 0.04427102. Within group and between group similarities are the same for bare and eelgrass habitats
eelgrassANOSIM$signif #p = 0.122. This result is not significant
plot(eelgrassANOSIM)

siteEelgrassANOSIM <- anosim(dissimArea4.t, grouping = ANOSIMReplicates[,3]) #Two-way ANOSIM by Site and Eelgrass
summary(siteEelgrassANOSIM)
siteEelgrassANOSIM$statistic #R = 0.08814604
siteEelgrassANOSIM$signif #p = 0.073
plot(siteEelgrassANOSIM)

regionANOSIM <- anosim(dissimArea4.t, grouping = NMDSColorShapeCustomization[,4]) #One-way ANOSIM by Region (Puget Sound vs. Willapa Bay)
summary(regionANOSIM)
regionANOSIM$statistic #0.2265681. Mild evidence for groupings
regionANOSIM$signif #0.053. Marginally significant
plot(regionANOSIM)

#### PAIRWISE ANOSIM ####

ANOSIMReplicates$Sample.Number <- rownames(ANOSIMReplicates) #Extract row names
attach(ANOSIMReplicates)
ANOSIMReplicates <- ANOSIMReplicates[order(Sample.Number),] #Reorder by sample number
head(ANOSIMReplicates) #Confirm reordering
detach(ANOSIMReplicates)

area4.tra$Sample.Number <- rownames(area4.tra) #Extract row names
head(area4.tra) #Confirm changes

pairwiseData <- merge(x = area4.tra, y = ANOSIMReplicates, by = "Sample.Number") #Merge dataframes
rownames(pairwiseData) <- pairwiseData$Sample.Number #Make sample number the rownames
pairwiseData <- pairwiseData[, -1] #Remove Sample.Number column
head(pairwiseData) #Confirm merge

#CI vs. FB

dataCIFB <- pairwiseData[pairwiseData$Site == "CI" | pairwiseData$Site == "FB", ] #Subset data
dataCIFB$Site <- factor(dataCIFB$Site) #Make sure only present factors are recognized
head(dataCIFB) #Confirm changes

siteCIFBANOSIM <- anosim(dataCIFB[,1:37], grouping = dataCIFB[,38])
summary(siteCIFBANOSIM)
siteCIFBANOSIM$statistic #R = 0.0196793
siteCIFBANOSIM$signif #p = 0.331
plot(siteCIFBANOSIM) #Obtain boxplots and permutation test histogram

#CI vs. PG

dataCIPG <- pairwiseData[pairwiseData$Site == "CI" | pairwiseData$Site == "PG", ] #Subset data
dataCIPG$Site <- factor(dataCIPG$Site) #Make sure only present factors are recognized
head(dataCIPG) #Confirm changes

siteCIPGANOSIM <- anosim(dataCIPG[,1:37], grouping = dataCIPG[,38])
summary(siteCIPGANOSIM)
siteCIPGANOSIM$statistic #R = -0.07069971
siteCIPGANOSIM$signif #p = 0.78
plot(siteCIPGANOSIM) #Obtain boxplots and permutation test histogram

#CI vs. SK

dataCISK <- pairwiseData[pairwiseData$Site == "CI" | pairwiseData$Site == "SK", ] #Subset data
dataCISK$Site <- factor(dataCISK$Site) #Make sure only present factors are recognized
head(dataCISK) #Confirm changes

siteCISKANOSIM <- anosim(dataCISK[,1:37], grouping = dataCISK[,38])
summary(siteCISKANOSIM)
siteCISKANOSIM$statistic #R = -0.05830904
siteCISKANOSIM$signif #p = 0.733
plot(siteCISKANOSIM) #Obtain boxplots and permutation test histogram

#CI vs. WB

dataCIWB <- pairwiseData[pairwiseData$Site == "CI" | pairwiseData$Site == "WB", ] #Subset data
dataCIWB$Site <- factor(dataCIWB$Site) #Make sure only present factors are recognized
head(dataCIWB) #Confirm changes

siteCIWBANOSIM <- anosim(dataCIWB[,1:37], grouping = dataCIWB[,38])
summary(siteCIWBANOSIM)
siteCIWBANOSIM$statistic #R = 0.09259259
siteCIWBANOSIM$signif #p = 0.16
plot(siteCIWBANOSIM) #Obtain boxplots and permutation test histogram

#FB vs. PG

dataFBPG <- pairwiseData[pairwiseData$Site == "FB" | pairwiseData$Site == "PG", ] #Subset data
dataFBPG$Site <- factor(dataFBPG$Site) #Make sure only present factors are recognized
head(dataFBPG) #Confirm changes

siteFBPGANOSIM <- anosim(dataFBPG[,1:37], grouping = dataFBPG[,38])
summary(siteFBPGANOSIM)
siteFBPGANOSIM$statistic #R = -0.03125
siteFBPGANOSIM$signif #p = 0.599
plot(siteFBPGANOSIM) #Obtain boxplots and permutation test histogram

#FB vs. SK

dataFBSK <- pairwiseData[pairwiseData$Site == "FB" | pairwiseData$Site == "SK", ] #Subset data
dataFBSK$Site <- factor(dataFBSK$Site) #Make sure only present factors are recognized
head(dataFBSK) #Confirm changes

siteFBSKANOSIM <- anosim(dataFBSK[,1:37], grouping = dataFBSK[,38])
summary(siteFBSKANOSIM)
siteFBSKANOSIM$statistic #R = 0.09486607
siteFBSKANOSIM$signif #p = 0.086
plot(siteFBPGANOSIM) #Obtain boxplots and permutation test histogram

#FB vs. WB

dataFBWB <- pairwiseData[pairwiseData$Site == "FB" | pairwiseData$Site == "WB", ] #Subset data
dataFBWB$Site <- factor(dataFBWB$Site) #Make sure only present factors are recognized
head(dataFBWB) #Confirm changes

siteFBWBANOSIM <- anosim(dataFBWB[,1:37], grouping = dataFBWB[,38])
summary(siteFBWBANOSIM)
siteFBWBANOSIM$statistic #R = 0.2567829
siteFBWBANOSIM$signif #p = 0.035
plot(siteFBWBANOSIM) #Obtain boxplots and permutation test histogram
simperFBWB <- simper(dataFBWB[,1:37], group = dataFBWB[,38])
summary(simperFBWB)

#PG vs. SK

dataPGSK <- pairwiseData[pairwiseData$Site == "PG" | pairwiseData$Site == "SK", ] #Subset data
dataPGSK$Site <- factor(dataPGSK$Site) #Make sure only present factors are recognized
head(dataPGSK) #Confirm changes

sitePGSKANOSIM <- anosim(dataPGSK[,1:37], grouping = dataPGSK[,38])
summary(sitePGSKANOSIM)
sitePGSKANOSIM$statistic #R = -0.006138393
sitePGSKANOSIM$signif #p = 0.382
plot(sitePGSKANOSIM) #Obtain boxplots and permutation test histogram

#PG vs. WB

dataPGWB <- pairwiseData[pairwiseData$Site == "PG" | pairwiseData$Site == "WB", ] #Subset data
dataPGWB$Site <- factor(dataPGWB$Site) #Make sure only present factors are recognized
head(dataPGWB) #Confirm changes

sitePGWBANOSIM <- anosim(dataPGWB[,1:37], grouping = dataPGWB[,38])
summary(sitePGWBANOSIM)
sitePGWBANOSIM$statistic #R = 0.07267442
sitePGWBANOSIM$signif #p = 0.192
plot(sitePGWBANOSIM) #Obtain boxplots and permutation test histogram

#SK vs. WB

dataSKWB <- pairwiseData[pairwiseData$Site == "SK" | pairwiseData$Site == "WB", ] #Subset data
dataSKWB$Site <- factor(dataSKWB$Site) #Make sure only present factors are recognized
head(dataSKWB) #Confirm changes

siteSKWBANOSIM <- anosim(dataSKWB[,1:37], grouping = dataSKWB[,38])
summary(siteSKWBANOSIM)
siteSKWBANOSIM$statistic #R = 0.1540698
siteSKWBANOSIM$signif #p = 0.079
plot(siteSKWBANOSIM) #Obtain boxplots and permutation test histogram
simperSKWB <- simper(dataSKWB[,1:37], group = dataSKWB[,38])
summary(simperSKWB)