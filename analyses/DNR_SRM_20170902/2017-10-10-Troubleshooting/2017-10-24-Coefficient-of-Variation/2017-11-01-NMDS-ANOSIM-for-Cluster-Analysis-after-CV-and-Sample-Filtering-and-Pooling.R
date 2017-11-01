#Steven went through my normalized data set of technical replicates and calculated a coefficient of variation for each transition. He removed any transitions that had CVs greater than 20. I then went through and removed any samples missing more than half of its transition data. Using this dataset, I'm going to pool my bare and eelgrass samples and remake my NMDS.

#### SET WORKING DIRECTORY ####

setwd("../..") #Change working directory to the master SRM folder, and not the folder where this script is hosted.
getwd() #Confirm changes

#### IMPORT DATA ####

SRMAveragedAreasPivoted <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-30-Protein-Areas-for-Boxplots.csv", header = TRUE) #Import modified dataset
head(SRMAveragedAreasPivoted) #Confirm import.
rownames(SRMAveragedAreasPivoted) <- SRMAveragedAreasPivoted[,1] #Set first column as rownames
SRMAveragedAreasPivoted <- SRMAveragedAreasPivoted[,-1] #Remove first column
head(SRMAveragedAreasPivoted) #Confirm changes

#### CORRECT DATA FOR NMDS ####

SRMAveragedAreasPivotedCorrected <- SRMAveragedAreasPivoted #Duplicate dataframe
SRMAveragedAreasPivotedCorrected[is.na(SRMAveragedAreasPivotedCorrected)] <- 0 #Replace NAs with 0s
head(SRMAveragedAreasPivotedCorrected) #Confirm there are no NAs

#### MAKE BASE NMDS PLOT ####

source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

area.protID <- SRMAveragedAreasPivotedCorrected #Save all area data as a new dataframe
head(area.protID) #Confirm change

area.t <- t(area.protID) #Transpose the file so that rows and columns are switched
head(area.t) #Confirm transposition
area.tra <- (area.t+1) #Add 1 to all values before transforming
area.tra <- data.trans(area.tra, method = 'log', plot = FALSE) #log(x+1) transformation

proc.nmds.norm.averaged.euclidean <- metaMDS(area.t, distance = 'euclidean', k = 2, trymax = 10000, autotransform = FALSE) #Make MDS dissimilarity matrix using euclidean distance. Julian confirmed that I should use euclidean distances, and not bray-curtis
#stressplot(proc.nmds.norm.averaged.euclidean) #Make Shepard plot
#vec.proc.nmds.norm.averaged.euclidean <- envfit(proc.nmds.norm.averaged.euclidean$points, area.t, perm = 1000) #Calculate loadings
ordiplot(proc.nmds.norm.averaged.euclidean, choices = c(1,2), type = "points", display = "sites") #Plot basic NMDS
#plot(vec.proc.nmds.norm.averaged.euclidean, p.max=.01, col='blue') #Plot eigenvectors

#### ASSIGN COLORS AND SHAPES ####

sampleIDs <- c("O01", "O04", "O06", "O08", "O10", "O100", "O101", "O102", "O103", "O106", "O118", "O121", "O122", "O124", "O128", "O131", "O137", "O14", "O140", "O145", "O147", "O17", "O21", "O22", "O24", "O26", "O30", "O31", "O32", "O35", "O40", "O43", "O46", "O49", "O51", "O52", "O56", "O60", "O64", "O66", "O71", "O78", "O90", "O91", "O96", "O99") #Create a sample ID vector

biologicalReplicates <- read.csv("2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-10-25-Biological-Replicate-Information-SampleID-Only.csv", header = TRUE) #Import biological replicate information
head(biologicalReplicates) #Confirm import

#Create a dataframe with biological replicate information for samples used in NMDS
temporaryData <- data.frame(sampleIDs = sampleIDs,
                            y = rep(x = 0, times = length(sampleIDs))) #Create a temporary dataframe with sample IDs used in NMDS
head(temporaryData) #Confirm dataframe creation

NMDSColorCustomization <- merge(x = temporaryData, y = biologicalReplicates, by = "sampleIDs") #Merge biological information with samples used
head(NMDSColorCustomization) #Confirm merge
tail(NMDSColorCustomization) #Confirm merge
NMDSColorCustomization <- NMDSColorCustomization[, -2] #Remove empty column
tail(NMDSColorCustomization) #Confirm removal
NMDSColorCustomization$sampleIDs #Confirm all sample IDs are there
colnames(NMDSColorCustomization) <- c("sample", "site", "eelgrassCondition") #Rename first sample column
head(NMDSColorCustomization) #Confirm changes

#Create a color and shape palette
attach(NMDSColorCustomization)
NMDSColorCustomization <- NMDSColorCustomization[order(site),] #Reorder so sites are sorted alphabetically
head(NMDSColorCustomization) #Confirm sorting
detach(NMDSColorCustomization)
NMDS.Colors <- c(rep(x = "red", times = sum(NMDSColorCustomization$site == "CI")),
                 rep(x = "blue", times = sum(NMDSColorCustomization$site == "FB")),
                 rep(x = "magenta", times = sum(NMDSColorCustomization$site == "PG")),
                 rep(x = "green", times = sum(NMDSColorCustomization$site == "SK")),
                 rep(x = "black", times = sum(NMDSColorCustomization$site == "WB"))) #Create a color vector
NMDSColorCustomization[,3] <- NMDS.Colors #Add the color vector to the dataframe
head(NMDSColorCustomization) #Confirm addition
attach(NMDSColorCustomization)
NMDSColorCustomization <- NMDSColorCustomization[order(sample),] #Reorder so samples are sorted alphabetically
head(NMDSColorCustomization) #Confirm sorting
detach(NMDSColorCustomization)
colnames(NMDSColorCustomization) <- c("Sample.Number", "Site", "Color") #Change column names
head(NMDSColorCustomization) #Confirm change

#### NMDS REFINEMENT ####

#jpeg(filename = "2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation/2017-11-01-NMDS-Norm-Analysis-Averaged-Pooled.jpeg", width = 1000, height = 1000)
fig.nmds <- ordiplot(proc.nmds.norm.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab= "Axis 1", ylab= "Axis 2", cex = 0.5) #Save NMDS as a new object

#Legend for NMDS plot:
#Case Inlet = Red
#Fidalgo Bay = Blue
#Willapa Bay = Black
#Skokomish River Delta = Green
#Port Gamble Bay = Magenta

points(fig.nmds, "sites", col = NMDSColorCustomization$Color, pch = 16)
legend("topright", cex = .5, pch = c(rep(x = 16, times = 5)), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble"), col=c('red', 'blue', 'black', 'green', 'magenta'))
title("Protein Expression Similarities between Sites and Habitats")
#dev.off()

#### ANOSIM ####

dissimArea.t <- vegdist(area.t, "euclidean") #Calculate dissimilarity matrix
ANOSIMReplicates <- biologicalReplicates #Duplicate dataframe
row.names(ANOSIMReplicates) <- ANOSIMReplicates[,1] #Assign sample numbers as row names
ANOSIMReplicates <- ANOSIMReplicates[,-1] #Remove Sample.Number column
head(ANOSIMReplicates) #Confirm changes

ANOSIMReplicates$site <- factor(ANOSIMReplicates$site) #Make sure residual factors are no longer present
str(ANOSIMReplicates) #Confirm new factors

siteNormANOSIM <- anosim(dat = dissimArea.t, grouping = ANOSIMReplicates[,1]) #One-way ANOSIM by Site presence
summary(siteNormANOSIM)
plot(siteNormANOSIM)
#simper(proc.nmds.norm.averaged.euclidean, ANOSIMReplicates$site) #Error in rowSums(comm, na.rm = TRUE): 'x' must be an array of at least two dimensions