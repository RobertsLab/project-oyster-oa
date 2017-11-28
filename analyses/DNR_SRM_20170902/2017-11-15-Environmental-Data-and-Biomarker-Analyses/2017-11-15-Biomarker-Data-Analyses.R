#In this script, I'll visualize biomarker data from Alex that corresponds to my samples.

#### SET WORKING DIRECTORY ####
#setwd("../../..") #Set working directory to the master SRM folder
#getwd()

#### IMPORT DATA ####

biomarkerData <- read.csv("../../data/DNR/2017-11-21-Alex-Data-Yaamini-Samples-Only.csv") #Import dataset
biomarkerData <- biomarkerData[,-c(1:3, 5, 8:11, 13:14)] #Remove empty columns and columns I don't need (Round, Spp, Bag, Rep.x, FAvial, topValve, bothValves)
head(biomarkerData) #Confirm changes

#FinalHeight = length of shell (mm) at longest point from umbo to edge of shell.
#TissueMass = dry weight (g) of the soft tissue of the oyster after gill was removed
#shellThick = thickness of the shell (mm) at the point we crushed in the shell stength test. measured with modified calipers.
#peakLoad = maximum force sustained by shell (N) measured with an MTS brand force gauge. 
#Strength = maximum force divided by thickness of shell (N/um), this is the main datapoint of interest for the shell strength analysis, and what I observed to be higher in eelgrass for Olympia oysters, and lower in eelgrass for Pacific oysters
#delC = ratio of C13/C12, the standard stable isotope unit for carbon. Often used to determine diet sources, range is typically -12 (enriched) per mil to -24 (depleted) per mil for marine-estuarine species. fresh healthy phytoplankton growing around here have a signature of ~-18 to -22, whereas resuspended detritus and terrestrial detritus is more depleted ~-22 to -26. Generally interpret more enriched signatures as reflecting a better quality diet in estuarine systems.
#C.perc = related to tissue chemistry, % carbon by mass
#N.perc = % Nitrogen by mass, related to protein content. relatively higher %N is usually a good thing.
#CN.ratio = ratio of %carbon/%Nitrogen. often used as an indicator of health where a lower ratio represents high protein content per unit mass.

#### ASSIGN FILENAMES ####

boxplotFilenames <- data.frame(biomarker = c("PRVial", "Site", "Habitat", "Final Shell Height", "Tissue Mass", "Final Shell Thickness", "Peak Load", "Strength", "delC", "Percent C", "delN", "Percent N", "C:N Ratio"),
                               modifier = rep(".jpeg", length(biomarkerData))) #Make filename sheet
boxplotFilenames$siteFilenames <- paste(boxplotFilenames$biomarker, boxplotFilenames$modifier) #Make a new column for the site only filenames
head(boxplotFilenames) #Confirm changes

#### CHANGE WORKING DIRECTORY ####

setwd("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-27-Biomarker-Boxplots")
getwd()

#### VISUALIZE BIOMARKER DATA ####

nBiomarkers <- (length(biomarkerData)) #The number of columns in the dataframe. The first 3 columns are descriptors.
for(i in 4:nBiomarkers) { #For all columns with biomarker data
  fileName <- boxplotFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(biomarkerData[,i] ~ biomarkerData$Site.x, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  stripchart(biomarkerData[,i] ~ biomarkerData$Site.x, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(biomarkerData[,i] ~ biomarkerData$Site.x) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
  title(boxplotFilenames$biomarker[i], cex.main = 3)
  dev.off() #Close file
}