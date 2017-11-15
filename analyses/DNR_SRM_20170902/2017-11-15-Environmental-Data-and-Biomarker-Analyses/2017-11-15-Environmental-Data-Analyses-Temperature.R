#In this script, I'll visualize temperature data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

environmentalData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA")
head(environmentalData)
environmentalData <- environmentalData[,-c(40:51)] #Get rid of empty columns at end of the dataset
head(environmentalData)

#### SUBSET DATA ####
#I only want the temperature data from the dissolved oxygen loggers.

colnames(environmentalData)
temperatureData <- environmentalData[,c(1,30:39)] #Save temperature data as a new dataframe
head(temperatureData)
colnames(temperatureData) <- c("DateTime", "WBE", "WBB", "SKE", "SKB", "PGE", "PGB", "CIE", "CIB", "FBE", "FBB") #Rename columns
head(temperatureData)

#### VISUALIZE DIURNAL FLUCTUATIONS ####

tempRange <- range(temperatureData$WBE, temperatureData$WBB, temperatureData$SKE, temperatureData$SKB, temperatureData$PGE, temperatureData$PGB, temperatureData$CIE, temperatureData$CIB, temperatureData$FBE, temperatureData$FBB, na.rm = TRUE) #Calculate range of temperature values
tempRange[1] <- 10 #Change minimum value to a round number
tempRange #Confirm changes

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-15-Diurnal-Temperature-Fluctuations.jpeg", height = 5000, width = 2000)
par(mfrow = c(5,2))
plot(x = temperatureData$DateTime, y = temperatureData$FBB, xlab = "", xaxt = "n", ylab = "Temperature (ºC)", ylim = tempRange, type = "l", col = "blue", main = "Fidalgo Bay Bare") #Set up plot with no x axis labels, but with y-axis that encompasses maximum and minimum values
plot(x = temperatureData$DateTime, y = temperatureData$FBE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", type = "o", col = "blue", main = "Fidalgo Bay Eelgrass") #Fidgalo Bay, Eelgrass
plot(x = temperatureData$DateTime, y = temperatureData$PGB, xlab = "", xaxt = "n", ylab = "Temperature (ºC)", ylim = tempRange, type = "l", col = "magenta", main = "Port Gamble Bay Bare") #Port Gamble Bay, Bare
plot(x = temperatureData$DateTime, y = temperatureData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", type = "l", col = "magenta", main = "Port Gamble Bay Eelgrass") #Port Gamble Bay, Eelgrass
plot(x = temperatureData$DateTime, y = temperatureData$SKB, xlab = "", xaxt = "n", ylab = "Temperature (ºC)", ylim = tempRange, type = "l", col = "green", main = "Skokomish River Delta Bare") #Skokomish River Delta, Bare
plot(x = temperatureData$DateTime, y = temperatureData$SKE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n",type = "l", col = "green", main = "Skokomish River Delta Eelgrass") #Skokomish River Delta, Eelgrass
plot(x = temperatureData$DateTime, y = temperatureData$CIB, xlab = "", xaxt = "n", ylab = "Temperature (ºC)", ylim = tempRange, type = "l", col = "red", main = "Case Inlet Bare") #Case Inlet, Bare
plot(x = temperatureData$DateTime, y = temperatureData$CIE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n",type = "l", col = "red", main = "Case Inlet Eelgrass") #Case Inlet, Eelgrass
plot(x = temperatureData$DateTime, y = temperatureData$WBB, xlab = "Date and Time", ylab = "Temperature (ºC)", ylim = tempRange, type = "l", col = "black", main = "Willapa Bay Bare") #Willapa Bay, Bare
plot(x = temperatureData$DateTime, y = temperatureData$WBE, xlab = "Date and Time", ylab = "", ylim = tempRange, yaxt = "n",type = "l", col = "black", main = "Willapa Bay Eelgrass") #Willapa Bay, Eelgrass
dev.off()

#### ASSIGN FILENAMES ####

boxplotFilenames <- data.frame(protein = colnames(boxplotData),
                               modifier = rep(".jpeg", length(boxplotData))) #Make filename sheet
boxplotFilenames$siteFilenames <- paste(boxplotFilenames$protein, boxplotFilenames$modifier) #Make a new column for the site only filenames
head(boxplotFilenames) #Confirm changes

#### CHANGE WORKING DIRECTORY ####

setwd("2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-06-Boxplots/")
getwd()

#### MAKE BOXPLOTS JUST BASED ON SITES ####

nPeptides <- (length(boxplotData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  fileName <- boxplotFilenames$siteFilenames[i] #Set the file name choices as the first column
  jpeg(filename = fileName, width = 1000, height = 750) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "", cex.lab = 2, cex.axis = 1.5) #Create the boxplot
  title(ylab = "Abundance", line = 2.3, cex.lab = 2) #Add the y-axis label
  stripchart(boxplotData[,i] ~ boxplotData$Site, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue') #Add each data point
  siteANOVA <- aov(boxplotData[,i] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
  legend("topleft", bty = "n", legend = paste("ANOVA p-value =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits=4))) #Plot p-value from ANOVA
  title(boxplotFilenames$protein[i])
  dev.off() #Close file
}

#### PERFORM TUKEY HSD POST-HOC TEST ####
#This test can be used to understand where significant ANOVA results come from

siteANOVATukeyResults <- data.frame("Protein.Peptide" = colnames(boxplotData),
                                    "ANOVA.pvalue" = rep(x = 0, times = length(boxplotData)),
                                    "FB-CI" = rep(x = 0, times = length(boxplotData)),
                                    "PG-CI" = rep(x = 0, times = length(boxplotData)),
                                    "SK-CI" = rep(x = 0, times = length(boxplotData)),
                                    "WB-CI" = rep(x = 0, times = length(boxplotData)),
                                    "PG-FB" = rep(x = 0, times = length(boxplotData)),
                                    "SK-FB" = rep(x = 0, times = length(boxplotData)),
                                    "WB-FB" = rep(x = 0, times = length(boxplotData)),
                                    "SK-PG" = rep(x = 0, times = length(boxplotData)),
                                    "WB-PG" = rep(x = 0, times = length(boxplotData)),
                                    "WB-SK" = rep(x = 0, times = length(boxplotData))) #Create a dataframe to hold all results
siteANOVATukeyResults <- siteANOVATukeyResults[-c(1:2),] #Remove the first two rows, since they are not peptides
head(siteANOVATukeyResults) #Confirm changes

#Perform Tukey HSD
for(i in 3:nPeptides) { #For all of my columns with peptide IDs
  siteANOVA <- aov(boxplotData[,i] ~ boxplotData$Site) #Perform an ANOVA to test for significant differences between sites
  siteANOVATukeyResults[(i-2), 2] <- summary(siteANOVA)[[1]][["Pr(>F)"]][[1]] #Paste ANOVA p-value in table
  siteTukeyHSD <- TukeyHSD(siteANOVA) #Perform Tukey Honest Significant Difference post-hoc test to determine where ANOVA significance is coming from
  siteANOVATukeyResults[(i-2),3:12] <- siteTukeyHSD$`boxplotData$Site`[,4] #Paste Tukey results into table
} #Add all ANOVA and Tukey HSD p-values to the table
head(siteANOVATukeyResults) #Confirm that tests were completed
#write.csv(siteANOVATukeyResults, "2017-11-06-OneWayANOVA-TukeyHSD-by-Site-pValues.csv") #Wrote out table for future analyses

#### POWER ANALYSIS ####

#Install dependencies
install.packages("pwr") #Install the power calculation package
library(pwr) #Load package

#Determine what kind of power I have for a small, medium or large effect sizes. I have k = 5 groups, roughly n = 7 observations per group, and a significance level of 0.05. A small effect size is denoted by f = 0.1, medium is f = 0.25, and large is f = 0.4. These values are suggested by the creator of the package.
pwr.anova.test(k = 5, n = 7, f = 0.1, sig.level = 0.05, power = NULL)[5] #Power = 0.06537487
pwr.anova.test(k = 5, n = 7, f = 0.25, sig.level = 0.05, power = NULL)[5] #Power = 0.163053
pwr.anova.test(k = 5, n = 7, f = 0.4, sig.level = 0.05, power = NULL)[5] #Power = 0.381159

#Determine what kind of effect size I can detect for my experimental design at different power levels. A small effect size is denoted by f = 0.1, medium is f = 0.25, and large is f = 0.4. These values are suggested by the creator of the package.
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 1)[3] #f = 1e+07
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.95)[3] #f = 0.7887625
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.90)[3] #f = 0.7180425
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.85)[3] #f = 0.6699974
pwr.anova.test(k = 5, n = 7, f = NULL, sig.level = 0.05, power = 0.80)[3] #f = 0.6316528