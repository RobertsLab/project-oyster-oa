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

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-15-Diurnal-Temperature-Fluctuations.jpeg", height = 5000, width = 4000)

par(mfrow = c(5,2)) #Create multipanel plot with 5 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(15, 15, 1, 1)) #Remove redundant white space

plot(temperatureData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, cex.axis = 5, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Bare") #FBB. Set up plot with no x axis labels, but with y-axis that encompasses maximum and minimum values
abline(h = median(temperatureData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$FBE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Eelgrass") #Fidgalo Bay, Eelgrass
abline(h = median(temperatureData$FBE, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$FBE, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$PGB, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, cex.axis = 5, cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Bare") #Port Gamble Bay, Bare
abline(h = median(temperatureData$PGB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$PGB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Eelgrass") #Port Gamble Bay, Eelgrass
abline(h = median(temperatureData$PGE, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$PGE, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, cex.axis = 5, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Bare") #Skokomish River Delta, Bare
abline(h = median(temperatureData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$SKE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Eelgrass") #Skokomish River Delta, Eelgrass
abline(h = median(temperatureData$SKE, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$SKE, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, cex.axis = 5, cex.main = 10, type = "l", col = "red", main = "Case Inlet Bare") #Case Inlet, Bare
abline(h = median(temperatureData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$CIE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", type = "l", cex.main = 10, col = "red", main = "Case Inlet Eelgrass") #Case Inlet, Eelgrass
abline(h = median(temperatureData$CIE, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$CIE, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, cex.axis = 5, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Bare") #Willapa Bay, Bare
abline(h = median(temperatureData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature
mtext(side = 2, text = "Temperature (ºC)", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels
axis(side = 1, at = seq(from = 1, to = length(temperatureData$DateTime), by = 144*5), lab = temperatureData$DateTime[seq(from = 1, to = length(temperatureData$DateTime), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis

plot(temperatureData$WBE, xlab = "", xaxt = "n", ylab = "", ylim = tempRange, yaxt = "n", cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Eelgrass") #Willapa Bay, Eelgrass
abline(h = median(temperatureData$WBE, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$WBE, na.rm = TRUE), lty = 2) #Add line depicting mean temperature
axis(side = 1, at = seq(from = 1, to = length(temperatureData$DateTime), by = 144*5), lab = temperatureData$DateTime[seq(from = 1, to = length(temperatureData$DateTime), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date and Time", line = 7, cex = 5, outer = TRUE) #Modify x-axis labels

dev.off()

#### REFORMAT DATA FOR BOXPLOT ####

temperatureBoxplotCIB <- data.frame(Date.Time = temperatureData$DateTime,
                                 Site = rep(x = "CI", times = length(temperatureData$DateTime)),
                                 Habitat = rep(x = "Bare", times = length(temperatureData$DateTime)),
                                 Temperature = temperatureData$CIB) #Create a dataframe for CIB data
temperatureBoxplotCIB <- temperatureBoxplotCIB[-c(7489:7490),] #Remove last two rows
temperatureBoxplotCIE <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "CI", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$CIE) #Create a dataframe for CIE data
temperatureBoxplotCIE <- temperatureBoxplotCIE[-c(7489:7490),] #Remove last two rows

temperatureBoxplotFBB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "FB", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$FBB) #Create a dataframe for FBB data
temperatureBoxplotFBB <- temperatureBoxplotFBB[-c(7489:7490),] #Remove last two rows
temperatureBoxplotFBE <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "FB", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$FBE) #Create a dataframe for FBE data
temperatureBoxplotFBE <- temperatureBoxplotFBE[-c(7489:7490),] #Remove last two rows

temperatureBoxplotPGB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "PG", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$PGB) #Create a dataframe for PGB data
temperatureBoxplotPGB <- temperatureBoxplotPGB[-c(7489:7490),] #Remove last two rows
temperatureBoxplotPGE <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "PG", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$PGE) #Create a dataframe for PGE data
temperatureBoxplotPGE <- temperatureBoxplotPGE[-c(7489:7490),] #Remove last two rows

temperatureBoxplotSKB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "SK", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$SKB) #Create a dataframe for SKB data
temperatureBoxplotSKB <- temperatureBoxplotSKB[-c(7489:7490),] #Remove last two rows
temperatureBoxplotSKE <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "SK", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$SKE) #Create a dataframe for SKE data
temperatureBoxplotSKE <- temperatureBoxplotSKE[-c(7489:7490),] #Remove last two rows

temperatureBoxplotWBB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "WB", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$WBB) #Create a dataframe for WBB data
temperatureBoxplotWBB <- temperatureBoxplotWBB[-c(7489:7490),] #Remove last two rows
temperatureBoxplotWBE <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "WB", times = length(temperatureData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(temperatureData$DateTime)),
                                    Temperature = temperatureData$WBE) #Create a dataframe for WBE data
temperatureBoxplotWBE <- temperatureBoxplotWBE[-c(7489:7490),] #Remove last two rows

temperatureBoxplot <- rbind(temperatureBoxplotCIB, temperatureBoxplotCIE, temperatureBoxplotFBB, temperatureBoxplotFBE, temperatureBoxplotPGB, temperatureBoxplotPGE, temperatureBoxplotSKB, temperatureBoxplotSKE, temperatureBoxplotWBB, temperatureBoxplotWBE)

#### MAKE BOXPLOTS BASED ON SITES AND HABITAT ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-15-Temperature-Boxplot-Site-Habitat.jpeg", height = 1000, width = 2000)
boxplot(temperatureBoxplot$Temperature ~ temperatureBoxplot$Site + temperatureBoxplot$Habitat, ylim = tempRange, names = c("CI.Bare", "FB.Bare", "PG.Bare", "SK.Bare", "WB.Bare", "CI.Eelgrass", "FB.Eelgrass", "PG.Eelgrass", "SK.Eelgrass", "WB.Eelgrass"), col = rep(c("red", "blue", "magenta", "green", "white"), times = 2), main = "Temperature at Site and Habitats", cex.main = 5, cex.axis = 1.5) #Make boxplot based on sites and habitat
title(xlab = "Site and Habitat", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Temperature (ºC)", cex.lab = 2.5, line = 2.2) #Add y-axis label
#dev.off()

#### MAKE BOXPLOT JUST BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-15-Temperature-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(temperatureBoxplot$Temperature ~ temperatureBoxplot$Site, ylim = tempRange, main = "Temperature at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(temperatureBoxplot$Temperature ~ temperatureBoxplot$Site) #Perform an ANOVA to test for significant differences in temperatures between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Temperature (ºC)", cex.lab = 2.5, line = 2.2) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for temperature differences between sites. All pairwise differences are significant at 0.05 level.