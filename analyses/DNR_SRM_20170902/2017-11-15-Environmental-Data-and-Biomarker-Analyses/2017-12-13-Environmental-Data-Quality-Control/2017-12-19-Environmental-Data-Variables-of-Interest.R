#In this script, I'll visualize temperature data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

temperatureData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA")
head(temperatureData)

#### SUBSET DATA ####
#I only want the temperature data from the dissolved oxygen loggers.

colnames(temperatureData)
temperatureData <- temperatureData[,c(1:3, seq(from = 33, to = 41, by = 2))] #Save temperature data from bare outplants as a new dataframe
head(temperatureData)
colnames(temperatureData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
head(temperatureData)

#### CALCULATE RANGE OF TEMPERATURES ####

temperatureRange <- range(temperatureData[, 4:8], na.rm = TRUE) #Calculate range of temperature values
temperatureRange[1] <- 10 #Change minimum value to a round number
temperatureRange[2] <- 40 #Change maximum value to a round number
temperatureRange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

temperatureBoxplotCIB <- data.frame(Date.Time = temperatureData$DateTime,
                                 Site = rep(x = "CI", times = length(temperatureData$DateTime)),
                                 temperature = temperatureData$CIB) #Create a dataframe for CIB data
temperatureBoxplotCIB <- temperatureBoxplotCIB[-c(7489:7490),] #Remove last two rows

temperatureBoxplotFBB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "FB", times = length(temperatureData$DateTime)),
                                    temperature = temperatureData$FBB) #Create a dataframe for FBB data
temperatureBoxplotFBB <- temperatureBoxplotFBB[-c(7489:7490),] #Remove last two rows

temperatureBoxplotPGB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "PG", times = length(temperatureData$DateTime)),
                                    temperature = temperatureData$PGB) #Create a dataframe for PGB data
temperatureBoxplotPGB <- temperatureBoxplotPGB[-c(7489:7490),] #Remove last two rows

temperatureBoxplotSKB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "SK", times = length(temperatureData$DateTime)),
                                    temperature = temperatureData$SKB) #Create a dataframe for SKB data
temperatureBoxplotSKB <- temperatureBoxplotSKB[-c(7489:7490),] #Remove last two rows

temperatureBoxplotWBB <- data.frame(Date.Time = temperatureData$DateTime,
                                    Site = rep(x = "WB", times = length(temperatureData$DateTime)),
                                    temperature = temperatureData$WBB) #Create a dataframe for WBB data
temperatureBoxplotWBB <- temperatureBoxplotWBB[-c(7489:7490),] #Remove last two rows

temperatureBoxplot <- rbind(temperatureBoxplotCIB, temperatureBoxplotFBB, temperatureBoxplotPGB, temperatureBoxplotSKB, temperatureBoxplotWBB)

#### MAKE BOXPLOT BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Temperature-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(temperatureBoxplot$temperature ~ temperatureBoxplot$Site, ylim = temperatureRange, main = "Temperature at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(temperatureBoxplot$temperature ~ temperatureBoxplot$Site) #Perform an ANOVA to test for significant differences in temperatures between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Temperature (ºC)", cex.lab = 2.5, line = 2.5) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for temperature differences between sites. All pairwise differences are significant at 0.05 level except for SK-PG (p = 0.9949722).

#### VISUALIZE DIURNAL FLUCTUATIONS AND BOXPLOT ####
#I'm going to make a 3x2 multipanel plot and put the site boxplot in the bottom right corner.

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Temperature-Fluctuations-and-Boxplot.jpeg", height = 5000, width = 4000)

par(mfrow = c(3,2)) #Create multipanel plot with 3 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(temperatureData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = temperatureRange, type = "l", cex.main = 10, cex.axis = 5, col = "red", main = "Case Inlet") #Case Inlet
abline(h = median(temperatureData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature
mtext(side = 2, text = "Temperature (ºC)", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(temperatureData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = temperatureRange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay") #Fidalgo Bay
abline(h = median(temperatureData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$PGB, xlab = "", xaxt = "n", ylab = "", ylim = temperatureRange, cex.main = 10, cex.axis = 5, type = "l", col = "magenta", main = "Port Gamble Bay") #Port Gamble Bay
abline(h = median(temperatureData$PGB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$PGB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = temperatureRange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta") #Skokomish River Delta
abline(h = median(temperatureData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature

plot(temperatureData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = temperatureRange, cex.main = 10, cex.axis = 5, type = "l", col = "dark grey", main = "Willapa Bay") #Willapa Bay
abline(h = median(temperatureData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median temperature
abline(h = mean(temperatureData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean temperature
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 7) #Modify x-axis labels

boxplot(temperatureBoxplot$temperature ~ temperatureBoxplot$Site, xaxt = "n", ylim = temperatureRange, yaxt = "n", main = "Temperature at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites
siteANOVA <- aov(temperatureBoxplot$temperature ~ temperatureBoxplot$Site) #Perform an ANOVA to test for significant differences in temperature between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 7) #Modify x-axis label

#dev.off()

#Be sure to clear all plot history to reset par.