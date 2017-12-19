#In this script, I'll visualize salinity data from Micah between sites. I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

salinityData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Salinity-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA")
head(salinityData)

#### SUBSET DATA ####
#I only want the salinity data.

colnames(salinityData)
salinityData <- salinityData[, -c(1, 10:14)] #Save salinity data as a new dataframe
head(salinityData)
colnames(salinityData) <- c("DateTime", "Date", "Time", "CIB", "FBB", "PGE", "SKB", "WBB") #Rename columns
head(salinityData)

#### CALCULATE RANGE OF SALINITY ####

salinityRange <- range(salinityData[, 4:8], na.rm = TRUE) #Calculate range of salinity values
salinityRange[1] <- 10 #Change minimum value to a smaller number to better visualize fluctuations
salinityRange[2] <- 45 #Change maximum value to a round number
salinityRange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

salinityBoxplotCIB <- data.frame(Date.Time = salinityData$DateTime,
                                 Site = rep(x = "CI", times = length(salinityData$DateTime)),
                                 salinity = salinityData$CIB) #Create a dataframe for CIB data
salinityBoxplotCIB <- salinityBoxplotCIB[-c(7456:7457),] #Remove last two rows

salinityBoxplotFBB <- data.frame(Date.Time = salinityData$DateTime,
                                    Site = rep(x = "FB", times = length(salinityData$DateTime)),
                                    salinity = salinityData$FBB) #Create a dataframe for FBB data
salinityBoxplotFBB <- salinityBoxplotFBB[-c(7489:7490),] #Remove last two rows

salinityBoxplotPGE <- data.frame(Date.Time = salinityData$DateTime,
                                    Site = rep(x = "PG", times = length(salinityData$DateTime)),
                                    salinity = salinityData$PGE) #Create a dataframe for PGE data
salinityBoxplotPGE <- salinityBoxplotPGE[-c(7489:7490),] #Remove last two rows

salinityBoxplotSKB <- data.frame(Date.Time = salinityData$DateTime,
                                    Site = rep(x = "SK", times = length(salinityData$DateTime)),
                                    salinity = salinityData$SKB) #Create a dataframe for SKB data
salinityBoxplotSKB <- salinityBoxplotSKB[-c(7489:7490),] #Remove last two rows

salinityBoxplotWBB <- data.frame(Date.Time = salinityData$DateTime,
                                    Site = rep(x = "WB", times = length(salinityData$DateTime)),
                                    salinity = salinityData$WBB) #Create a dataframe for WBB data
salinityBoxplotWBB <- salinityBoxplotWBB[-c(7489:7490),] #Remove last two rows

salinityBoxplot <- rbind(salinityBoxplotCIB, salinityBoxplotFBB, salinityBoxplotPGE, salinityBoxplotSKB, salinityBoxplotWBB) #Bind together all values.

#### MAKE BOXPLOT JUST BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Salinity-QC-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(salinityBoxplot$salinity ~ salinityBoxplot$Site, ylim = salinityRange, main = "Salinity at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites
siteANOVA <- aov(salinityBoxplot$salinity ~ salinityBoxplot$Site) #Perform an ANOVA to test for significant differences in salinitys between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Salinity", cex.lab = 2.5, line = 2.2) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for salinity differences between sites. All pairwise differences are significant at 0.05 level.

#### VISUALIZE DIURNAL FLUCTUATIONS AND BOXPLOT ####

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Diurnal-Salinity-QC-Fluctuations-and-Boxplot.jpeg", height = 5000, width = 4000)

par(mfrow = c(3,2)) #Create multipanel plot with 3 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(salinityData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, type = "l", cex.main = 10, cex.axis = 5, col = "red", main = "Case Inlet") #Case Inlet
abline(h = median(salinityData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity
mtext(side = 2, text = "Salinity", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(salinityData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay") #Fidalgo Bay
abline(h = median(salinityData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, cex.axis = 5, type = "l", col = "magenta", main = "Port Gamble Bay") #Port Gamble Bay
abline(h = median(salinityData$PGE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$PGE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta") #Skokomish River Delta
abline(h = median(salinityData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, cex.axis = 5, type = "l", col = "dark grey", main = "Willapa Bay") #Willapa Bay
abline(h = median(salinityData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity
axis(side = 1, at = seq(from = 1, to = length(salinityData$Date), by = 144*5), lab = salinityData$Date[seq(from = 1, to = length(salinityData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 7) #Modify x-axis labels

boxplot(salinityBoxplot$salinity ~ salinityBoxplot$Site, xaxt = "n", ylim = salinityRange, yaxt = "n", main = "Salinity at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites
siteANOVA <- aov(salinityBoxplot$salinity ~ salinityBoxplot$Site) #Perform an ANOVA to test for significant differences in salinity between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 7) #Modify x-axis label

dev.off()

#Be sure to clear all plot history to reset par.