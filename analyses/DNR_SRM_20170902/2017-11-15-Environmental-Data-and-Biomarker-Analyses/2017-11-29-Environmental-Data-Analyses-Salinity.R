#In this script, I'll visualize salinity data from Micah between sites. I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

environmentalData <- read.csv("../../data/DNR/2017-11-25-Calculated-Salinity-Output-from-Micah.csv", header = TRUE, na.strings = "NA")
head(environmentalData)

#### SUBSET DATA ####
#I only want the salinity data.

colnames(environmentalData)
salinityData <- environmentalData[,c(1:2, seq(from = 3, to = 20, by = 2))] #Save salinity data as a new dataframe
head(salinityData)
colnames(salinityData) <- c("Date", "Time", "CIB", "CIE", "FBB", "FBE", "PGE", "SKE", "SKB", "WBB", "WBE") #Rename columns
head(salinityData)

#### REFORMAT DATA ####

salinityData$Date.Time <- paste(salinityData$Date, salinityData$Time) #Create a Date.Time column
head(salinityData) #Confirm changes

#### CALCULATE RANGE OF SALINITY ####

salinityRange <- range(salinityData$WBE, salinityData$WBB, salinityData$SKE, salinityData$SKB, salinityData$PGE, salinityData$CIE, salinityData$CIB, salinityData$FBE, salinityData$FBB, na.rm = TRUE) #Calculate range of salinity values
salinityRange[2] <- 45 #Change maximum value to a round number
salinityRange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

salinityBoxplotCIB <- data.frame(Date.Time = salinityData$Date.Time,
                                 Site = rep(x = "CI", times = length(salinityData$Date.Time)),
                                 Habitat = rep(x = "Bare", times = length(salinityData$Date.Time)),
                                 salinity = salinityData$CIB) #Create a dataframe for CIB data
salinityBoxplotCIB <- salinityBoxplotCIB[-c(7456:7457),] #Remove last two rows
salinityBoxplotCIE <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "CI", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Eelgrass", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$CIE) #Create a dataframe for CIE data
salinityBoxplotCIE <- salinityBoxplotCIE[-c(7489:7490),] #Remove last two rows

salinityBoxplotFBB <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "FB", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Bare", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$FBB) #Create a dataframe for FBB data
salinityBoxplotFBB <- salinityBoxplotFBB[-c(7489:7490),] #Remove last two rows
salinityBoxplotFBE <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "FB", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Eelgrass", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$FBE) #Create a dataframe for FBE data
salinityBoxplotFBE <- salinityBoxplotFBE[-c(7489:7490),] #Remove last two rows

salinityBoxplotPGE <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "PG", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Eelgrass", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$PGE) #Create a dataframe for PGE data
salinityBoxplotPGE <- salinityBoxplotPGE[-c(7489:7490),] #Remove last two rows

salinityBoxplotSKB <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "SK", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Bare", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$SKB) #Create a dataframe for SKB data
salinityBoxplotSKB <- salinityBoxplotSKB[-c(7489:7490),] #Remove last two rows
salinityBoxplotSKE <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "SK", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Eelgrass", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$SKE) #Create a dataframe for SKE data
salinityBoxplotSKE <- salinityBoxplotSKE[-c(7489:7490),] #Remove last two rows

salinityBoxplotWBB <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "WB", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Bare", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$WBB) #Create a dataframe for WBB data
salinityBoxplotWBB <- salinityBoxplotWBB[-c(7489:7490),] #Remove last two rows
salinityBoxplotWBE <- data.frame(Date.Time = salinityData$Date.Time,
                                    Site = rep(x = "WB", times = length(salinityData$Date.Time)),
                                    Habitat = rep(x = "Eelgrass", times = length(salinityData$Date.Time)),
                                    salinity = salinityData$WBE) #Create a dataframe for WBE data
salinityBoxplotWBE <- salinityBoxplotWBE[-c(7489:7490),] #Remove last two rows

salinityBoxplot <- rbind(salinityBoxplotCIB, salinityBoxplotCIE, salinityBoxplotFBB, salinityBoxplotFBE, salinityBoxplotPGE, salinityBoxplotSKB, salinityBoxplotSKE, salinityBoxplotWBB, salinityBoxplotWBE) #Bind together all values.

#### MAKE BOXPLOT JUST BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Salinity-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(salinityBoxplot$salinity ~ salinityBoxplot$Site, ylim = salinityRange, main = "Salinity at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(salinityBoxplot$salinity ~ salinityBoxplot$Site) #Perform an ANOVA to test for significant differences in salinitys between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Salinity (ºC)", cex.lab = 2.5, line = 2.2) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for salinity differences between sites. All pairwise differences are significant at 0.05 level.


#### VISUALIZE DIURNAL FLUCTUATIONS AND BOXPLOT ####
#Because I don't have PGB data, I'm going to make a 5x2 multipanel plot and put the site boxplot where PGB would have gone.

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Diurnal-Salinity-Fluctuations-and-Boxplot.jpeg", height = 6000, width = 4000)

par(mfrow = c(5,2)) #Create multipanel plot with 5 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(salinityData$WBE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Eelgrass") #Willapa Bay, Eelgrass
abline(h = median(salinityData$WBE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$WBE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity
mtext(side = 2, text = "Salinity (ºC)", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(salinityData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Bare") #Willapa Bay, Bare
abline(h = median(salinityData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$FBE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Eelgrass") #Fidgalo Bay, Eelgrass
abline(h = median(salinityData$FBE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$FBE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Bare") #FBB
abline(h = median(salinityData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$SKE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Eelgrass") #Skokomish River Delta, Eelgrass
abline(h = median(salinityData$SKE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$SKE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Bare") #Skokomish River Delta, Bare
abline(h = median(salinityData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$CIE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, type = "l", cex.main = 10, col = "red", main = "Case Inlet Eelgrass") #Case Inlet, Eelgrass
abline(h = median(salinityData$CIE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$CIE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange,  yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "red", main = "Case Inlet Bare") #Case Inlet, Bare
abline(h = median(salinityData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean salinity

plot(salinityData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = salinityRange, cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Eelgrass") #Port Gamble Bay, Eelgrass
abline(h = median(salinityData$PGE, na.rm = TRUE), lty = 1) #Add line depicting median salinity
abline(h = mean(salinityData$PGE, na.rm = TRUE), lty = 2) #Add line depicting mean salinity
axis(side = 1, at = seq(from = 1, to = length(salinityData$Date), by = 144*5), lab = salinityData$Date[seq(from = 1, to = length(salinityData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 3) #Modify x-axis labels

boxplot(salinityBoxplot$salinity ~ salinityBoxplot$Site, xaxt = "n", ylim = salinityRange, yaxt = "n", main = "Salinity at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites and habitat
siteANOVA <- aov(salinityBoxplot$salinity ~ salinityBoxplot$Site) #Perform an ANOVA to test for significant differences in salinity between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 3) #Modify x-axis label

dev.off()

#Be sure to clear all plot history to reset par.