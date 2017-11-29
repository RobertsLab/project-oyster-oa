#In this script, I'll visualize pH data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

environmentalData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA")
head(environmentalData)

#### SUBSET DATA ####
#I only want the pH data.

colnames(environmentalData)
pHData <- environmentalData[,c(1:12)] #Save pH data as a new dataframe
head(pHData)
colnames(pHData) <- c("DateTime", "Date", "Time", "WBE", "WBB", "SKE", "SKB", "PGB", "CIE", "CIB", "FBE", "FBB") #Rename columns
head(pHData)

#### CALCULATE RANGE OF pHS ####

pHRange <- range(pHData$WBE, pHData$WBB, pHData$SKE, pHData$SKB, pHData$PGB, pHData$CIE, pHData$CIB, pHData$FBE, pHData$FBB, na.rm = TRUE) #Calculate range of pH values
pHRange[1] <- 6.4 #Change minimum value to a round number
pHRange[2] <- 8.8 #Change maximum value to a round number
pHRange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

pHBoxplotCIB <- data.frame(Date.Time = pHData$DateTime,
                                 Site = rep(x = "CI", times = length(pHData$DateTime)),
                                 Habitat = rep(x = "Bare", times = length(pHData$DateTime)),
                                 pH = pHData$CIB) #Create a dataframe for CIB data
pHBoxplotCIB <- pHBoxplotCIB[-c(7489:7490),] #Remove last two rows
pHBoxplotCIE <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "CI", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(pHData$DateTime)),
                                    pH = pHData$CIE) #Create a dataframe for CIE data
pHBoxplotCIE <- pHBoxplotCIE[-c(7489:7490),] #Remove last two rows

pHBoxplotFBB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "FB", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(pHData$DateTime)),
                                    pH = pHData$FBB) #Create a dataframe for FBB data
pHBoxplotFBB <- pHBoxplotFBB[-c(7489:7490),] #Remove last two rows
pHBoxplotFBE <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "FB", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(pHData$DateTime)),
                                    pH = pHData$FBE) #Create a dataframe for FBE data
pHBoxplotFBE <- pHBoxplotFBE[-c(7489:7490),] #Remove last two rows

pHBoxplotPGB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "PG", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(pHData$DateTime)),
                                    pH = pHData$PGB) #Create a dataframe for PGB data
pHBoxplotPGB <- pHBoxplotPGB[-c(7489:7490),] #Remove last two rows

pHBoxplotSKB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "SK", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(pHData$DateTime)),
                                    pH = pHData$SKB) #Create a dataframe for SKB data
pHBoxplotSKB <- pHBoxplotSKB[-c(7489:7490),] #Remove last two rows
pHBoxplotSKE <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "SK", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(pHData$DateTime)),
                                    pH = pHData$SKE) #Create a dataframe for SKE data
pHBoxplotSKE <- pHBoxplotSKE[-c(7489:7490),] #Remove last two rows

pHBoxplotWBB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "WB", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(pHData$DateTime)),
                                    pH = pHData$WBB) #Create a dataframe for WBB data
pHBoxplotWBB <- pHBoxplotWBB[-c(7489:7490),] #Remove last two rows
pHBoxplotWBE <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "WB", times = length(pHData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(pHData$DateTime)),
                                    pH = pHData$WBE) #Create a dataframe for WBE data
pHBoxplotWBE <- pHBoxplotWBE[-c(7489:7490),] #Remove last two rows

pHBoxplot <- rbind(pHBoxplotCIB, pHBoxplotCIE, pHBoxplotFBB, pHBoxplotFBE, pHBoxplotPGB, pHBoxplotSKB, pHBoxplotSKE, pHBoxplotWBB, pHBoxplotWBE)

#### MAKE BOXPLOT JUST BASED ON SITES ####

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-pH-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(pHBoxplot$pH ~ pHBoxplot$Site, ylim = pHRange, main = "pH at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(pHBoxplot$pH ~ pHBoxplot$Site) #Perform an ANOVA to test for significant differences in pHs between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "pH", cex.lab = 2.5, line = 2.5) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for pH differences between sites. All pairwise differences are significant at 0.05 level.

#### VISUALIZE DIURNAL FLUCTUATIONS AND BOXPLOT ####
#Because I don't have PGB data, I'm going to make a 5x2 multipanel plot and put the site boxplot where PGB would have gone.

jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Diurnal-pH-Fluctuations-and-Boxplot.jpeg", height = 6000, width = 4000)

par(mfrow = c(5,2)) #Create multipanel plot with 5 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(pHData$WBE, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Eelgrass") #Willapa Bay, Eelgrass
abline(h = median(pHData$WBE, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$WBE, na.rm = TRUE), lty = 2) #Add line depicting mean pH
mtext(side = 2, text = "pH (ºC)", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(pHData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Bare") #Willapa Bay, Bare
abline(h = median(pHData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$FBE, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Eelgrass") #Fidgalo Bay, Eelgrass
abline(h = median(pHData$FBE, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$FBE, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Bare") #FBB
abline(h = median(pHData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$SKE, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Eelgrass") #Skokomish River Delta, Eelgrass
abline(h = median(pHData$SKE, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$SKE, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Bare") #Skokomish River Delta, Bare
abline(h = median(pHData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$CIE, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, type = "l", cex.main = 10, col = "red", main = "Case Inlet Eelgrass") #Case Inlet, Eelgrass
abline(h = median(pHData$CIE, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$CIE, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange,  yaxt = "n", cex.axis = 5, cex.main = 10, type = "l", col = "red", main = "Case Inlet Bare") #Case Inlet, Bare
abline(h = median(pHData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Eelgrass") #Port Gamble Bay, Eelgrass
abline(h = median(pHData$PGE, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$PGE, na.rm = TRUE), lty = 2) #Add line depicting mean pH
axis(side = 1, at = seq(from = 1, to = length(pHData$Date), by = 144*5), lab = pHData$Date[seq(from = 1, to = length(pHData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 7) #Modify x-axis labels

boxplot(pHBoxplot$pH ~ pHBoxplot$Site, xaxt = "n", ylim = pHRange, yaxt = "n", main = "pH at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites and habitat
siteANOVA <- aov(pHBoxplot$pH ~ pHBoxplot$Site) #Perform an ANOVA to test for significant differences in pH between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 7) #Modify x-axis label

dev.off()

#Be sure to clear all plot history to reset par.