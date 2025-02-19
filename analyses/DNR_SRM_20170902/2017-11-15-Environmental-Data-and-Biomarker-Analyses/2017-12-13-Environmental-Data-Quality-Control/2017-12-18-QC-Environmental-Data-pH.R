#In this script, I'll visualize pH data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

pHData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-pH-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA")
head(pHData)

#### SUBSET DATA ####
#I only want the pH data.

colnames(pHData)
pHData <- pHData[,-c(1, 10:14)] #Save pH data as a new dataframe
head(pHData)
colnames(pHData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
pHData$Date <- as.Date(pHData$Date) #Recognize dates
pHData <- pHData[pHData$Date >= "2016-06-19", ] #Subset only data after outplant start date
head(pHData) #Confirm changes

#### CALCULATE RANGE OF pHS ####

pHRange <- range(pHData$WBB, pHData$SKB, pHData$PGB, pHData$CIB, pHData$FBB, na.rm = TRUE) #Calculate range of pH values
pHRange[1] <- 6.5 #Change minimum value to a round number
pHRange[2] <- 8.5 #Change maximum value to a round number
pHRange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

pHBoxplotCIB <- data.frame(Date.Time = pHData$DateTime,
                                 Site = rep(x = "CI", times = length(pHData$DateTime)),
                                 pH = pHData$CIB) #Create a dataframe for CIB data
pHBoxplotCIB <- pHBoxplotCIB[-c(7489:7490),] #Remove last two rows

pHBoxplotFBB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "FB", times = length(pHData$DateTime)),
                                    pH = pHData$FBB) #Create a dataframe for FBB data
pHBoxplotFBB <- pHBoxplotFBB[-c(7489:7490),] #Remove last two rows

pHBoxplotPGB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "PG", times = length(pHData$DateTime)),
                                    pH = pHData$PGB) #Create a dataframe for PGB data
pHBoxplotPGB <- pHBoxplotPGB[-c(7489:7490),] #Remove last two rows

pHBoxplotSKB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "SK", times = length(pHData$DateTime)),
                                    pH = pHData$SKB) #Create a dataframe for SKB data
pHBoxplotSKB <- pHBoxplotSKB[-c(7489:7490),] #Remove last two rows

pHBoxplotWBB <- data.frame(Date.Time = pHData$DateTime,
                                    Site = rep(x = "WB", times = length(pHData$DateTime)),
                                    pH = pHData$WBB) #Create a dataframe for WBB data
pHBoxplotWBB <- pHBoxplotWBB[-c(7489:7490),] #Remove last two rows

pHBoxplot <- rbind(pHBoxplotCIB, pHBoxplotFBB, pHBoxplotPGB, pHBoxplotSKB, pHBoxplotWBB)

#### MAKE BOXPLOT BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-pH-QC-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(pHBoxplot$pH ~ pHBoxplot$Site, ylim = pHRange, main = "pH at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(pHBoxplot$pH ~ pHBoxplot$Site) #Perform an ANOVA to test for significant differences in pHs between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "pH", cex.lab = 2.5, line = 2.5) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for pH differences between sites. All pairwise differences are significant at 0.05 level.

#### VISUALIZE DIURNAL FLUCTUATIONS AND BOXPLOT ####
#I'm going to make a 3x2 multipanel plot and put the site boxplot in the bottom right corner.

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-pH-QC-Fluctuations-and-Boxplot.jpeg", height = 5000, width = 4000)

par(mfrow = c(3,2)) #Create multipanel plot with 3 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(pHData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, type = "l", cex.main = 10, cex.axis = 5, col = "red", main = "Case Inlet") #Case Inlet
abline(h = median(pHData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean pH
mtext(side = 2, text = "pH", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(pHData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay") #Fidalgo Bay
abline(h = median(pHData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$PGB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, cex.axis = 5, type = "l", col = "magenta", main = "Port Gamble Bay") #Port Gamble Bay
abline(h = median(pHData$PGB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$PGB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta") #Skokomish River Delta
abline(h = median(pHData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean pH

plot(pHData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = pHRange, cex.main = 10, cex.axis = 5, type = "l", col = "dark grey", main = "Willapa Bay") #Willapa Bay
abline(h = median(pHData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median pH
abline(h = mean(pHData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean pH
axis(side = 1, at = seq(from = 1, to = length(pHData$Date), by = 144*5), lab = pHData$Date[seq(from = 1, to = length(pHData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 7) #Modify x-axis labels

boxplot(pHBoxplot$pH ~ pHBoxplot$Site, xaxt = "n", ylim = pHRange, yaxt = "n", main = "pH at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites
siteANOVA <- aov(pHBoxplot$pH ~ pHBoxplot$Site) #Perform an ANOVA to test for significant differences in pH between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 7) #Modify x-axis label

#dev.off()

#Be sure to clear all plot history to reset par.