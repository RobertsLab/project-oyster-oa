#In this script, I'll visualize DO data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("../..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

DOData <- read.csv("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-DO-Data-QC-with-Tide-Data.csv", header = TRUE, na.strings = "NA")
head(DOData)

#### SUBSET DATA ####
#I only want the DO data.

colnames(DOData)
DOData <- DOData[,-c(1, 10:14)] #Save DO data as a new dataframe
head(DOData)
colnames(DOData) <- c("DateTime", "Date", "Time", "WBB", "SKB", "PGB", "CIB", "FBB") #Rename columns
DOData$Date <- as.Date(DOData$Date) #Recognize dates
DOData <- DOData[DOData$Date >= "2016-06-19", ] #Subset data from after outplant start date
head(DOData)

#### CALCULATE RANGE OF DOS ####

DORange <- range(DOData[, 4:8], na.rm = TRUE) #Calculate range of DO values
DORange[1] <- 0 #Change minimum value
DORange[2] <- 25 #Change maximum value to a round number
DORange #Confirm changes

#### REFORMAT DATA FOR BOXPLOT ####

DOBoxplotCIB <- data.frame(Date.Time = DOData$DateTime,
                           Site = rep(x = "CI", times = length(DOData$DateTime)),
                           DO = DOData$CIB) #Create a dataframe for CIB data
DOBoxplotCIB <- DOBoxplotCIB[-c(7489:7490),] #Remove last two rows

DOBoxplotFBB <- data.frame(Date.Time = DOData$DateTime,
                           Site = rep(x = "FB", times = length(DOData$DateTime)),
                           DO = DOData$FBB) #Create a dataframe for FBB data
DOBoxplotFBB <- DOBoxplotFBB[-c(7489:7490),] #Remove last two rows

DOBoxplotPGB <- data.frame(Date.Time = DOData$DateTime,
                           Site = rep(x = "PG", times = length(DOData$DateTime)),
                           DO = DOData$PGB) #Create a dataframe for PGB data
DOBoxplotPGB <- DOBoxplotPGB[-c(7489:7490),] #Remove last two rows

DOBoxplotSKB <- data.frame(Date.Time = DOData$DateTime,
                           Site = rep(x = "SK", times = length(DOData$DateTime)),
                           DO = DOData$SKB) #Create a dataframe for SKB data
DOBoxplotSKB <- DOBoxplotSKB[-c(7489:7490),] #Remove last two rows

DOBoxplotWBB <- data.frame(Date.Time = DOData$DateTime,
                           Site = rep(x = "WB", times = length(DOData$DateTime)),
                           DO = DOData$WBB) #Create a dataframe for WBB data
DOBoxplotWBB <- DOBoxplotWBB[-c(7489:7490),] #Remove last two rows

DOBoxplot <- rbind(DOBoxplotCIB, DOBoxplotFBB, DOBoxplotPGB, DOBoxplotSKB, DOBoxplotWBB)

#### MAKE BOXPLOT JUST BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-DO-QC-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(DOBoxplot$DO ~ DOBoxplot$Site, ylim = DORange, main = "DO at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(DOBoxplot$DO ~ DOBoxplot$Site) #Perform an ANOVA to test for significant differences in DOs between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Dissolved Oxygen Content", cex.lab = 2.5, line = 2.5) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for DO differences between sites. Significant differences between all pairwise combinations

#### VISUALIZE DIURNAL FLUCTUATIONS ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control/2017-12-18-Diurnal-DO-QC-Fluctuations.jpeg", height = 5000, width = 4000)

par(mfrow = c(3,2)) #Create multipanel plot with 3 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(40, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(DOData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, type = "l", cex.main = 10, cex.axis = 5, col = "red", main = "Case Inlet") #Case Inlet
abline(h = median(DOData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean DO
mtext(side = 2, text = "DO", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels

plot(DOData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay") #Fidalgo Bay
abline(h = median(DOData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$PGB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.main = 10, cex.axis = 5, type = "l", col = "magenta", main = "Port Gamble Bay") #Port Gamble Bay
abline(h = median(DOData$PGB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$PGB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta") #Skokomish River Delta
abline(h = median(DOData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.main = 10, cex.axis = 5, type = "l", col = "dark grey", main = "Willapa Bay") #Willapa Bay
abline(h = median(DOData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean DO
axis(side = 1, at = seq(from = 1, to = length(DOData$Date), by = 144*5), lab = DOData$Date[seq(from = 1, to = length(DOData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 35, cex = 7) #Modify x-axis labels

boxplot(DOBoxplot$DO ~ DOBoxplot$Site, xaxt = "n", ylim = DORange, yaxt = "n", main = "Dissolved Oxygen at Sites", cex.main = 10, cex.axis = 5, line.axis = 2, col = c("red", "blue", "magenta", "green", "dark gray")) #Make boxplot based on sites
siteANOVA <- aov(DOBoxplot$DO ~ DOBoxplot$Site) #Perform an ANOVA to test for significant differences in DO between sites
legend("topright", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
axis(side = 1, at = 1:5, lab = c("CI", "FB", "PG", "SK", "WB"), cex.axis = 5, line = 10, lwd = 0, lwd.ticks = 0) #Make x-axis
mtext(side = 1, text = "Site", line = 35, cex = 7) #Modify x-axis label

#dev.off()

#Be sure to clear all plot history to reset par.