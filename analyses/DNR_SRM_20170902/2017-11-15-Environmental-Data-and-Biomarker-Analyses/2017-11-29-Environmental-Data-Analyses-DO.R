#In this script, I'll visualize DO data from Micah between sites and bare and eelgrass habitats I'll examine diurnal fluctuations as well as an overall boxplot to see variation between sites and habitats over the course of the outplant. 

#### SET WORKING DIRECTORY ####
setwd("..") #Set working directory to the master SRM folder
getwd()

#### IMPORT DATA ####

environmentalData <- read.csv("../../data/DNR/2017-11-14-Environmental-Data-from-Micah.csv", header = TRUE, na.strings = "NA")
head(environmentalData)

#### SUBSET DATA ####
#I only want the DO data.

colnames(environmentalData)
DOData <- environmentalData[,c(1:3,22:31)] #Save DO data as a new dataframe
head(DOData)
colnames(DOData) <- c("DateTime", "Date", "Time", "WBE", "WBB", "SKE", "SKB", "PGE", "PGB", "CIE", "CIB", "FBE", "FBB") #Rename columns
head(DOData)

#### CALCULATE RANGE OF DOS ####

DORange <- range(DOData$WBE, DOData$WBB, DOData$SKE, DOData$SKB, DOData$PGE, DOData$PGB, DOData$CIE, DOData$CIB, DOData$FBE, DOData$FBB, na.rm = TRUE) #Calculate range of DO values
DORange[1] <- 0 #Change minimum value
DORange[2] <- 30 #Change maximum value to a believeable number
DORange #Confirm changes

#### VISUALIZE DIURNAL FLUCTUATIONS ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Diurnal-DO-Fluctuations.jpeg", height = 6000, width = 4000)

par(mfrow = c(5,2)) #Create multipanel plot with 5 rows and 2 columns
par(mar = c(0, 0, 10, 0), oma = c(30, 15, 1, 1)) #Remove redundant white space and change outer margins

plot(DOData$FBB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.axis = 5, cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Bare") #FBB. Set up plot with no x axis labels, but with y-axis that encompasses maximum and minimum values
abline(h = median(DOData$FBB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$FBB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$FBE, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "blue", main = "Fidalgo Bay Eelgrass") #Fidgalo Bay, Eelgrass
abline(h = median(DOData$FBE, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$FBE, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$PGB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.axis = 5, cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Bare") #Port Gamble Bay, Bare
abline(h = median(DOData$PGB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$PGB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$PGE, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "magenta", main = "Port Gamble Bay Eelgrass") #Port Gamble Bay, Eelgrass
abline(h = median(DOData$PGE, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$PGE, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$SKB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.axis = 5, cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Bare") #Skokomish River Delta, Bare
abline(h = median(DOData$SKB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$SKB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$SKE, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "green", main = "Skokomish River Delta Eelgrass") #Skokomish River Delta, Eelgrass
abline(h = median(DOData$SKE, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$SKE, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$CIB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.axis = 5, cex.main = 10, type = "l", col = "red", main = "Case Inlet Bare") #Case Inlet, Bare
abline(h = median(DOData$CIB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$CIB, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$CIE, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", type = "l", cex.main = 10, col = "red", main = "Case Inlet Eelgrass") #Case Inlet, Eelgrass
abline(h = median(DOData$CIE, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$CIE, na.rm = TRUE), lty = 2) #Add line depicting mean DO

plot(DOData$WBB, xlab = "", xaxt = "n", ylab = "", ylim = DORange, cex.axis = 5, cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Bare") #Willapa Bay, Bare
abline(h = median(DOData$WBB, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$WBB, na.rm = TRUE), lty = 2) #Add line depicting mean DO
mtext(side = 2, text = "Dissolved Oxygen Content", line = 7, cex = 5, outer = TRUE) #Modify y-axis labels
axis(side = 1, at = seq(from = 1, to = length(DOData$Date), by = 144*5), lab = DOData$Date[seq(from = 1, to = length(DOData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis

plot(DOData$WBE, xlab = "", xaxt = "n", ylab = "", ylim = DORange, yaxt = "n", cex.main = 10, type = "l", col = "dark grey", main = "Willapa Bay Eelgrass") #Willapa Bay, Eelgrass
abline(h = median(DOData$WBE, na.rm = TRUE), lty = 1) #Add line depicting median DO
abline(h = mean(DOData$WBE, na.rm = TRUE), lty = 2) #Add line depicting mean DO
axis(side = 1, at = seq(from = 1, to = length(DOData$Date), by = 144*5), lab = DOData$Date[seq(from = 1, to = length(DOData$Date), by = 144*5)], las = 3, cex.axis = 5, line = 2) #Make x-axis
mtext(side = 1, text = "Date", line = 25, cex = 5, outer = TRUE) #Modify x-axis labels

#dev.off()

#Be sure to clear all plot history to reset par.

#### REFORMAT DATA FOR BOXPLOT ####

DOBoxplotCIB <- data.frame(Date.Time = DOData$DateTime,
                                 Site = rep(x = "CI", times = length(DOData$DateTime)),
                                 Habitat = rep(x = "Bare", times = length(DOData$DateTime)),
                                 DO = DOData$CIB) #Create a dataframe for CIB data
DOBoxplotCIB <- DOBoxplotCIB[-c(7489:7490),] #Remove last two rows
DOBoxplotCIE <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "CI", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(DOData$DateTime)),
                                    DO = DOData$CIE) #Create a dataframe for CIE data
DOBoxplotCIE <- DOBoxplotCIE[-c(7489:7490),] #Remove last two rows

DOBoxplotFBB <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "FB", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(DOData$DateTime)),
                                    DO = DOData$FBB) #Create a dataframe for FBB data
DOBoxplotFBB <- DOBoxplotFBB[-c(7489:7490),] #Remove last two rows
DOBoxplotFBE <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "FB", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(DOData$DateTime)),
                                    DO = DOData$FBE) #Create a dataframe for FBE data
DOBoxplotFBE <- DOBoxplotFBE[-c(7489:7490),] #Remove last two rows

DOBoxplotPGB <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "PG", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(DOData$DateTime)),
                                    DO = DOData$PGB) #Create a dataframe for PGB data
DOBoxplotPGB <- DOBoxplotPGB[-c(7489:7490),] #Remove last two rows
DOBoxplotPGE <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "PG", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(DOData$DateTime)),
                                    DO = DOData$PGE) #Create a dataframe for PGE data
DOBoxplotPGE <- DOBoxplotPGE[-c(7489:7490),] #Remove last two rows

DOBoxplotSKB <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "SK", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(DOData$DateTime)),
                                    DO = DOData$SKB) #Create a dataframe for SKB data
DOBoxplotSKB <- DOBoxplotSKB[-c(7489:7490),] #Remove last two rows
DOBoxplotSKE <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "SK", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(DOData$DateTime)),
                                    DO = DOData$SKE) #Create a dataframe for SKE data
DOBoxplotSKE <- DOBoxplotSKE[-c(7489:7490),] #Remove last two rows

DOBoxplotWBB <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "WB", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Bare", times = length(DOData$DateTime)),
                                    DO = DOData$WBB) #Create a dataframe for WBB data
DOBoxplotWBB <- DOBoxplotWBB[-c(7489:7490),] #Remove last two rows
DOBoxplotWBE <- data.frame(Date.Time = DOData$DateTime,
                                    Site = rep(x = "WB", times = length(DOData$DateTime)),
                                    Habitat = rep(x = "Eelgrass", times = length(DOData$DateTime)),
                                    DO = DOData$WBE) #Create a dataframe for WBE data
DOBoxplotWBE <- DOBoxplotWBE[-c(7489:7490),] #Remove last two rows

DOBoxplot <- rbind(DOBoxplotCIB, DOBoxplotCIE, DOBoxplotFBB, DOBoxplotFBE, DOBoxplotPGB, DOBoxplotPGE, DOBoxplotSKB, DOBoxplotSKE, DOBoxplotWBB, DOBoxplotWBE)

#### MAKE BOXPLOT JUST BASED ON SITES ####

#jpeg("2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-15-DO-Boxplot-Site-Only.jpeg", height = 1000, width = 2000)
boxplot(DOBoxplot$DO ~ DOBoxplot$Site, ylim = DORange, main = "DO at Sites", cex.main = 3, cex.axis = 1.5) #Make boxplot based on sites and habitat
siteANOVA <- aov(DOBoxplot$DO ~ DOBoxplot$Site) #Perform an ANOVA to test for significant differences in DOs between sites
legend("topleft", bty = "n", legend = paste("F =", format(summary(siteANOVA)[[1]][["F value"]][[1]], digits = 4), "p =", format(summary(siteANOVA)[[1]][["Pr(>F)"]][[1]], digits = 4))) #Add F and p-value from ANOVA
title(xlab = "Site", cex.lab = 2.5, line = 3.5) #Add x-axis label
title(ylab = "Dissolved Oxygen Content", cex.lab = 2.5, line = 2.2) #Add y-axis label
#dev.off()

TukeyHSD(siteANOVA) #Tukey HSD post-hoc test for DO differences between sites. All pairwise differences are significant at 0.05 level.