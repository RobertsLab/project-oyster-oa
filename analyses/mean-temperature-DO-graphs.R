# In this script, I will create graphs to compare mean temperature and mean DO at each sample site.

# Step 1: Import data
waterChemistry <- read.csv("../data/Environmental Summary Data for Proteomics Project.csv", header = TRUE, na.strings = NA) #Read in the water chemistry data and save as waterChemistry. Interpret "NA" as NA. The column descriptions are as follows: Outplant (Round 1 (1) or Round 2 (2)), Site (Fidalgo Bay (FB), Port Gamble Bay (PG), Case Inlet (CI), Skokomish River Delta (SK) or Willapa Bay (WB)), Habitat (Bare (B) or Eelgrass (E), Mean.Temp. (mean temperature in celsius), Mean..DO. (mean dissolved oxygen))
waterChemistry #view data

#Step 2: Plot mean temperature for all sites (Round 1 and Round 2)

allDataColors <- c(rep(c("lightskyblue1", "darkseagreen1"), 5), rep(c("skyblue1", "seagreen1"), 5)) #create a vector of colors alternating blue (no eelgrass) and green (eelgrass). first 5 repetitions are for Round 1, second 5 are for Round 2.
barplot(waterChemistry$Mean.Temp., names.arg = waterChemistry$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), col = allDataColors, main = "Mean Temperature for all Sites and Both Rounds") #plot data
legend(x = -1, y = 20.5, fill = c("lightskyblue1", "darkseagreen1", "skyblue1", "seagreen1"), legend = c("Round 1 Bare", "Round 1 Eelgrass", "Round 2 Bare", "Round 2 Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 3: Plot mean temperature for Round 1 sites

round1Colors <- rep(c("lightskyblue1", "darkseagreen1"), 5) #create a vector of colors for round 1
# waterChemistry[waterChemistry$Outplant == 1,] subsets only Round 1 (Outplant 1) from overall waterChemistry data
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), main = "Mean Temperature for all Sites during Round 1", col = round1Colors) #Barplot for data from only Round 1
legend(x = "topleft", fill = c("lightskyblue1", "darkseagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 4: Plot mean temperature for Round 2 sites

round2Colors <- rep(c("skyblue1", "seagreen1"), 5)
# waterChemistry[waterChemistry$Outplant == 2,] subsets only Round 2 from overall waterChemistry data
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), main = "Mean Temperature for all Sites during Round 2", col = round2Colors) #Barplot for data from only Round 2
legend(x = "topleft", fill = c("skyblue1", "seagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 5: Plot mean dissolved oxygen (DO) for all sites (Round 1 and Round 2)

allDataColors <- c(rep(c("lightskyblue1", "darkseagreen1"), 5), rep(c("skyblue1", "seagreen1"), 5)) #create a vector of colors alternating blue (no eelgrass) and green (eelgrass). first 5 repetitions are for Round 1, second 5 are for Round 2.
barplot(waterChemistry$Mean..DO., names.arg = waterChemistry$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen", ylim = c(0,20), col = allDataColors, main = "Mean Dissolved Oxygen for all Sites and Both Rounds") #plot data
legend(x = -1, y = 20.5, fill = c("lightskyblue1", "darkseagreen1", "skyblue1", "seagreen1"), legend = c("Round 1 Bare", "Round 1 Eelgrass", "Round 2 Bare", "Round 2 Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 6: Plot mean DO for Round 1 sites

round1Colors <- rep(c("lightskyblue1", "darkseagreen1"), 5) #create a vector of colors for round 1
# waterChemistry[waterChemistry$Outplant == 1,] subsets only Round 1 (Outplant 1) from overall waterChemistry data
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen", ylim = c(0,20), main = "Mean Dissolved Oxygen for all Sites during Round 1", col = round1Colors) #Barplot for data from only Round 1
legend(x = "topleft", fill = c("lightskyblue1", "darkseagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 7: Plot mean DO for Round 2 sites

round2Colors <- rep(c("skyblue1", "seagreen1"), 5)
# waterChemistry[waterChemistry$Outplant == 2,] subsets only Round 2 from overall waterChemistry data
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen", ylim = c(0,20), main = "Mean Dissolved Oxygen for all Sites during Round 2", col = round2Colors) #Barplot for data from only Round 2
legend(x = "topleft", fill = c("skyblue1", "seagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend

#Step 8: Call Steps 3, 4, 6 and 7 in one composite graph
par(mfrow = c(2,2), mar = c(2,1,2,0), oma = c(5, 5, 4, 1))
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), col = round1Colors, main = "Mean Temperatures for Round 1", xaxt = "n") #Barplot for temperature data from only Round 1
mtext(text = "Mean Temperature (ºC)", side = 2, line = 3) #add y-axis label
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), col = round2Colors, xaxt = "n", yaxt = "n", main = "Mean Temperatures for Round 2") #Barplot for temperature data from only Round 2
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen", ylim = c(0,20), col = round1Colors, main = "Mean Dissolved Oxygen for Round 1") #Barplot for dissolved oxygen data from only Round 1
legend(x = 7, y = 23, fill = c("lightskyblue1", "darkseagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend
mtext(text = "Mean Dissolved Oxygen", side = 2, line = 3) #add y-axis label
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen", ylim = c(0,20), col = round2Colors, main = "Mean Dissolved Oxygen for Round 2", yaxt = "n") #Barplot for dissolve oxygen data from only Round 2
legend(x = 0, y = 23, fill = c("skyblue1", "seagreen1"), legend = c("Bare", "Eelgrass"), bty = "n", y.intersp = .5) #add legend
mtext(text = "Sample Site", side = 1, outer =  TRUE, line = 1.5) #add x-axis label

# Step 9: Save composite graph as a .pdf
pdf(file = "Mean-Temperature-DO.pdf", width = 11, height = 8) #create .pdf to save graph
par(mfrow = c(2,2), mar = c(2,1,2,0), oma = c(5, 5, 4, 1))
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), col = round1Colors, main = "Mean Temperatures for Round 1", xaxt = "n") #Barplot for temperature data from only Round 1
mtext(text = "Mean Temperature (ºC)", side = 2, line = 3) #add y-axis label
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean.Temp., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Temperature (ºC)", ylim = c(0,20), col = round2Colors, xaxt = "n", yaxt = "n", main = "Mean Temperatures for Round 2") #Barplot for temperature data from only Round 2
barplot(waterChemistry[waterChemistry$Outplant == 1,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 1,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen (mg/L)", ylim = c(0,20), col = round1Colors, main = "Mean Dissolved Oxygen for Round 1") #Barplot for dissolved oxygen data from only Round 1
legend(x = "topright", fill = c("lightskyblue1", "darkseagreen1"), legend = c("Bare", "Eelgrass"), bty = "n") #add legend
mtext(text = "Mean Dissolved Oxygen (mg/L)", side = 2, line = 3) #add y-axis label
barplot(waterChemistry[waterChemistry$Outplant == 2,]$Mean..DO., names.arg = waterChemistry[waterChemistry$Outplant == 2,]$Site, xlab = "Sample Sites", ylab = "Mean Dissolved Oxygen (mg/L)", ylim = c(0,20), col = round2Colors, main = "Mean Dissolved Oxygen for Round 2", yaxt = "n") #Barplot for dissolve oxygen data from only Round 2
legend(x = "topleft", fill = c("skyblue1", "seagreen1"), legend = c("Bare", "Eelgrass"), bty = "n") #add legend
mtext(text = "Sample Site", side = 1, outer =  TRUE, line = 1.5) #add x-axis label
dev.off() #stop saving as .pdf
