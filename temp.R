par(mfrow = c(4, 5), oma = c(10, 5, 1, 1), mar = c(0, 2, 0, 0)) #Create a 4x5 multipanel plot, filling in the rows first. Add space along certain outer and inner margins.

#Temperature
plot(temperatureData$CIB, xaxs = "i", yaxs = "i", axes = F, ann = F, xlim = c(1, 4896), ylim = temperatureRange, pch = 16, cex = 0.2, col = "#00A9BD") #Case Inlet bare
lines(temperatureData$CIE, lty = 3, col = "grey80") #Case Inlet eelgrass
mtext(side = 3, line = -1, adj = 0, text = "Case Inlet", cex = 0.5) #Add site name. Use adj = 0 to left-justify the text
mtext(side = 3, line = -1, adj = 1, text = "Unvegetated", cex = 0.5, col = "#00A9BD") #Add habitat specificiation name. Use adj = 1 to right-justify the text
mtext(side = 3, line = -1.75, adj = 1, text = "Eelgrass", cex = 0.5, col = "grey80") #Add habitat specificiation name. Use adj = 1 to right-justify the text
box(col = "grey80")
axis(side = 2, las = 1, at = c(15, 35), col = "grey80") #Put the axis labels at the values specified
mtext(side = 2, line = 3, "Temperature (ºC)") #Add environmental variable indication

plot(temperatureData$FBB, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = temperatureRange, pch = 16, cex = 0.2, col = "#38001C") #Fidalgo Bay bare
lines(temperatureData$FBE, lty = 3, col = "grey80") #Eelgrass
mtext(side = 3, line = -1, adj = 0, text = "Fidalgo Bay", cex = 0.5) #Add site name. Use adj = 0 to left-justify the text
mtext(side = 3, line = -1, adj = 1, text = "Unvegetated", cex = 0.5, col = "#38001C") #Add habitat specificiation name. Use adj = 1 to right-justify the text
mtext(side = 3, line = -1.75, adj = 1, text = "Eelgrass", cex = 0.5, col = "grey80") #Add habitat specificiation name. Use adj = 1 to right-justify the text
box(col = "grey80")

plot(temperatureData$PGB, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = temperatureRange, pch = 16, cex = 0.2, col = "#440D82") #Port Gamble bare
lines(temperatureData$PGE, lty = 3, col = "grey80") #Eelgrass
mtext(side = 3, line = -1, adj = 0, text = "Port Gamble Bay", cex = 0.5) #Add site name. Use adj = 0 to left-justify the text
mtext(side = 3, line = -1, adj = 1, text = "Unvegetated", cex = 0.5, col = "#440D82") #Add habitat specificiation name. Use adj = 1 to right-justify the text
mtext(side = 3, line = -1.75, adj = 1, text = "Eelgrass", cex = 0.5, col = "grey80") #Add habitat specificiation name. Use adj = 1 to right-justify the text
box(col = "grey80")

plot(temperatureData$SKB, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = temperatureRange, pch = 16, cex = 0.2, col = "#017A74") #Skokomish River bare
lines(temperatureData$SKE, lty = 3, col = "grey80") #Eelgrass
mtext(side = 3, line = -1, adj = 0, text = "Skokomish River Delta", cex = 0.5) #Add site name. Use adj = 0 to left-justify the text
mtext(side = 3, line = -1, adj = 1, text = "Unvegetated", cex = 0.5, col = "#017A74") #Add habitat specificiation name. Use adj = 1 to right-justify the text
mtext(side = 3, line = -1.75, adj = 1, text = "Eelgrass", cex = 0.5, col = "grey80") #Add habitat specificiation name. Use adj = 1 to right-justify the text
box(col = "grey80")

plot(temperatureData$WBB, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = temperatureRange, pch = 16, cex = 0.2, col = "#EB8B0C") #Willapa Bay bare
lines(temperatureData$WBE, lty = 3, col = "grey80") #Eelgrass
mtext(side = 3, line = -1, adj = 0, text = "Willapa Bay", cex = 0.5) #Add site name. Use adj = 0 to left-justify the text
mtext(side = 3, line = -1, adj = 1, text = "Unvegetated", cex = 0.5, col = "#EB8B0C") #Add habitat specificiation name. Use adj = 1 to right-justify the text
mtext(side = 3, line = -1.75, adj = 1, text = "Eelgrass", cex = 0.5, col = "grey80") #Add habitat specificiation name. Use adj = 1 to right-justify the text
box(col = "grey80")

#pH data
plot(pHData$CIB.pH, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = pHRange, pch = 16, cex = 0.2, col = "#00A9BD") #Case Inlet bare
lines(pHData$CIE.pH, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 2, las = 1, at = c(7.5, 8.5), col = "grey80") #Put the axis labels at the values specified
mtext(side = 2, line = 3, "pH") #Add environmental variable indication

plot(pHData$FBB.pH, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = pHRange, pch = 16, cex = 0.2, col = "#38001C") #Fidalgo Bay bare
lines(pHData$FBE.pH, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

plot(pHData$PGB.pH, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = pHRange, pch = 16, cex = 0.2, col = "#440D82") #Port Gamble bare
box(col = "grey80")

plot(pHData$SKB.pH, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = pHRange, pch = 16, cex = 0.2, col = "#017A74") #Skokomish River bare
lines(pHData$SKE.pH, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

plot(pHData$WBB.pH, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = pHRange, pch = 16, cex = 0.2, col = "#EB8B0C") #Willapa Bay bare
lines(pHData$WBE.pH, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

#Salinity
plot(salinityData$CIB.Salinity, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = salinityRange, pch = 16, cex = 0.2, col = "#00A9BD") #Case Inlet bare
lines(salinityData$CIE.Salinity, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 2, las = 1, at = c(15, 25), col = "grey80") #Put the axis labels at the values specified
mtext(side = 2, line = 3, "Salinity (PSU)") #Add environmental variable indication

plot(salinityData$FBB.Salinity, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = salinityRange, pch = 16, cex = 0.2, col = "#38001C") #Fidalgo Bay bare
lines(salinityData$FBE.Salinity, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

plot(salinityData$PGE.Salinity, type = "l", lty = 3, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = salinityRange, col = "grey80") #Port Gamble Bay bare
lines(salinityData$PGE.Salinity, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

plot(salinityData$SKB.Salinity, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = salinityRange, pch = 16, cex = 0.2, col = "#017A74") #Skokomish River bare
lines(salinityData$SKE.Salinity, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")

plot(salinityData$WBB.Salinity, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = salinityRange, pch = 16, cex = 0.2, col = "#EB8B0C") #Willapa Bay bare
box(col = "grey80")

#Dissolved oxygen
plot(DOData$CIB.DO, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = DORange, pch = 16, cex = 0.2, col = "#00A9BD") #Case Inlet bare
lines(DOData$CIE.DO, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 2, las = 1, at = c(25, 45), col = "grey80") #Put the axis labels at the values specified
mtext(side = 2, line = 3, "DO (mg/L)") #Add environmental variable indication
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, col = "grey80") #Make x-axis

plot(DOData$FBB.DO, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = DORange, pch = 16, cex = 0.2, col = "#38001C") #Fidalgo Bay bare
lines(DOData$FBE.DO, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, col = "grey80") #Make x-axis

plot(DOData$PGB.DO, xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = DORange, pch = 16, cex = 0.2, col = "#440D82") #Port Gamble bare
lines(DOData$PGE.DO, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, col = "grey80") #Make x-axis

plot(DOData$SKB.DO, type = "l", xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = DORange, pch = 16, cex = 0.2, col = "#017A74") #Skokomish River bare
lines(DOData$SKE.DO, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, col = "grey80") #Make x-axis

plot(DOData$WBB.DO, type = "l", xaxs = "i", yaxs = "i", axes = F, ann = F, ylim = DORange, pch = 16, cex = 0.2, col = "#EB8B0C") #Willapa Bay bare
lines(DOData$WBE.DO, lty = 3, col = "grey80") #Eelgrass
box(col = "grey80")
axis(side = 1, at = seq(from = 1, to = length(temperatureData$Date), by = 144*5), lab = temperatureData$Date[seq(from = 1, to = length(temperatureData$Date), by = 144*5)], las = 3, col = "grey80") #Make x-axis
mtext(side = 1, outer = TRUE, line = 7, "Date")