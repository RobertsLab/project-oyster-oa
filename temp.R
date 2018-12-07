plotMatrix <- matrix(c(1, 2, 3,
                       1, 4, 5,
                       1, 6, 7), nrow = 3, ncol = 3, byrow = TRUE) #Create a matrix and fill it in by row.
plotMatrix #Confirm matrix creation

pdf("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-12-01-Multipanel-Ordination.pdf", width = 11, height = 8.5)

par(mar = c(3, 0, 0, 2.5), oma = c(0, 3.5, 1, 0)) #Specify inner and outer margins
layout(mat = plotMatrix, width = c(45, 15, 15)) #Create a layout based on the plot matrix. Column 1s width should be 3x as large as column 2. 
#layout.show(n = 7) #Confirm the plot layout is what is desired

#### PLOT 1 ####

plot(protEnvRDA, choices=c(1,2), type = 'none', scaling = 2, xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(-1, 1)) #Create an empty plot based on RDA dimensions
points(protEnvRDA, choices = c(1,2), display = 'wa', pch = NMDSColorShapeCustomization$Shape, cex = 2, scaling = 2, col = NMDSColorShapeCustomization$Color) #Plot objects as points
points(protEnvRDA, choices = c(1,2), display = 'sp', col = 'grey20', pch = 4, cex = 2, scaling = 2, select = c(4:5, 12:14, 19, 22, 30, 33:34)) #Plot significant proteins
text(protEnvRDA, choices = c(1,2), display = 'bp', col = 'grey20', cex = 1, select = 1:2) #Plot only marginally significant predictors
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -11, at = c(-1, -10), text = "RDA") #Add test name

#### PLOT 2 ####

fig.nmds <- ordiplot(proc.nmds.averaged.euclidean, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object
points(fig.nmds, "sites", col = NMDSColorShapeCustomization$Color, pch = NMDSColorShapeCustomization$Shape, cex = 2)
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -2, adj = -1, text = "  NMDS: Stress = 0.075") #Add stress value

ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "CI", col = "#00A9BD88") #Add confidence ellipse around the oyster samples from Case Inlet
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "FB", col = "#38001C88") #Add confidence ellipse around the oyster samples from Fidalgo Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "PG", col = "#440D8288") #Add confidence ellipse around the oyster samples from Port Gamble Bay
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "SK", col = "#017A7488") #Add confidence ellipse around the oyster samples from Skokomish River Delta
ordiellipse(proc.nmds.averaged.euclidean, NMDSColorShapeCustomization$Site, show.groups = "WB", col = "#EB8B0C88") #Add confidence ellipse around the oyster samples from Willapa Bay

#### PLOT 3 ####

ordiplot(proc.nmds.averaged.euclidean, choices = c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Create an empty plot
plot(sigProtLoadings, col = 'grey20', labels = c("  1", "1", "    2", "2", "", "4", "5 ", "6", "   6", "3 6"), cex = 1) #Plot loadings that simper determined were significant
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)

#### PLOT 4 ####

fig.nmds2 <- ordiplot(meanData.log.gower.NMDS, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object
points(fig.nmds2, "sites", col = plotCustomization$Color, pch = plotCustomization$Shape2, bg = plotCustomization$Color, cex = 2)
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -2, adj = -1, text = "  NMDS: Stress = 0.017") #Add stress value

ordiellipse(meanData.log.gower.NMDS, plotCustomization$Site, show.groups = "CI", col = "#00A9BD88") #Add confidence ellipse around the data from Case Inlet
ordiellipse(meanData.log.gower.NMDS, plotCustomization$Site, show.groups = "FB", col = "#38001C88") #Add confidence ellipse around the data  from Fidalgo Bay
ordiellipse(meanData.log.gower.NMDS, plotCustomization$Site, show.groups = "PG", col = "#440D8288") #Add confidence ellipse around the data from Port Gamble Bay
ordiellipse(meanData.log.gower.NMDS, plotCustomization$Site, show.groups = "SK", col = "#017A7488") #Add confidence ellipse around the data from Skokomish River Delta
ordiellipse(meanData.log.gower.NMDS, plotCustomization$Site, show.groups = "WB", col = "#EB8B0C88") #Add confidence ellipse around the data from Willapa Bay

#### PLOT 5 ####

fig.nmds2 <- ordiplot(meanData.log.gower.NMDS, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object
plot(sigMeanLoadings, col = 'grey20', labels = c("1", "2", "   3", "", "5", "64", "  7", "18", "19")) #Plot loadings that simper determined were significant
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)

#### PLOT 6 ####

fig.nmds3 <- ordiplot(varData.log.gower.NMDS, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n") #Save NMDS as a new object
points(fig.nmds3, "sites", col = plotCustomization$Color, pch = plotCustomization$Shape2, bg = plotCustomization$Color, cex = 2)
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -2, adj = -1, text = "  NMDS: Stress = 0.034") #Add stress value

ordiellipse(varData.log.gower.NMDS, plotCustomization$Site, show.groups = "CI", col = "#00A9BD88") #Add confidence ellipse around the data from Case Inlet
ordiellipse(varData.log.gower.NMDS, plotCustomization$Site, show.groups = "FB", col = "#38001C88") #Add confidence ellipse around the data  from Fidalgo Bay
ordiellipse(varData.log.gower.NMDS, plotCustomization$Site, show.groups = "PG", col = "#440D8288") #Add confidence ellipse around the data from Port Gamble Bay
ordiellipse(varData.log.gower.NMDS, plotCustomization$Site, show.groups = "SK", col = "#017A7488") #Add confidence ellipse around the data from Skokomish River Delta
ordiellipse(varData.log.gower.NMDS, plotCustomization$Site, show.groups = "WB", col = "#EB8B0C88") #Add confidence ellipse around the data from Willapa Bay

#### PLOT 7 ####

fig.nmds3 <- ordiplot(varData.log.gower.NMDS, choices=c(1,2), type = "none", display = "sites", xlab = "", ylab = "", cex = 0.5, xaxt = "n", yaxt = "n", xlim = c(-0.4, 0.9)) #Save NMDS as a new object
plot(sigVarLoadings, col = 'grey20', labels = c("1 3", "2", "", "4 18", "5", "", "     21", "29", "30", "31")) #Plot loadings that simper determined were significant
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
box(col = "white")
axis(side = 1, labels = TRUE, col = "grey80", cex.axis = 0.75)

dev.off()
