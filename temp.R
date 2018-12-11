plotMatrix <- matrix(c(1, 2,
                       1, 3), nrow = 2, ncol = 2, byrow = TRUE) #Create a matrix and fill it in by row.
plotMatrix #Confirm matrix creation

pdf("analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-12-01-Multipanel-Ordination.pdf", width = 11, height = 8.5)

par(mar = c(0, 0, 0, 2.5), oma = c(3, 3.5, 1, 0)) #Specify inner and outer margins
layout(mat = plotMatrix, width = c(45, 25)) #Create a layout based on the plot matrix. Column 1s width should be 3x as large as column 2. 
#layout.show(n = 3) #Confirm the plot layout is what is desired

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
mtext(side = 3, line = -8, at = c(-1, -9), text = "     a. RDA") #Add test name

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
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -4, adj = -1, text = "                 b. NMDS: Stress = 0.075") #Add stress value

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
axis(side = 2, labels = TRUE, col = "grey80", cex.axis = 0.75)
mtext(side = 3, line = -4, adj = -1, text = "                    c. Influential peptides") #Add subplot description

#### LEGEND ####

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE) #Solution from KLo's blog! http://dr-k-lo.blogspot.com/2014/03/the-simplest-way-to-plot-legend-outside.html. Overlay a larger plot onto the already-created plot
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n") #Create an empty plot

rect(xleft = -0.80, ybottom = 0.7782, xright = 0.10, ytop = 1.1, col = "white", border = "grey80") #Add a box with a grey80 border to section off legend. The top of the box will bleed off the page.

#Legend with PS site specification
legend(x = -0.78, y = 1.05, xpd = TRUE, inset = c(0, 0),
       legend = c("Case Inlet", "Fidalgo Bay", "Port Gamble Bay", "Skokomish River Delta"),
       pch = c(rep(16, times = 4)), 
       col = c('#00A9BD', '#38001C', '#440D82', '#017A74'),
       cex = rep(1.2, times = 4),
       bg = "white", box.col = "white") #Create a horizontal legend (horiz = TRUE) that can be plotted outside of the plot boundaries (xpd = TRUE). Place the legend at x = -1, y = 1.

#Legend with WB, peptide and habitat specification
legend(x = -0.30, y = 1.05, xpd = TRUE, inset = c(0, 0),
       legend = c("Willapa Bay", "Unvegetated", "Eelgrass", "Influential Peptides"),
       pch = c(16, 1, 16, 4), 
       col = c('#EB8B0C', "grey20", "grey20", "grey20"),
       cex = rep(1.2, times = 4),
       bg = "white", box.col = "white") #Create a horizontal legend (horiz = TRUE) that can be plotted outside of the plot boundaries (xpd = TRUE). Place the legend at x = -1, y = 1.

dev.off()