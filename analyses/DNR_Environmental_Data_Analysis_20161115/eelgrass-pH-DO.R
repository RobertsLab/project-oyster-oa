#In this script, I will create a diagram explaining how pH and dissolved oxygen cycle with eelgrass photosynthesis and respiration.

pdf(file = "Eelgrass-pH-DO.pdf", width = 9, height = 4) #create .pdf to save graph

curve(expr = sin(((pi*x)/24)), from = 0, to = 24, xlab = "Time (Hours)", ylab = "", xaxt = "n", yaxt = "n", col = "skyblue", lwd = 3) #create base curve
axis(side = 1, at = seq(0, 24, 2), labels = seq(0, 24, 2)) #add x-axis tick marks
mtext(text = "Dissolved Oxygen (mg/L)", side = 2, line = 1) #add y-axis label

dev.off() #close pdf file
