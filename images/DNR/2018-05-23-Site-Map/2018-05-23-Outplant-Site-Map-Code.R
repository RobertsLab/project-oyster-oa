#In this script, I will create a map of the five outplant locations

#### LOAD PACKAGES AND DEPENDENCIES ####

rm(list = ls()) #Remove list

require(maps) #Basic mapping functions and some data
require(mapdata) #Some additional HiRes data
require(maptools) #Useful tools such as reading shapefiles
require(mapproj) #Various mapping projections
require(PBSmapping) #Powerful mapping functions developed by Pacific Biological Station

# =======================================================
# Creating Map of locations
data(nepacLLhigh)
locationCords <- read.csv("WSGpopLocations.csv")
head(locationCords)

colnames(locationCords) <- c("location", "lat", "long")
locationCords <- locationCords[order(locationCords$lat),]
head(locationCords)
locationCords

graphics.off()
nrow(locationCords)

broodLoc <- locationCords[c(1,4,5,7,8,9,10,11,13),]

#================================================================
# Map for OYSTER ECOLOGY
# all populations sampled for reproduction 



# windows(8,7)

#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-124.8, -121.5), ylim=c(46.9, 49.1), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE, axis=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="", 
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, col.axis='grey20')
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, col.axis='grey20')
# add points for each location
# symbol colors
sym.col <- c(rep("#00000088", 4), ##grey, most sites
             rep("#00EE00", 1), #green, Mud
             rep("#00000088", 4), #grey, most sites
             rep("#00EE00", 1), #green, Fidalgo
             rep("#00000088", 1)) # las site - lopez

symbols(locationCords$long, locationCords$lat,
        circles=c(rep(1, times=11)), add=T, inches=0.07,
        fg="grey10", bg=sym.col )
# add labels to points
text(x=locationCords$long, y=locationCords$lat,
     labels=locationCords$location,
     pos=c(2,4,2,2,4,2,4,2,2,1,4,2), col="black",
     cex=1.2)

# label Vancouver Island
text(x=-124.18, y=48.68, labels="Vancouver Island", col="grey60")

# label Washington State
text(x=-124, y=47.85, labels="Washington \n State", col="grey60")



# -------------------------------------------
#================================================================
# Map for OYSTER ECOLOGY & TRACE ELEMENTAL FINGERPRINTING 
# all populations sampled for brooding larvae

windows(8,7)


# pdf('Map_BroodLocationsSampled2.pdf', 8,7)


#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-124.2, -121.7), ylim=c(46.8, 48.9), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE, axis=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="", 
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, col.axis='grey20')
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, col.axis='grey20')





# add points for each location
# symbol colors
sym.col <- c(rep("#00000088", 13)) ##grey, most sites


symbols(locationCords$long, locationCords$lat,
        circles=c(rep(1, times=13)), add=T, inches=0.07,
        fg="grey10", bg=sym.col )
# add labels to points
text(x=locationCords$long, y=locationCords$lat,
     labels=locationCords$location,
     pos=c(2,4,2,2,4,2,4,2,2,1,4,2, 3), col="black",
     cex=1.2)

# label Vancouver Island
text(x=-123.9, y=48.68, labels="Vancouver \n Island", col="grey60")

# label Washington State
text(x=-123.8, y=47.85, labels="Washington \n State", col="grey60")

map.scale(relwidth = .2 ,ratio=FALSE)


# seattle
# symbols(-122.328831, 47.600731, squares=1, add=T, 
#         inches=0.1, fg='darkblue', bg='blue')
# 
# text(-122.328831, 47.600731, labels='Seattle', 
#      pos=4, col='black', cex=1.2, offset=1.2)

dev.off()

###############################################################
################################################################

#================================================================
# Map for OYSTER ECOLOGY & TRACE ELEMENTAL FINGERPRINTING 
# all populations sampled for brooding larvae
# populations where late stage larvae were collected

# windows(8,7)


# pdf('Map_BroodLocationsSamped+Collection.pdf', 8,7)


#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-124.2, -121.7), ylim=c(46.8, 48.9), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE, axis=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="", 
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, col.axis='grey20')
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, col.axis='grey20')

# add points for each location
# symbol colors
sym.col <- c(rep("#00000088", 13)) ##grey, most sites


symbols(locationCords$long, locationCords$lat,
        circles=c(rep(1, times=13)), add=T, inches=0.07,
        fg="grey10", bg=sym.col )
# add labels to points
text(x=locationCords$long, y=locationCords$lat,
     labels=locationCords$location,
     pos=c(2,4,2,2,4,2,4,2,2,1,4,2, 3), col="black",
     cex=1.2)

# highlight sampled populations
symbols(broodLoc$long, broodLoc$lat,
        circles=c(rep(1, times=9)), add=T, inches=0.07,
        fg="grey10", bg='gold2' )

# label Vancouver Island
text(x=-123.9, y=48.68, labels="Vancouver \n Island", col="grey60")

# label Washington State
text(x=-123.8, y=47.85, labels="Washington \n State", col="grey60")

map.scale(relwidth = .2 ,ratio=FALSE)



# dev.off()


# dev.off()




#================================================================
# Map for TRACE ELEMENTAL FINGERPRINTING 
# all populations brooding larvae were collected
# ************* ONLY COLLECTED 

#windows(8,7)

# pdf('Map_BroodCollection.pdf', 8,7)


#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-124.2, -121.7), ylim=c(46.8, 48.9), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="", 
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, col.axis='grey20')
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, col.axis='grey20')

# add points for each location
# symbol colors
sym.col <- c(rep("#00000088", 9)) ##grey, most sites


symbols(broodLoc$long, broodLoc$lat,
        circles=c(rep(1, times=9)), add=T, inches=0.05,
        fg="grey10", bg=sym.col )

# add labels to points
text(x=broodLoc$long, y=broodLoc$lat,
     labels=broodLoc$location,
     pos=c(2,2,4,4,2,2,1,4,3), col="black",
     cex=1)

# label Vancouver Island
text(x=-123.9, y=48.68, labels="Vancouver \n Island", col="grey60")

# label Washington State
text(x=-123.8, y=47.85, labels="Washington \n State", col="grey60")

map.scale(relwidth = .2 ,ratio=FALSE)



dev.off()


#================================================================
# Map for TRACE ELEMENTAL FINGERPRINTING 
# all populations brooding larvae were collected
# ************* ONLY COLLECTED 
# PLOT 2 - no names

#windows(8,7)


# pdf('Map_BroodCollection_nolabels.pdf', 8,7)

# remove LSI
broodLoc
broodLoc2 <- broodLoc[-c(broodLoc$location=='Little Skookum'), ]


#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-124.2, -121.7), ylim=c(46.8, 48.9), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="", 
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, col.axis='grey20')
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, col.axis='grey20')

# add points for each location
# symbol colors
sym.col <- c(rep("#00000088", 8)) ##grey, most sites


symbols(broodLoc2$long, broodLoc2$lat,
        circles=c(rep(1, times=8)), add=T, inches=0.05,
        fg="grey10", bg=sym.col )


# label Vancouver Island
text(x=-123.9, y=48.68, labels="Vancouver \n Island", col="grey60")

# label Washington State
text(x=-123.8, y=47.85, labels="Washington \n State", col="grey60")

map.scale(relwidth = .2 ,ratio=FALSE)



dev.off()

#####################################################################
#################################################################


# WEST COAST BEST COAST
#west coast
# plotMap(nepacLLhigh, xlim=c(-130, -110), ylim=c(30, 55), 
#         border="lightskyblue2", bg='lightskyblue1', #ann=FALSE, axis=FALSE,
#         #xaxt="n", yaxt="n", ylab="", xlab="", 
#         col = "snow1")
# map("state", region = c("washington", "oregon", "california"),
#     xlim=c(-125, -117), ylim=c(32, 50))



#==================================================
# PUGET SOUND MAP

windows(5,5)
#Puget Sound limits
plotMap(nepacLLhigh, xlim=c(-125, -121), ylim=c(46.8, 49.1), 
        border="lightskyblue2", bg='lightskyblue1', ann=FALSE, axis=FALSE,
        xaxt="n", yaxt="n", ylab="", xlab="",
        col = "snow1")
# adding axis labels
axis(side=1, at=c(-124,-123,-122), 
     labels=c("124°W", "123°W", "122°W") , 
     tick=TRUE, cex.axis=.8)
axis(side=2, at=c(47,47.5,48,48.5, 49), 
     labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
     tick=TRUE, cex.axis=.8)

# add points for the two location
mainL <- locationCords[c(5,10), ]

# seattle
symbols(-122.328831, 47.600731, squares=1, add=T, 
        inches=0.1, fg='black', bg='black')
text(-122.328831, 47.600731, labels='Seattle', 
     pos=4, col='black', cex=1.2, offset=.7)

# # add symbols
# symbols(mainL$long, mainL$lat,
#         circles=c(rep(1, times=2)), add=T, inches=0.05,
#         fg="grey10", bg="grey10")
# # add labels to points
# text(x=mainL$long, y=mainL$lat,
#      labels=c("Dyes Inlet", "Fidalgo Bay"),
#      pos=c(2,4), col="black",
#      cex=1.2)

# label Vancouver Island
text(x=-124.26, y=48.75, labels="Vancouver Island", col="grey60", cex=1)

# label Washington State
text(x=-123.8, y=47.78, labels="Washington \n State", col="grey60", cex=1)

# add compass
# addCompass(-121.6,48.65, cex=.8)




# 
# 
# #=================================================================
# # Color coding puget sound regions (North, central, & south)
# 
# #pdf("MapofPopEMILY.PDF", 8,8)
# #windows(8,8)
# #Puget Sound limits
# plotMap(nepacLLhigh, xlim=c(-125, -121.8), ylim=c(47, 48.8), 
#         border="grey", bg='lightblue', col = "seashell1")
# # add points for each location
# # symbol colors
# sym.col <- c(rep("#FF000088", 3), #red, south sound
#              rep("#FFFF0088", 5), #yellow, central sound
#              rep("#008B0088", 3)) #green, north sound
# 
# symbols(locationCords$long, locationCords$lat, 
#         circles=c(rep(1, times=11)), add=T, inches=0.06,
#         fg="grey10", bg=sym.col )
# # add labels to points
# text(x=locationCords$long, y=locationCords$lat, 
#      labels=locationCords$location,
#      pos=c(2,4,2,2,4,2,4,2,2,1,2), col="black",
#      cex=.8)
# 
# # label Vancouver Island
# text(x=-124.18, y=48.68, labels="Vancouver Island", col="grey60")
# 
# # label Washington State
# text(x=-124, y=47.85, labels="Washington \n State", col="grey60")
# 
# dev.off()
# 
# #legend
# 
# #=================================================================
# # Map of main locations
# # Dyes Inlet and Fidalgo Bay
# 
# 
# #windows(5,6)
# #Puget Sound limits
# plotMap(nepacLLhigh, xlim=c(-125, -121.5), ylim=c(46.9, 49.1), 
#         border="snow", bg='lightskyblue1', ann=FALSE, axis=FALSE,
#         xaxt="n", yaxt="n", ylab="", xlab="",
#         col = "beige")
# # adding axis labels
# axis(side=1, at=c(-124,-123,-122), 
#      labels=c("124°W", "123°W", "122°W") , 
#      tick=TRUE)
# axis(side=2, at=c(47,47.5,48,48.5, 49), 
#      labels=c("47°N", "47.5°N", "48°N", "48.5°N", "49°N") , 
#      tick=TRUE)
# 
# # add points for the two location
# mainL <- locationCords[c(5,11), ]
# 
# # other locations
# otherlocal <- locationCords[-c(5,11),]
# nrow(otherlocal)
# 
# # other locations
# symbols(otherlocal$long, otherlocal$lat,
#         circles=c(rep(1, times=9)), add=T, inches=0.05,
#         fg="grey65", bg="grey65")
# 
# # add symbols
# symbols(mainL$long, mainL$lat,
#         circles=c(rep(1, times=2)), add=T, inches=0.07,
#         fg="grey10", bg="grey10")
# # add labels to points
# text(x=mainL$long, y=mainL$lat,
#      labels=c("Dyes Inlet", "Fidalgo Bay"),
#      pos=c(4,2), col="black",
#      cex=1.5)
# 
# # label Vancouver Island
# text(x=-124.18, y=48.75, labels="Vancouver Island", col="grey50", cex=1.2)
# 
# # label Washington State
# text(x=-124, y=47.85, labels="Washington \n State", col="grey50", cex=1.2)
# 
# # END
# 
# 
# 
# # 
# # #Puget Sound limits
# # plotMap(nepacLLhigh, xlim=c(-125, -121.5), ylim=c(46.9, 49), 
# #         border="grey", bg='lightblue', col = "seashell1")
# 
# 
# # trying map colors
# # grey90
# # "cornsilk2"
# # cornsilk1
# # seashell
# 
# # green #00FF00 
# # dark green #008B00 
# # yellow 1 #FFFF00 
# 
