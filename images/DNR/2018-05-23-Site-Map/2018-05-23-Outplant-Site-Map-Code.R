#In this script, I will create a map of the five outplant locations

#### SET WORKING DIRECTORY ####
getwd()
setwd("../../../") #Set directory to top level folder of repository

#### LOAD PACKAGES AND DEPENDENCIES ####

rm(list = ls()) #Remove list

#Install packages
install.packages("maps") #Basic mapping functions and some data
install.packages("mapdata") #Some additional HiRes data
install.packages("maptools") #Useful tools such as reading shapefiles
install.packages("mapproj") #Various mapping projections
install.packages("PBSmapping") #Powerful mapping functions developed by Pacific Biological Station

#Load packages
library(maps)
library(mapdata) 
library(maptools)
library(mapproj)
library(PBSmapping)

#### IMPORT LATITUDE AND LONGITUDE INFORMATION ####
locationCords <- read.csv("data/DNR/2018-05-23-Outplant-Coordinates.csv") #Import outplant coordinate information
head(locationCords) #Confirm import

locationCords <- locationCords[,1:4] #Only keep columns with data
locationCords <- locationCords[order(locationCords$Latitude),] #Reorder location coordinates by latittude
head(locationCords) #Confirm changes

#### CREATE BASE MAP ####
data(nepacLLhigh) #Load set of polygons for the NE Pacific Ocean in high resolution from PBSmapping
plotMap(nepacLLhigh, xlim = c(-125, -121.9), ylim = c(46, 48.9), col = "snow1", bg = "lightskyblue1", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ann = FALSE, axis = FALSE) #Create a map with high resolution NE Pacific Ocean data. Remove axes since those will be manually added

#### MODIFY BASE MAP ####


#================================================================

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