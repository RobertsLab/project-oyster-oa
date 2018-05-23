#In this script, I will create a map of the five outplant locations.

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
#jpeg(filename = "images/DNR/2018-05-23-Site-Map/2018-05-23-Outplant-Site-Map.jpeg", height = 1000) #Create a new file to save map

data(nepacLLhigh) #Load set of polygons for the NE Pacific Ocean in high resolution from PBSmapping
plotMap(nepacLLhigh, xlim = c(-125, -121.9), ylim = c(46, 48.9), col = "snow1", bg = "slategray1", xaxt = "n", yaxt = "n", xlab = "", ylab = "", ann = FALSE) #Create a map with high resolution NE Pacific Ocean data. Remove axes since those will be manually added

#### MODIFY BASE MAP ####
axis(side = 1, at = c(-124.5, -124, -123.5, -123, -122.5, -122), labels=c("124.5°W", "124°W", "123.5°W", "123°W", "122.5°W", "122°W"), tick = TRUE, col.axis = "grey20") #Add longitude axis
axis(side = 2, at = c(46.5, 47, 47.5, 48, 48.5), labels=c("46.5ºN", "47°N", "47.5°N", "48°N", "48.5°N"), tick = TRUE, col.axis = "grey20") #Add latitude axis

textCoordinates <- locationCords #Duplicate dataframe
textCoordinates$Latitude[5] <- 48.4 #Move FB text down so it's easier to read
head(textCoordinates) #Confirm changes

symbols(x = locationCords$Longitude, y = locationCords$Latitude, circles = c(rep(1, times = 5)), add = TRUE, inches = 0.05, bg = "grey30") #Add points to map
text(x = textCoordinates$Longitude, y = textCoordinates$Latitude, labels = locationCords$Abbreviation, pos = c(2, 2, 2, 1, 4), col = "black", cex = 1.8) #Add site abbreviations as labels

#dev.off() #Turn off plotting device