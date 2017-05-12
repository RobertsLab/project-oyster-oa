#### Normalizing Peak Areas with TIC content ####

peakAreas <- read.csv("2017-05-11-peak-areas.csv", na.strings = "#N/A") #Import peak area data
TICvalues <- data.frame(oyster = seq(from = 1, to = 25, by = 1),
                        TIC = c(1.52E+09, 1.48E+09, NA, NA, 1.48E+09, 2.67E+09, 1.53E+09, 1.71E+09, 1.75E+09, 2.70E+09, 1.65E+09, 1.63E+09, 1.54E+09, 2.63E+09, 1.52E+09, 2.22E+09, 1.56E+09, 1.88E+09, 1.37E+09, 2.47E+09, 6.91E+08, 2.26E+08, 5.98E+07, 5.66E+08, 6.12E+09)) #Manually create dataframe for TIC data

#Delete extraneous rows from TIC value spreadsheet
TICvalues <- TICvalues[-c(17, 21, 22, 23),] #Delete the rows with TIC values for samples not included in Skyline report
TICvalues #Confirm changes

#Divide the Peak Area values with the appropriate TIC values
peakAreas[,4] <- peakAreas$X1.Area/TICvalues$oyster[1]
peakAreas[,5] <- peakAreas$X1.Area/TICvalues$oyster[2]
peakAreas[,6] <- peakAreas$X1.Area/TICvalues$oyster[3]
peakAreas[,7] <- peakAreas$X1.Area/TICvalues$oyster[4]
peakAreas[,8] <- peakAreas$X1.Area/TICvalues$oyster[5]
peakAreas[,9] <- peakAreas$X1.Area/TICvalues$oyster[6]
peakAreas[,10] <- peakAreas$X1.Area/TICvalues$oyster[7]
peakAreas[,11] <- peakAreas$X1.Area/TICvalues$oyster[8]
peakAreas[,12] <- peakAreas$X1.Area/TICvalues$oyster[9]
peakAreas[,13] <- peakAreas$X1.Area/TICvalues$oyster[10]
peakAreas[,14] <- peakAreas$X1.Area/TICvalues$oyster[11]
peakAreas[,15] <- peakAreas$X1.Area/TICvalues$oyster[12]
peakAreas[,16] <- peakAreas$X1.Area/TICvalues$oyster[13]
peakAreas[,17] <- peakAreas$X1.Area/TICvalues$oyster[14]
peakAreas[,18] <- peakAreas$X1.Area/TICvalues$oyster[15]
peakAreas[,19] <- peakAreas$X1.Area/TICvalues$oyster[16]
peakAreas[,20] <- peakAreas$X1.Area/TICvalues$oyster[17]
peakAreas[,21] <- peakAreas$X1.Area/TICvalues$oyster[18]
peakAreas[,22] <- peakAreas$X1.Area/TICvalues$oyster[19]
peakAreas[,23] <- peakAreas$X1.Area/TICvalues$oyster[20]
peakAreas[,24] <- peakAreas$X1.Area/TICvalues$oyster[21]

#Rename column headers to include DNR vial number, eelgrass condition and sample site.

names(peakAreas) <- c("protein-name", "peptide", "peptide-sequence", "O127-wb-bare-1",
"O107-sk-eelgrass-1", "O07-ci-eelgrass-1", "O77-pg-eelgrass-1", "O47-fb-bare-1", "O55-pg-bare-1", "O37-fb-eelgrass-1", "O15-ci-eelgrass-1", "O142-wb-eelgrass-1", "O119-sk-bare-1", "O47-fb-bare-2", "O127-wb-bare-2", "O37-fb-eelgrass-2", "O55-pg-bare-2", "O15-ci-bare-2", "O119-sk-bare-2", "O142-wb-eelgrass-2", "O07-ci-eelgrass-2", "O77-pg-eelgrass-2", "O107-2-sk-eegrass-1", "O107-2-sk-eelgrass-2")
head(peakAreas) #Confirm changes

#Save peakAreas as a new CSV file
write.csv(x = peakAreas, file = "2017-05-11-normalized-peak-areas.csv")

#### Renaming Column Headers without Normalizing ####

peakAreas <- read.csv("2017-05-11-peak-areas.csv", na.strings = "#N/A") #Import peak area data

#Rename column headers to include DNR vial number, eelgrass condition and sample site.

names(peakAreas) <- c("protein-name", "peptide", "peptide-sequence", "O127-wb-bare-1",
                      "O107-sk-eelgrass-1", "O07-ci-eelgrass-1", "O77-pg-eelgrass-1", "O47-fb-bare-1", "O55-pg-bare-1", "O37-fb-eelgrass-1", "O15-ci-bare-1", "O142-wb-eelgrass-1", "O119-sk-bare-1", "O47-fb-bare-2", "O127-wb-bare-2", "O37-fb-eelgrass-2", "O55-pg-bare-2", "O15-ci-bare-2", "O119-sk-bare-2", "O142-wb-eelgrass-2", "O07-ci-eelgrass-2", "O77-pg-eelgrass-2", "O107-2-sk-eegrass-1", "O107-2-sk-eelgrass-2")
head(peakAreas) #Confirm changes

#Save peakAreas as a new CSV file
write.csv(x = peakAreas, file = "2017-05-11-nonnormalized-peak-areas.csv")