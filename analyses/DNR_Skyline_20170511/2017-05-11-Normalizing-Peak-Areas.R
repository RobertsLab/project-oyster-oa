#### Normalizing Peak Areas with TIC content ####

peakAreas <- read.csv("2017-05-11-peak-areas.csv", na.strings = "#N/A") #Import peak area data
TICvalues <- data.frame(oyster = seq(from = 1, to = 25, by = 1),
                        TIC = c(1.52E+09, 1.48E+09, NA, NA, 1.48E+09, 2.67E+09, 1.53E+09, 1.71E+09, 1.75E+09, 2.70E+09, 1.65E+09, 1.63E+09, 1.54E+09, 2.63E+09, 1.52E+09, 2.22E+09, 1.56E+09, 1.88E+09, 1.37E+09, 2.47E+09, 6.91E+08, 2.26E+08, 5.98E+07, 5.66E+08, 6.12E+09)) #Manually create dataframe for TIC data

#Delete extraneous rows from TIC value spreadsheet
TICvalues <- TICvalues[-c(17, 21, 22, 23),] #Delete the rows with TIC values for samples not included in Skyline report
TICvalues #Confirm changes

#Divide the Peak Area values with the appropriate TIC values
