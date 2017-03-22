#In this script, I will reformat my data table so I can use it in MSStats.

#Step 1: Import data
AverageArea <- read.csv(file = "/Users/yaaminivenkataraman/Documents/School/project-oyster-oa/analyses/DNR_Skyline_20170314/Oyster-AverageArea-Proteinbased.csv", header = TRUE)
AverageArea #View data

#Step 2: Remove columns from poor O107 runs
averageAreaAdjusted <- within(AverageArea, rm("X2017_January_23_envtstress_oyster2.raw", "X2017_January_23_envtstress_oyster17.raw", "X2017_January_23_envtstress_oyster21.raw", "X2017_January_23_envtstress_oyster22.raw", "X2017_January_23_envtstress_oyster23.raw")) #Based on my sequence file, I need to remove samples oyster 2, 17, 22 and 23. I'm also going to remove the blank.
averageAreaAdjusted #View data to confirm changes

#Step 3: Change column names
names(averageAreaAdjusted) <- c("proteins", "bare-willapa-bay-1", "bare-skokomish-river-1", "bare-fidalgo-bay-1", "bare-willapa-bay-2", "eelgrass-fidalgo-bay-1", "bare-port-gamble-1", "bare-case-inlet-1", "bare-skokomish-river-2", "eelgrass-willapa-bay-1", "eelgrass-case-inlet-1", "eelgrass-port-gamble-1", "eelgrass-skokomish-river-1", "eelgrass-skokomish-river-2", "eelgrass-case-inlet-2", "eelgrass-port-gamble-2", "bare-fidalgo-bay-2", "bare-port-gamble-2", "eelgrass-fidalgo-bay-2", "bare-case-inlet-2", "eelgrass-willapa-bay-2") #Change the column names to indicate condition, location, and replicate
averageAreaAdjusted #View data to ensure the column names changed

#Step 4: Merge columns together. I want to average the values between my replicates so I have a total of 10 columns with data to analyze.

#First, I need to extract all of the columns
bareCaseInlet1 <- averageAreaAdjusted$`bare-case-inlet-1`
bareCaseInlet2 <- averageAreaAdjusted$`bare-case-inlet-2`
bareFidalgoBay1 <- averageAreaAdjusted$`bare-fidalgo-bay-1`
bareFidalgoBay2 <- averageAreaAdjusted$`bare-fidalgo-bay-1`
bareWillapaBay1 <- averageAreaAdjusted$`bare-willapa-bay-1`
bareWillapaBay2 <- averageAreaAdjusted$`bare-willapa-bay-2`
bareSkokomishRiver1 <- averageAreaAdjusted$`bare-skokomish-river-1`
bareSkokomishRiver2 <- averageAreaAdjusted$`bare-skokomish-river-2`
barePortGamble1 <- averageAreaAdjusted$`bare-port-gamble-1`
barePortGamble2 <- averageAreaAdjusted$`bare-port-gamble-2`

eelgrassCaseInlet1 <- averageAreaAdjusted$`eelgrass-case-inlet-1`
eelgrassCaseInlet2 <- averageAreaAdjusted$`eelgrass-case-inlet-2`
eelgrassFidalgoBay1 <- averageAreaAdjusted$`eelgrass-fidalgo-bay-1`
eelgrassFidalgoBay2 <- averageAreaAdjusted$`eelgrass-fidalgo-bay-2`
eelgrassWillapaBay1 <- averageAreaAdjusted$`eelgrass-willapa-bay-1`
eelgrassWillapaBay2 <- averageAreaAdjusted$`eelgrass-willapa-bay-2`
eelgrassSkokomishRiver1 <- averageAreaAdjusted$`eelgrass-skokomish-river-1`
eelgrassSkokomishRiver2 <- averageAreaAdjusted$`eelgrass-skokomish-river-2`
eelgrassPortGamble1 <- averageAreaAdjusted$`eelgrass-port-gamble-1`
eelgrassPortGamble2 <- averageAreaAdjusted$`eelgrass-port-gamble-2`

#Then I will average all of the pairs

bareCaseInlet <- ave(bareCaseInlet1, bareCaseInlet2)
bareFidalgoBay <- ave(bareFidalgoBay1, bareFidalgoBay2)
bareWillapaBay <- ave(bareWillapaBay1, bareWillapaBay2)
bareSkokomishRiver <- ave(bareSkokomishRiver1, bareSkokomishRiver2)
barePortGamble <- ave(barePortGamble1, barePortGamble2)

eelgrassCaseInlet <- ave(eelgrassCaseInlet1, eelgrassCaseInlet2)
eelgrassFidalgoBay <- ave(eelgrassFidalgoBay1, eelgrassFidalgoBay2)
eelgrassWillapaBay <- ave(eelgrassWillapaBay1, eelgrassWillapaBay2)
eelgrassSkokomishRiver <- ave(eelgrassSkokomishRiver1, eelgrassSkokomishRiver2)
eelgrassPortGamble <- ave(eelgrassPortGamble1, eelgrassPortGamble2)

#And finally I'll create a new data frame

averageAreaAdjustedMerged <- data.frame(averageAreaAdjusted$proteins, bareCaseInlet, bareFidalgoBay, bareWillapaBay, bareSkokomishRiver, barePortGamble, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassWillapaBay, eelgrassSkokomishRiver, eelgrassPortGamble)

#Step 5: Export the new .csv

write.csv(x = averageAreaAdjustedMerged, file = "/Users/yaaminivenkataraman/Documents/School/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Oyster-AverageAdjustedMergedArea.csv")


######


#Now, I'll do the same with my maximum area spreadsheet

#Step 1: Import data
MaxArea <- read.csv(file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Skyline_20170314/Oyster-MaxArea-Proteinbased.csv", header = TRUE)
MaxArea #View data

#Step 2: Remove columns from poor O107 runs
maxAreaAdjusted <- within(MaxArea, rm("X2017_January_23_envtstress_oyster2.raw", "X2017_January_23_envtstress_oyster17.raw", "X2017_January_23_envtstress_oyster21.raw", "X2017_January_23_envtstress_oyster22.raw", "X2017_January_23_envtstress_oyster23.raw")) #Based on my sequence file, I need to remove samples oyster 2, 17, 22 and 23. I'm also going to remove the blank.
maxAreaAdjusted #View data to confirm changes

#Step 3: Change column names
names(maxAreaAdjusted) <- c("proteins", "bare-willapa-bay-1", "bare-skokomish-river-1", "bare-fidalgo-bay-1", "bare-willapa-bay-2", "eelgrass-fidalgo-bay-1", "bare-port-gamble-1", "bare-case-inlet-1", "bare-skokomish-river-2", "eelgrass-willapa-bay-1", "eelgrass-case-inlet-1", "eelgrass-port-gamble-1", "eelgrass-skokomish-river-1", "eelgrass-skokomish-river-2", "eelgrass-case-inlet-2", "eelgrass-port-gamble-2", "bare-fidalgo-bay-2", "bare-port-gamble-2", "eelgrass-fidalgo-bay-2", "bare-case-inlet-2", "eelgrass-willapa-bay-2") #Change the column names to indicate condition, location, and replicate
maxAreaAdjusted #View data to ensure the column names changed

#Step 4: Merge columns together. I want to average the values between my replicates so I have a total of 10 columns with data to analyze.

#First, I need to extract all of the columns
bareCaseInlet1 <- maxAreaAdjusted$`bare-case-inlet-1`
bareCaseInlet2 <- maxAreaAdjusted$`bare-case-inlet-2`
bareFidalgoBay1 <- maxAreaAdjusted$`bare-fidalgo-bay-1`
bareFidalgoBay2 <- maxAreaAdjusted$`bare-fidalgo-bay-1`
bareWillapaBay1 <- maxAreaAdjusted$`bare-willapa-bay-1`
bareWillapaBay2 <- maxAreaAdjusted$`bare-willapa-bay-2`
bareSkokomishRiver1 <- maxAreaAdjusted$`bare-skokomish-river-1`
bareSkokomishRiver2 <- maxAreaAdjusted$`bare-skokomish-river-2`
barePortGamble1 <- maxAreaAdjusted$`bare-port-gamble-1`
barePortGamble2 <- maxAreaAdjusted$`bare-port-gamble-2`

eelgrassCaseInlet1 <- maxAreaAdjusted$`eelgrass-case-inlet-1`
eelgrassCaseInlet2 <- maxAreaAdjusted$`eelgrass-case-inlet-2`
eelgrassFidalgoBay1 <- maxAreaAdjusted$`eelgrass-fidalgo-bay-1`
eelgrassFidalgoBay2 <- maxAreaAdjusted$`eelgrass-fidalgo-bay-2`
eelgrassWillapaBay1 <- maxAreaAdjusted$`eelgrass-willapa-bay-1`
eelgrassWillapaBay2 <- maxAreaAdjusted$`eelgrass-willapa-bay-2`
eelgrassSkokomishRiver1 <- maxAreaAdjusted$`eelgrass-skokomish-river-1`
eelgrassSkokomishRiver2 <- maxAreaAdjusted$`eelgrass-skokomish-river-2`
eelgrassPortGamble1 <- maxAreaAdjusted$`eelgrass-port-gamble-1`
eelgrassPortGamble2 <- maxAreaAdjusted$`eelgrass-port-gamble-2`

#Then I will average all of the pairs

bareCaseInlet <- ave(bareCaseInlet1, bareCaseInlet2)
bareFidalgoBay <- ave(bareFidalgoBay1, bareFidalgoBay2)
bareWillapaBay <- ave(bareWillapaBay1, bareWillapaBay2)
bareSkokomishRiver <- ave(bareSkokomishRiver1, bareSkokomishRiver2)
barePortGamble <- ave(barePortGamble1, barePortGamble2)

eelgrassCaseInlet <- ave(eelgrassCaseInlet1, eelgrassCaseInlet2)
eelgrassFidalgoBay <- ave(eelgrassFidalgoBay1, eelgrassFidalgoBay2)
eelgrassWillapaBay <- ave(eelgrassWillapaBay1, eelgrassWillapaBay2)
eelgrassSkokomishRiver <- ave(eelgrassSkokomishRiver1, eelgrassSkokomishRiver2)
eelgrassPortGamble <- ave(eelgrassPortGamble1, eelgrassPortGamble2)

#And finally I'll create a new data frame

maxAreaAdjustedMerged <- data.frame(maxAreaAdjusted$proteins, bareCaseInlet, bareFidalgoBay, bareWillapaBay, bareSkokomishRiver, barePortGamble, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassWillapaBay, eelgrassSkokomishRiver, eelgrassPortGamble)

#Step 5: Export the new .csv

write.csv(x = maxAreaAdjustedMerged, file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Oyster-MaxAdjustedMergedArea.csv")