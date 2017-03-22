#In this script, I'll do a preliminary comparison to determine which proteins are differentially expressed between sites and eelgrass conditions.

#Step 1: Import data
averageAreaAdjustedMerged <- read.csv(file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Skyline_20170314/Oyster-AverageAdjustedMergedArea.csv", header = TRUE)
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Removing the extra column "X"

#Step 2: Eelgrass vs. Bare

#Extract columns
proteins <- averageAreaAdjustedMerged$averageAreaAdjusted.proteins

bareCaseInlet <- averageAreaAdjustedMerged$bareCaseInlet
bareFidalgoBay <- averageAreaAdjustedMerged$bareFidalgoBay
bareWillapaBay <- averageAreaAdjustedMerged$bareWillapaBay
bareSkokomisRiver <- averageAreaAdjustedMerged$bareSkokomishRiver
barePortGamble <- averageAreaAdjustedMerged$barePortGamble

eelgrassCaseInlet <- averageAreaAdjustedMerged$eelgrassCaseInlet
eelgrassFidalgoBay <- averageAreaAdjustedMerged$eelgrassFidalgoBay
eelgrassWillapaBay <- averageAreaAdjustedMerged$eelgrassWillapaBay
eelgrassSkokomishRiver <- averageAreaAdjustedMerged$eelgrassSkokomishRiver
eelgrassPortGamble <- averageAreaAdjustedMerged$eelgrassPortGamble

#Find average of all bare and eelgrass sites
averageBare <- ave(, )
averageEelgrass

#Calculate ratios between eelgrass and bare sites

#Create dataframe with protein names and ratios
ratiosEelgrassBare <- data.frame(proteins, ratiosCaseInlet, ratiosFidalgoBay, ratiosWillapaBay, ratiosSkokomishRiver, ratiosPortGamble)

#Select rows 

#Step 3: Case Inlet vs. Fidalgo Bay vs. Willapa Bay vs. Skokomish River vs. Port Gamble
