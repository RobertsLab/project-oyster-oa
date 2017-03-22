#In this script, I'll do a preliminary comparison to determine which proteins are differentially expressed between sites and eelgrass conditions.

#Step 1: Import data
averageAreaAdjustedMerged <- read.csv(file = "/Users/yaaminivenkataraman/Documents/School/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Oyster-AverageAdjustedMergedArea.csv", header = TRUE)
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Removing the extra column "X"

#Step 2: Eelgrass vs. Bare

#Extract columns
proteins <- averageAreaAdjustedMerged$averageAreaAdjusted.proteins


#Find average of all bare and eelgrass sites

#Calculate ratios between eelgrass and bare sites
ratiosCaseInlet <- averageAreaAdjustedMerged$eelgrassCaseInlet/averageAreaAdjustedMerged$bareCaseInlet
ratiosFidalgoBay <- averageAreaAdjustedMerged$eelgrassFidalgoBay/averageAreaAdjustedMerged$bareFidalgoBay
ratiosWillapaBay <- averageAreaAdjustedMerged$eelgrassWillapaBay/averageAreaAdjustedMerged$bareWillapaBay
ratiosSkokomishRiver <- averageAreaAdjustedMerged$eelgrassSkokomishRiver/averageAreaAdjustedMerged$bareSkokomishRiver
ratiosPortGamble <- averageAreaAdjustedMerged$eelgrassPortGamble/averageAreaAdjustedMerged$barePortGamble

#Create dataframe with protein names and ratios
ratiosEelgrassBare <- data.frame(proteins, ratiosCaseInlet, ratiosFidalgoBay, ratiosWillapaBay, ratiosSkokomishRiver, ratiosPortGamble)

#Select rows 

#Step 3: Case Inlet vs. Fidalgo Bay vs. Willapa Bay vs. Skokomish River vs. Port Gamble
