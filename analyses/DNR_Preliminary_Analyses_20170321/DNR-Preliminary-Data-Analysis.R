#In this script, I'll do a preliminary comparison to determine which proteins are differentially expressed between sites and eelgrass conditions.

#Step 1: Import data

averageAreaAdjustedMerged <- read.csv(file = "/Users/yaaminivenkataraman/Documents/School/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Oyster-AverageAdjustedMergedArea.csv", header = TRUE)
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Removing the extra column "X"

#Step 2: Eelgrass vs. Bare

#I didn't do any type of robust statistical analysis (MSstats or otherwise), I just did a few simple pivot tables showing sum area, max area, & average area, then calculated ratios of the above stats with eelgrass/bare and then pulled any protein that had a ratio > 1.5 (more expressed in eelgrass). Again, not robust, but with my limited time and limited data that's good enough for a poster! 

#Step 3: Case Inlet vs. Fidalgo Bay vs. Willapa Bay vs. Skokomish River vs. Port Gamble
