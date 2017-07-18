#### IMPORT FULL PROTEIN AREA DATA ####

proteinAreas <- read.csv("2017-06-10-protein-areas-only-error-checked.csv", na.strings = "#N/A") #Specify Skyline's special way of designating N/A values
head(proteinAreas) #Confirm successful import
averageProteinAreas <- aggregate(proteinAreas[-1], proteinAreas[1], mean, na.action = na.omit, na.rm = TRUE) #Average protein areas across transitions and peptides
head(averageProteinAreas) #Confirm changes 

#### REFORMAT FULL PROTEIN DATA ####

names(averageProteinAreas) <- c("protein", "wb-bare-1", "sk-eelgrass-1", "ci-eelgrass-1", "pg-eelgrass-1", "fb-bare-1", "pg-bare-1", "fb-eelgrass-1", "ci-bare-1", "wb-eelgrass-1", "sk-bare-1", "fb-bare-2", "wb-bare-2", "fb-eelgrass-2", "pg-bare-2", "ci-bare-2", "sk-bare-2", "wb-eelgrass-2", "ci-eelgrass-2", "pg-eelgrass-2", "sk-eelgrass-2", "sk-eelgrass-3")
head(averageProteinAreas) #Confirm changes

#### AVERAGE REPLICATES ####

bareCaseInlet <- ave(averageProteinAreas$`ci-bare-1`, averageProteinAreas$`ci-bare-2`)
bareFidalgoBay <- ave(averageProteinAreas$`fb-bare-1`, averageProteinAreas$`fb-bare-2`)
bareWillapaBay <- ave(averageProteinAreas$`wb-bare-1`, averageProteinAreas$`wb-bare-2`)
bareSkokomishRiver <- ave(averageProteinAreas$`sk-bare-1`, averageProteinAreas$`sk-bare-2`)
barePortGamble <- ave(averageProteinAreas$`pg-bare-1`, averageProteinAreas$`pg-bare-2`)

eelgrassCaseInlet <- ave(averageProteinAreas$`ci-eelgrass-1`, averageProteinAreas$`ci-eelgrass-2`)
eelgrassFidalgoBay <- ave(averageProteinAreas$`fb-eelgrass-1`, averageProteinAreas$`fb-eelgrass-2`)
eelgrassWillapaBay <- ave(averageProteinAreas$`wb-eelgrass-1`, averageProteinAreas$`wb-eelgrass-2`)
eelgrassSkokomishRiver <- ave(averageProteinAreas$`sk-eelgrass-2`, averageProteinAreas$`sk-eelgrass-3`)
eelgrassPortGamble <- ave(averageProteinAreas$`pg-eelgrass-1`, averageProteinAreas$`pg-eelgrass-2`)

averageProteinAreasMerged <- data.frame(averageProteinAreas$protein, bareCaseInlet, bareFidalgoBay, barePortGamble, bareSkokomishRiver, bareWillapaBay, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassPortGamble, eelgrassSkokomishRiver, eelgrassWillapaBay) #Create new dataframe
names(averageProteinAreasMerged) <- c("protein", "bareCaseInlet", "bareFidalgoBay", "barePortGamble", "bareSkokomishRiver", "bareWillapaBay", "eelgrassCaseInlet", "eelgrassFidalgoBay", "eelgrassPortGamble", "eelgrassSkokomishRiver", "eelgrassWillapaBay") #Rename columns
head(averageProteinAreasMerged) #Confirm changes

#### ISOLATE PROTEINS FROM FINAL TRANSITION LIST ####

transitionList <- read.csv("2017-07-10-SRM-Transitions-With-PRTC.csv", header = FALSE) #Import final transition list
names(transitionList) <- c("V1", "V2", "retentionTime", "peptide", "protein", "V6", "V7") #Rename columns
transitionListTargetsOnly <- transitionList[-(124:153),] #Remove PRTC peptides
tail(transitionListTargetsOnly) #Confirm removal
transitionProteins <- transitionListTargetsOnly[,c(4,5)] #Since I want to merge by protein names to examine differential expression across sample sites, I don't need the transitions or peptides. I just need protein names and peak area. Therefore, I removed everything from the final transition list except for the column with protein names and peptides
head(transitionProteins) #Confirm changes
uniqueTransitionProteins <- unique(transitionProteins) #Removing duplicate entries from transitionProteins
head(uniqueTransitionProteins) #Confirm changes...didn't work. Will try again post-merge

#### MERGE FULL PROTEIN DATA WITH FINAL TRANSITION PROTEINS ####

uniqueTransitionProteinAreas <- merge(uniqueTransitionProteins, averageProteinAreasMerged, by = "protein")
head(uniqueTransitionProteinAreas) #Confirm merge
uniqueTransitionProteinAreas <- uniqueTransitionProteinAreas[,-2] #Remove peptide column
head(uniqueTransitionProteinAreas) #Confirm changes
uniqueTransitionProteinAreas <- unique(uniqueTransitionProteinAreas) #Remove duplicate rows.
head(uniqueTransitionProteinAreas) #Confirm changes. It works!

#### NMDS PLOT ####

uniqueTransitionProteinAreas[is.na(uniqueTransitionProteinAreas)] <- 0 #Replace NAs with 0s
head(uniqueTransitionProteinAreas) #Confirm changes

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- uniqueTransitionProteinAreas[-1]
rownames(area.protID2) <- uniqueTransitionProteinAreas[,1]

#Transpose the file so that rows and columns are switched and normalized by log(x+1)
area2.t <- t(area.protID2[,1:10])
area2.tra <- (area2.t+1)
area2.tra <- data.trans(area2.tra, method = 'log', plot = FALSE)

#Make MDS dissimilarity matrix
proc.nmds <- metaMDS(area2.tra, distance = 'bray', k = 2, trymax = 100, autotransform = FALSE)

#Make figure
fig.nmds <- ordiplot(proc.nmds, choices=c(1,2), type='none', display='sites', xlab='Axis 1', ylab='Axis 2', cex=0.5)
#bare=circle
#eelgrass=triangle
#case=red
#fidalgo=blue
#willapa=black
#skokomish=green
#gamble=magenta

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5))) #Looks like all of the points are on top of eachother.......
legend(0,0.06, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black')) 

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area2.tra dataset.

#Create heatmap
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .jpeg
jpeg(filename = "transitionHeatmap.jpeg")
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)
dev.off()

#### MERGE WITH GOTERMS AND PVALUES ####
sitesAndEelgrass <- read.csv(file = "/Users/yaaminivenkataraman/Documents/project-oyster-oa/analyses/DNR_TransitionSelection_20170707/2017-07-07-Preliminary-Transitions/2017-07-05-SitesEelgrass-Accession-nohead.csv")
head(sitesAndEelgrass)
goTermspValues <- sitesAndEelgrass[,c(1, 8, 18)] #Isolate GOterms and pvalues
head(goTermspValues) #Confirm isolation
names(goTermspValues) <- c("protein", "adj.pvalue", "goterm") #Rename columns
head(goTermspValues) #Confirm changes
uniqueTransitionProteinAreasGO <- merge(x = uniqueTransitionProteinAreas, y = goTermspValues, by = "protein")
head(uniqueTransitionProteinAreasGO) #This list gives me the p-values for each replicate. I can sort through this manually before I use the information in REVIGO.
write.csv(x = uniqueTransitionProteinAreasGO, file = "2017-07-18-Transition-Proteins-for-REVIGO.csv")