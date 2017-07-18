#### IMPORT DATA ####

proteinAreas <- read.csv("2017-06-10-protein-areas-only-error-checked.csv", na.strings = "#N/A") #Specify Skyline's special way of designating N/A values
head(proteinAreas) #Confirm successful import
averageProteinAreas <- aggregate(proteinAreas[-1], proteinAreas[1], mean, na.action = na.omit, na.rm = TRUE) #Average protein areas across transitions and peptides
head(averageProteinAreas) #Confirm changes 

#### REFORMAT DATA ####

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

averageProteinAreasMerged <- data.frame(averageProteinAreas$protein, bareCaseInlet, bareFidalgoBay, barePortGamble, bareSkokomishRiver, bareWillapaBay, eelgrassCaseInlet, eelgrassFidalgoBay, eelgrassPortGamble, eelgrassSkokomishRiver, eelgrassWillapaBay)

#### FULL NMDS PLOT ####

averageProteinAreasMerged[is.na(averageProteinAreasMerged)] <- 0 #Replace NAs with 0s
head(averageProteinAreasMerged) #Confirm changes

#Load the source file for the biostats package
source("biostats.R") #Either load the source R script or copy paste. Must run this code before NMDS.
install.packages("vegan") #Install vegan package
library(vegan)

#Make sure first column of protein names is recognized as row names instead of values
area.protID2 <- averageProteinAreasMerged[-1]
rownames(area.protID2) <- averageProteinAreasMerged[,1]

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

points(fig.nmds, 'sites', col=c('red', 'blue', 'black', 'green', 'magenta','red', 'blue', 'black', 'green', 'magenta'), pch=c(rep(16,5), rep(17,5)))
legend(-0.045,0.025, pch=c(rep(16,5), 1, 2), legend=c('Case Inlet', "Fidalgo Bay", "Willapa Bay", "Skokomish", "Port Gamble", "Bare", "Eelgrass"), col=c('red', 'blue', 'black', 'green', 'magenta', 'black', 'black'))

#### FULL HEATMAP ####

#Install package
install.packages("pheatmap")
library(pheatmap)

#Data should be log(x) or log(x+1) transformed for this analysis, so I'll use my area2.tra dataset.

#Create heatmap
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)

#Export preliminary heatmap as a .png
png(filename = "fullHeatmap.png")
pheatmap(area2.tra, cluster_rows = T, cluster_cols = T, clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'average', show_rownames = T, show_colnames = F)
dev.off()

#### MERGE WITH GO TERMS ####

#Upload table with proteins and accession codes.
proteinAccessionCodes <- read.table(file = "background-proteome-accession.txt", header = FALSE, col.names = c("averageAreaAdjusted.proteins", "accession", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"))
proteinAccessionCodes <- within(proteinAccessionCodes, rm("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")) #Removing the extra columns
names(proteinAccessionCodes) <- c("averageProteinAreas.protein", "accession")
head(proteinAccessionCodes) #View uploaded data and confirm changes
head(averageProteinAreasMerged)

#Merge Skyline output with Uniprot information
skylineProteinAccession <- merge(x = averageProteinAreasMerged, y = proteinAccessionCodes, by = "averageProteinAreas.protein")
head(skylineProteinAccession) #confirm merge

#### WRITE OUT MERGED TABLE ####

#Write out alltreatments_DEG_Uniprot as a tab file. Remove row and column names using "row.names" and "col.names" arguments
write.table(skylineProteinAccession, "2017-06-13-Skyline-ProteinAccession-nohead.txt", col.names = F, row.names = F)

#### MERGE WITH ANNOTATIONS ####
proteinAnnotationsEvalues <- read.csv(file = "/Users/yaaminivenkataraman/Documents/project-oyster-oa/analyses/DNR_TransitionSelection_20170707/2017-07-07-Preliminary-Transitions/2017-07-07-Gigas-Annotations-Evalues.csv", header = FALSE, col.names = c("C1", "averageProteinAreas.protein", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "Evalue", "C14", "C15", "reviewed", "Annotation", "C18", "Species", "C20", "BiologicalProcess", "GOTerm", "Pathway", "C24", "C25", "C26")) #Import full annotation
fullSkylineAnnotations <- merge(x = averageProteinAreasMerged, y = proteinAnnotationsEvalues, by = "averageProteinAreas.protein")
write.csv(fullSkylineAnnotations, "2017-07-09-Full-Skyline-Output-Annotations.csv")
