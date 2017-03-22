#Step 1: Upload data
averageAreaAdjustedMerged <- read.csv("/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Oyster-AverageAdjustedMergedArea.csv", row.names = NULL, na.strings = "NA")
averageAreaAdjustedMerged <- within(averageAreaAdjustedMerged, rm("X")) #Remove the extra column "X"
averageAreaAdjustedMerged #View uploaded data and confirm changes

#Step 2: Upload table with proteins and accession codes.
proteinAccessionCodes <- read.table(file = "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/background-proteome-accession.txt", header = FALSE, col.names = c("averageAreaAdjusted.proteins", "accession", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"))
proteinAccessionCodes <- within(proteinAccessionCodes, rm("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")) #Removing the extra columns
proteinAccessionCodes #View uploaded data and confirm changes

#Step 3: Merge Skyline output with Uniprot information
skylineProteinAccession <- merge(x = averageAreaAdjustedMerged, y = proteinAccessionCodes, by = "averageAreaAdjusted.proteins")
head(skylineProteinAccession) #confirm merge

#Step 4: Write out alltreatments_DEG_Uniprot as a tab file. Remove row and column names using "row.names" and "col.names" arguments
write.table(skylineProteinAccession, "/Users/yaamini/Documents/project-oyster-oa/analyses/DNR_Preliminary_Analyses_20170321/Skyline-ProteinAccession-nohead.txt", col.names = F, row.names = F)