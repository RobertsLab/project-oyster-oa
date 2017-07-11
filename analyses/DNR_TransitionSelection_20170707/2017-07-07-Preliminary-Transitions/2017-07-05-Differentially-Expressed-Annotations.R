##### DIFFERENTIAL PROTEIN ANNOTATIONS #####

#### IMPORT ACCESSION CODES ####

proteinAccessionCodes <- read.table(file = "background-proteome-accession.txt", header = FALSE, col.names = c("averageAreaAdjusted.proteins", "accession", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")) #Import accession codes
proteinAccessionCodes <- within(proteinAccessionCodes, rm("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")) #Removing the extra columns
names(proteinAccessionCodes) <- c("Protein", "accession")
head(proteinAccessionCodes) #View uploaded data and confirm changes

#### IMPORT ANNOTATIONS ####

proteinAnnotations <- read.csv(file = "2017-07-05-Gigas-Annotations.csv", header = FALSE, col.names = c("Protein", "Code1", "Annotation", "Code2", "Species", "BiologicalProcess", "GOTerms", "Pathway", "ProteinCode")) #Import annotations
head(proteinAnnotations)

#### MERGE ANNOTATIONS AND ACCESSION CODES ####

proteinAnnotationsAccessions <- merge(x = proteinAccessionCodes, y = proteinAnnotations, by = "Protein")

#### BARE VS. EELGRASS ####

bareEelgrass <- read.csv("2017-06-23-MSstats-BarevEelgrass-Differential-Expression.csv", na.strings = "NA") #Import data
head(bareEelgrass) #Confirm import. protein = "Protein"
bareEelgrassProteinAccession <- merge(x = bareEelgrass, y = proteinAnnotationsAccessions, by = "Protein") #Merge Skyline output with Uniprot and accession information
head(bareEelgrassProteinAccession) #Confirm merge
write.csv(bareEelgrassProteinAccession, "2017-07-05-BarevEelgrass-Accession-nohead.csv", col.names = F, row.names = F) #Write out csv file

bareEelgrassAnnotations <- merge(x = bareEelgrass, y = proteinAnnotations, by = "Protein") #Merge Skyline with annotations only
write.csv(bareEelgrassAnnotations, "2017-07-05-BarevEelgrass-Annotations-nohead.csv", col.names = F, row.names = F) #Write out csv file

#### SITES ONLY ####

sitesOnly <- read.csv("2017-06-30-MSstats-Sites-Differential-Expression.csv", na.strings = "NA") #Import data
head(bareEelgrass) #Confirm import. protein = "Protein"
sitesOnlyProteinAccession <- merge(x = sitesOnly, y = proteinAnnotationsAccessions, by = "Protein") #Merge Skyline output with Uniprot and accession information
head(sitesOnlyProteinAccession) #Confirm merge
write.csv(sitesOnlyProteinAccession, "2017-07-05-SitesOnly-Accession-nohead.csv", col.names = F, row.names = F) #Write out csv file

sitesOnlyAnnotations <- merge(x = sitesOnly, y = proteinAnnotations, by = "Protein") #Merge Skyline with annotations only
write.csv(sitesOnlyAnnotations, "2017-07-05-SitesOnly-Annotations-nohead.csv", col.names = F, row.names = F) #Write out csv file

#### SITES AND EELGRASS ####

sitesEelgrass <- read.csv("2017-06-30-MSstats-Sites-Eelgrass-Differential-Expression.csv", na.strings = "NA") #Import data
head(sitesEelgrass) #Confirm import. protein = "Protein"
sitesEelgrassProteinAccession <- merge(x = sitesEelgrass, y = proteinAnnotationsAccessions, by = "Protein") #Merge Skyline output with Uniprot and accession information
head(sitesEelgrassProteinAccession) #Confirm merge
write.csv(sitesEelgrassProteinAccession, "2017-07-05-SitesEelgrass-Accession-nohead.csv", col.names = F, row.names = F) #Write out csv file

sitesEelgrassAnnotations <- merge(x = sitesEelgrass, y = proteinAnnotations, by = "Protein") #Merge Skyline with annotations only
write.csv(sitesEelgrassAnnotations, "2017-07-05-SitesEelgrass-Annotations-nohead.csv", col.names = F, row.names = F) #Write out csv file

#### IMPORT ANNOTATIONS WITH EVALUES ####

proteinAnnotationsEvalues <- read.csv(file = "2017-07-07-Gigas-Annotations-Evalues.csv", header = FALSE, col.names = c("C1", "Protein", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "Evalue", "C14", "C15", "reviewed", "Annotation", "C18", "Species", "C20", "BiologicalProcess", "GOTerm", "Pathway", "C24", "C25", "C26"))

#### IMPORT PRELIMINARY SHORTLIST ####

preliminaryShortlist <- read.csv("2017-07-07-Protein-Shortlist.csv", header = TRUE)

#### MERGE PRELIMINARY SHORTLIST WITH EVALUES ####

preliminaryShortlistEvalues <- merge(x = preliminaryShortlist, y = proteinAnnotationsEvalues, by = "Protein")
write.csv(preliminaryShortlistEvalues, "2017-07-07-Protein-Shortlist-Evalues.csv", col.names = F, row.names = F) #Write out csv file

#### IMPORT PRELIMINARY TRANSITION TARGETS ####

preliminaryTransitions <- read.csv("2017-07-07-Preliminary-Target-Transitions.csv")
preliminaryTransitions <- within(preliminaryTransitions, rm("X")) #Remove extra column
preliminaryTransitions <- preliminaryTransitions[-28,] #Remove extra row
head(preliminaryTransitions) #Confirm changes
tail(preliminaryTransitions) #Confirm changes

#### MERGE PRELIMINARY TRANSITION LIST WITH EVALUES ####

preliminaryTransitionsEvalue <- merge(x = preliminaryTransitions, y = proteinAnnotationsEvalues, by = "Protein")
write.csv(preliminaryTransitionsEvalue, "2017-07-07-Preliminary-Target-Transitions-Evalues.csv", col.names = F, row.names = F) #Write out csv file