#### INSTALL MSSTATS ####

source("http://bioconductor.org/biocLite.R")
biocLite("MSstats")
library(MSstats)

#### IMPORT AND PROCESS DATA ####
peakAreas <- read.csv("2017-05-13-peptide-results-MSstats.csv", na.strings = "#N/A") #Import peak area data
names(peakAreas) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", "Run", "BioReplicate", "Intensity") #Rename columns, primarily to change "FileName" and "Area" to "Run" and "Intensity"
head(peakAreas) #Confirm changes
peakAreasUnique <- peakAreas[!duplicated(peakAreas), ] #Remove duplicated rows
QuantData <- dataProcess(peakAreasUnique) #Process data using default settings of log2 transformation and constant normalization based on reference signals (equalize means)
