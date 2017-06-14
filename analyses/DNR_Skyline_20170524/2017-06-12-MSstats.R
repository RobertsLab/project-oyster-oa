#### INSTALL MSSTATS ####
source("http://bioconductor.org/biocLite.R")
biocLite("MSstats")
library(MSstats)

#### IMPORT AND PROCESS DATA ####
peakAreas <- read.csv("2017-06-10-peptide-transition-results-MSstats-no-pivot-error-checked.csv", na.strings = "#N/A") #Import peak area data
names(peakAreas) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", "Run", "BioReplicate", "Intensity") #Rename columns, primarily to change "FileName" and "Area" to "Run" and "Intensity"
head(peakAreas) #Confirm changes
QuantData <- dataProcess(peakAreas) #Process data using default settings of log2 transformation and constant normalization based on reference signals (equalize means). Message: ** There are 1295 intensities which are zero or less than 1. These intensities are replaced with 1. Error in dataProcess(peakAreas) : ** MSstats suspects that there are fractionations and potentially technical replicates too. Please add Fraction column in the input.
peakAreasUnique <- peakAreas[!duplicated(peakAreas), ] #Remove duplicated rows
QuantData <- dataProcess(peakAreasUnique) #** There are 1295 intensities which are zero or less than 1. These intensities are replaced with 1. Error in dataProcess(peakAreasUnique) : ** MSstats suspects that there are fractionations and potentially technical replicates too. Please add Fraction column in the input.
head(QuantData$ProcessedData)
