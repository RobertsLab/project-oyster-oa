#In this script, I'll get preliminary methylation information from diploid and triploid oyster samples. Additionally, identify differentially methylated loci (DML) between treatments and ploidy status using pairwise comparisons.

# Load packages

setwd("/gscratch/scrubbed/yaaminiv/Hawes/analyses/methylKit/") #Set root directory

require(methylKit) #Load methylKit
require(dplyr) #Load dplyr

# Obtain session information

sessionInfo()

# Process methylation data

analysisFiles <- list("zr3644_1_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_2_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_3_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_4_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_5_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_6_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_7_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_8_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_9_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_10_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_11_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_12_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_13_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_14_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_15_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_16_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_17_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_18_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_19_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_20_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_21_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_22_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_23_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov",
                      "zr3644_24_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov") #Put all .cov files into a list for analysis.

sampleMetadata <- data.frame("sampleID" = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                                            "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                                            "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                                            "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6"),
                             "ploidyTreatment" = c(rep(0, times = 12),
                                                    rep(1, times = 12)),
                             "pHTreatment" = c(rep(0, times = 6),
                                               rep(1, times = 6),
                                               rep(0, times = 6),
                                               rep(1, times = 6))) #Create dataframe with metadata
head(sampleMetadata) #Confirm dataframe creation

#I'll use `methRead` to create a methylation object from the coverage files, and include sample ID and treatment information.

processedFiles <- methylKit::methRead(analysisFiles,
                                      sample.id = list("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                                                       "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                                                       "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                                                       "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6"),
                                      assembly = "roslin",
                                      treatment = sampleMetadata$ploidyTreatment,
                                      pipeline = "bismarkCoverage",
                                      mincov = 2) #Process files. Treatment specified based on ploidy status. Use mincov = 2 to quickly process reads.

processedFilteredFilesCov5 <- methylKit::filterByCoverage(processedFiles,
                                                          lo.count = 5, lo.perc = NULL,
                                                          high.count = NULL, high.perc = 99.9) %>%
  methylKit::normalizeCoverage(.) #Filter coverage information for minimum 5x coverage, and remove PCR duplicates by excluding data in the 99.9th percentile of coverage with hi.perc = 99.9. Normalize coverage between samples to avoid over-sampling reads from one sample during statistical testing

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

# Characterize general methylation

## Sample-specific descriptive statistics

nFiles <- 24 #Count number of samples
fileName <- data.frame("nameBase" = rep("general-stats/percent-CpG-methylation", times = nFiles),
                       "nameBase2" = rep("general-stats/percent-CpG-coverage", times = nFiles),
                       "sample.ID" = 1:24) #Create new dataframe for filenames
head(fileName) #Confirm dataframe creation

fileName$actualFileName1 <- paste(fileName$nameBase, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG methylation plot
fileName$actualFileName2 <- paste(fileName$nameBase2, "-Filtered", "-5xCoverage", "-Sample", fileName$sample.ID, ".jpeg", sep = "") #Create a new column for the full filename for filtered + 5x coverage + specific sample's percent CpG coverage plot
head(fileName) #Confirm column creation

### Create plots

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName1[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getMethylationStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
  dev.off() #Turn off plotting device
} #Plot and save %CpG methylation information

for(i in 1:nFiles) { #For each data file
  jpeg(filename = fileName$actualFileName2[i], height = 1000, width = 1000) #Save file with designated name
  methylKit::getCoverageStats(processedFilteredFilesCov5[[i]], plot = TRUE, both.strands = FALSE) #Get CpG coverage information
  dev.off() #Turn off plotting device
} #Plot and save CpG coverage information

## Comparative analysis

# methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
#                                                        destrand = FALSE) #Combine all processed files into a single table. Use destrand = TRUE to not destrand. By default only bases with data in all samples will be kept

# methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
#                                                        destrand = FALSE,
#                                                        min.per.group = 9L) #Combine all processed files into a single table. Use destrand = TRUE to not destrand. Based with data in at least 9/12 samples peer treatment will be included

methylationInformationFilteredCov5 <- methylKit::unite(processedFilteredFilesCov5,
                                                       destrand = FALSE,
                                                       min.per.group = 8L) #Combine all processed files into a single table. Use destrand = TRUE to not destrand. Based with data in at least 9/12 samples peer treatment will be included

head(methylationInformationFilteredCov5) #Confirm unite
length(methylationInformationFilteredCov5$chr) #1984530 CpG loci with data in all samples. 4557452 CpG loci with data when using min.per.group = 9, 5103729 with min.per.group = 8.

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

clusteringInformationFilteredCov5 <- methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = FALSE) #Save cluster information as a new object

jpeg(filename = "general-stats/Full-Sample-Pearson-Correlation-Plot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::getCorrelation(methylationInformationFilteredCov5, plot = TRUE) #Understand correlation between methylation patterns in different samples
dev.off()

jpeg(filename = "general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::clusterSamples(methylationInformationFilteredCov5, dist = "correlation", method = "ward", plot = TRUE) #Cluster samples based on correlation coefficients
dev.off()

jpeg(filename = "general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5) #Run a PCA analysis on percent methylation for all samples
dev.off() #Turn off plotting device

jpeg(filename = "general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand.jpeg", height = 1000, width = 1000) #Save file with designated name
methylKit::PCASamples(methylationInformationFilteredCov5, screeplot = TRUE) #Run the PCA analysis and plot variances against PC number in a screeplot
dev.off()

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

# Differentially methylated loci

#I will examine ploidy and treatment differences separately, using the other variable as a covariate.

## Ploidy differences

### Create covariate matrix

covariatepH <- data.frame("pH" = c(rep("H", times = 6),
                                   rep("L", times = 6),
                                   rep("H", times = 6),
                                   rep("L", times = 6))) #Create dataframe with pH covariate information
head(covariatepH) #Check dataframe format

### Identify DML

differentialMethylationStatsPloidy <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5,
                                                                   covariates = covariatepH,
                                                                   overdispersion = "MN", test = "Chisq") #Calculate differential methylation statistics based on treatment indication from methRead. Include pH as a covariate.
head(differentialMethylationStatsPloidy) #Look at differential methylation statistics

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

diffMethStatsPloidy25 <- methylKit::getMethylDiff(differentialMethylationStatsPloidy, difference = 25, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsPloidy25$chr) #Count the number of DML: 29
head(diffMethStatsPloidy25) #Confirm creation

diffMethStatsPloidy50 <- methylKit::getMethylDiff(differentialMethylationStatsPloidy, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsPloidy50$chr) #Count the number of DML: 1
head(diffMethStatsPloidy50) #Confirm creation

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

write.csv(diffMethStatsPloidy25, "DML/DML-ploidy-25-Cov5.csv", quote = FALSE) #Save table as .csv
write.csv(diffMethStatsPloidy50, "DML/DML-ploidy-50-Cov5.csv", quote = FALSE) #Save table as .csv

## pH differences

### Assign treatment information

methylationInformationFilteredCov5T <- methylKit::reorganize(methylationInformationFilteredCov5,
                                                             sample.id = sampleMetadata$sampleID,
                                                             treatment = sampleMetadata$pHTreatment) #Reorganize methylationInformationFilteredCov5 to provide pH treatment specification.
head(methylationInformationFilteredCov5T) #Confirm unite
length(methylationInformationFilteredCov5T$chr) #1984530 CpGs with data in all samples. 4557452 CpGs with data in 9/12 samples using min.per.group = 9, 5103729 with min.per.group = 8.

### Create covariate matrix

covariatePloidy <- data.frame("ploidy" = c(rep("2N", times = 12),
                                           rep("3N", times = 12))) #Create dataframe with ploidy covariate information
head(covariatePloidy) #Check dataframe format

### Identify DML

differentialMethylationStatsTreatment <- methylKit::calculateDiffMeth(methylationInformationFilteredCov5T,
                                                                      covariates = covariatePloidy,
                                                                      overdispersion = "MN", test = "Chisq") #Calculate differential methylation statistics based on treatment indication from methRead. Include pH as a covariate.
head(differentialMethylationStatsTreatment) #Look at differential methylation statistics

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

diffMethStatsTreatment25 <- methylKit::getMethylDiff(differentialMethylationStatsTreatment, difference = 25, qvalue = 0.01) #Identify loci that are at least 75% different
length(diffMethStatsTreatment25$chr) #Count the number of DML: 42
head(diffMethStatsTreatment25) #Confirm creation

diffMethStatsTreatment50 <- methylKit::getMethylDiff(differentialMethylationStatsTreatment, difference = 50, qvalue = 0.01) #Identify loci that are at least 50% different
length(diffMethStatsTreatment50$chr) #Count the number of DML: 3
head(diffMethStatsTreatment50) #Confirm creation

save.image("methylKit.RData") #Save R Data in case R crashes
#load("methylKit.RData") #Load R Data

write.csv(diffMethStatsTreatment25, "DML/DML-pH-25-Cov5.csv", quote = FALSE) #Save table as .csv
write.csv(diffMethStatsTreatment50, "DML/DML-pH-50-Cov5.csv", quote = FALSE) #Save table as .csv
