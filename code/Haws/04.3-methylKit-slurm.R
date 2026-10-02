#!/usr/bin/env Rscript
# Slurm version of 04.2-methylKit-parameter-testing.Rmd. Same methylKit calls, split into stages so the
# 12 DML tests (ploidy and pH x min.per.group All/11/10/9/8 and outliers removed) run as a job array.
#
# Run from analyses/Haws_04.3-methylKit-slurm/ (04.3-methylKit-slurm-submit.sh handles this):
#   Rscript 04.3-methylKit-slurm.R prep
#   Rscript 04.3-methylKit-slurm.R dml <none|MN|shrinkMN>      #Task 1-12 taken from SLURM_ARRAY_TASK_ID
#   Rscript 04.3-methylKit-slurm.R summary <none|MN|shrinkMN>
#
# Environment variables:
#   HI_PERC   Upper coverage percentile filter (default 99.9). Use "none" to reproduce 04.2 as written: 04.2 passes
#             high.perc, which filterByCoverage does not recognize, so no upper filter is applied there.

suppressPackageStartupMessages({
  library(tidyverse)
  library(methylKit)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
overdispersion <- ifelse(length(args) >= 2, args[2], "none")
stopifnot(mode %in% c("prep", "dml", "summary"),
          overdispersion %in% c("none", "MN", "shrinkMN"))

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) #Use all cores allocated by Slurm
hiPercSetting <- Sys.getenv("HI_PERC", "99.9")
hiPerc <- if (hiPercSetting == "none") NULL else as.numeric(hiPercSetting)

dataDir <- "../../data/Haws"
prepDir <- paste0("prep-hiperc-", hiPercSetting) #Separate prep output for each upper coverage filter
dmlDir <- paste0("DML-hiperc-", hiPercSetting, "-overdispersion-", overdispersion)

outliers <- c(3, 14) #Samples 2H-3 and 3H-2

tasks <- expand.grid(setting = c("All", "11", "10", "9", "8", "OutRM"),
                     test = c("ploidy", "pH"),
                     stringsAsFactors = FALSE)[, c("test", "setting")] #Tasks 1-6 are ploidy, 7-12 are pH

# Sample metadata, as in 04.2

sampleMetadata <- read.csv(file.path(dataDir, "sample_metadata.csv")) %>%
  dplyr::select(-c(1)) %>%
  arrange(sample_number) %>%
  mutate("sampleID" = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                        "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                        "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                        "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6")) %>%
  dplyr::rename(pH = ph) %>%
  mutate(ploidyTreatment = case_when(ploidy == "3N" ~ 1,
                                     ploidy == "2N" ~ 0)) %>%
  mutate(pHTreatment = case_when(pH == "high" ~ 0,
                                 pH == "low" ~ 1))
stopifnot(identical(sampleMetadata$sample_number, 1:24)) #File zr3644_N must match row N

filterNormalize <- function(methylObj) {
  methylKit::filterByCoverage(methylObj,
                              lo.count = 5, lo.perc = NULL,
                              hi.count = NULL, hi.perc = hiPerc) %>%
    methylKit::normalizeCoverage(.) #Minimum 5x coverage, optional upper percentile filter, then median normalization
}

# Stage 1: read, filter, and normalize coverage files once for all tasks

if (mode == "prep") {
  dir.create(prepDir, showWarnings = FALSE)
  analysisFiles <- as.list(file.path(dataDir, paste0("zr3644_", 1:24, "_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov")))
  stopifnot(all(file.exists(unlist(analysisFiles))))

  processedFiles <- methylKit::methRead(analysisFiles,
                                        sample.id = as.list(sampleMetadata$sampleID),
                                        assembly = "oyster_v9",
                                        treatment = sampleMetadata$ploidyTreatment,
                                        pipeline = "bismarkCoverage",
                                        mincov = 2)

  processedFilesOutRM <- methylKit::reorganize(processedFiles,
                                               sample.ids = sampleMetadata$sampleID[-outliers],
                                               treatment = sampleMetadata$ploidyTreatment[-outliers]) #Drop outliers before filtering, since normalization depends on which samples are included

  saveRDS(filterNormalize(processedFiles), file.path(prepDir, "filtered-All.rds"), compress = FALSE)
  saveRDS(filterNormalize(processedFilesOutRM), file.path(prepDir, "filtered-OutRM.rds"), compress = FALSE)
  file.create(file.path(prepDir, "prep-complete")) #Lets the submit script skip prep on later runs
}

# Stage 2: one DML test per array task

if (mode == "dml") {
  taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  stopifnot(taskID %in% seq_len(nrow(tasks)))
  test <- tasks$test[taskID]
  setting <- tasks$setting[taskID]
  message("Task ", taskID, ": ", test, ", min.per.group = ", setting, ", overdispersion = ", overdispersion, ", hi.perc = ", hiPercSetting, ", cores = ", nCores)

  for (d in file.path(dmlDir, c("counts", "rds"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dir.create("general-stats", showWarnings = FALSE)

  keep <- if (setting == "OutRM") setdiff(1:24, outliers) else 1:24
  metadata <- sampleMetadata[keep, ]
  filteredFiles <- readRDS(file.path(prepDir, ifelse(setting == "OutRM", "filtered-OutRM.rds", "filtered-All.rds")))

  treatment <- if (test == "ploidy") metadata$ploidyTreatment else metadata$pHTreatment
  covariates <- if (test == "ploidy") data.frame(pH = metadata$pH) else data.frame(ploidy = metadata$ploidy) #Other factor is the covariate
  filteredFiles <- methylKit::reorganize(filteredFiles,
                                         sample.ids = metadata$sampleID,
                                         treatment = treatment) #Set treatment for this test

  if (setting %in% c("All", "OutRM")) {
    methylationInformation <- methylKit::unite(filteredFiles, destrand = FALSE, mc.cores = 2) #Keep bases with data in all samples
  } else {
    methylationInformation <- methylKit::unite(filteredFiles, destrand = FALSE, mc.cores = 2,
                                               min.per.group = as.integer(setting))
  }

  if (test == "ploidy" & setting == "OutRM") { #Clustering and PCA plots from the outlier removal section of 04.2
    jpeg(filename = "general-stats/Full-Sample-Pearson-Correlation-Plot-FilteredCov5Destrand-OutRM.jpeg", height = 1000, width = 1000)
    methylKit::getCorrelation(methylationInformation, plot = TRUE)
    dev.off()
    jpeg(filename = "general-stats/Full-Sample-CpG-Methylation-Clustering-FilteredCov5Destrand-OutRM.jpeg", height = 1000, width = 1000)
    methylKit::clusterSamples(methylationInformation, dist = "correlation", method = "ward", plot = TRUE)
    dev.off()
    jpeg(filename = "general-stats/Full-Sample-Methylation-PCA-FilteredCov5Destrand-OutRM.jpeg", height = 1000, width = 1000)
    methylKit::PCASamples(methylationInformation)
    dev.off()
    jpeg(filename = "general-stats/Full-Sample-Methylation-Screeplot-FilteredCov5Destrand-OutRM.jpeg", height = 1000, width = 1000)
    methylKit::PCASamples(methylationInformation, screeplot = TRUE)
    dev.off()
  }

  differentialMethylationStats <- methylKit::calculateDiffMeth(methylationInformation,
                                                               covariates = covariates,
                                                               overdispersion = overdispersion,
                                                               mc.cores = nCores) #With overdispersion = "none" this uses a Chisq test; with "MN" or "shrinkMN" it uses the F-test
  saveRDS(differentialMethylationStats, file.path(dmlDir, "rds", paste0("diffMeth-", test, "-", setting, ".rds")), compress = FALSE) #Keep all CpG results for other thresholds

  settingLabel <- ifelse(setting == "OutRM", "All-OutRM", setting) #Match 04.2 file names
  nDML <- c()
  for (difference in c(25, 50)) {
    diffMethStats <- methylKit::getMethylDiff(differentialMethylationStats, difference = difference, qvalue = 0.01) #Loci at least 25% or 50% different
    write_delim(diffMethStats, file.path(dmlDir, paste0("DML-", test, "-", difference, "-Covar-Cov5-", settingLabel, ".bed")), delim = "\t", quote = "none")
    nDML[as.character(difference)] <- nrow(diffMethStats)
  }

  write_tsv(tibble(test = test,
                   samples = ifelse(setting == "OutRM", "Outliers removed", "All"),
                   min.per.group = ifelse(setting == "OutRM", "All", setting),
                   CpGs.tested = nrow(methylationInformation),
                   DML.25 = nDML[["25"]],
                   DML.50 = nDML[["50"]]),
            file.path(dmlDir, "counts", paste0(test, "-", setting, ".tsv")))
}

# Stage 3: combine per-task counts into the 04.2 summary table

if (mode == "summary") {
  countFiles <- list.files(file.path(dmlDir, "counts"), pattern = "\\.tsv$", full.names = TRUE)
  if (length(countFiles) < nrow(tasks)) warning("Only ", length(countFiles), " of ", nrow(tasks), " tasks have results")

  counts <- map_dfr(countFiles, read_tsv, col_types = cols(.default = "c")) %>%
    mutate(across(c(CpGs.tested, DML.25, DML.50), as.integer),
           samples = factor(samples, levels = c("All", "Outliers removed")),
           min.per.group = factor(min.per.group, levels = c("All", "11", "10", "9", "8"))) %>%
    arrange(test, samples, min.per.group)
  write_csv(counts, file.path(dmlDir, "DML-counts-long.csv")) #Includes number of CpGs tested

  summaryDML <- counts %>%
    dplyr::select(-CpGs.tested) %>%
    pivot_wider(names_from = test, values_from = c(DML.25, DML.50), names_glue = "{test}-DML-{sub('DML.', '', .value)}") %>%
    arrange(samples, min.per.group) %>%
    dplyr::select(samples, min.per.group, any_of(c("ploidy-DML-25", "ploidy-DML-50", "pH-DML-25", "pH-DML-50")))
  print(as.data.frame(summaryDML))
  write.csv(summaryDML, file.path(dmlDir, "DML-summary-table.csv"), row.names = FALSE, quote = FALSE) #Same layout as 04.2
}

sessionInfo()
