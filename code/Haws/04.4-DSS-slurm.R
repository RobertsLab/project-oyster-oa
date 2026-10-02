#!/usr/bin/env Rscript
# DSS differential methylation for the Haws ploidy x pH experiment. Plan: 04.4-DSS-plan.md.
# Stages follow 04.3-methylKit-slurm.R. 04.4-DSS-slurm-submit.sh runs them from analyses/Haws_04.4-DSS/:
#   Rscript 04.4-DSS-slurm.R prep-raw   #Read the 24 coverage files once into count matrices (all CpGs, >= 1x)
#   Rscript 04.4-DSS-slurm.R qc         #Sample QC table, PCA, and flags (needs prep-raw)
#   Rscript 04.4-DSS-slurm.R prep       #Apply coverage and presence filters for one setting (needs prep-raw)
#   Rscript 04.4-DSS-slurm.R fit        #DMLfit.multiFactor + DMLtest.multiFactor on one block of CpGs; block from SLURM_ARRAY_TASK_ID
#   Rscript 04.4-DSS-slurm.R summary    #Combine blocks, recompute FDR over all CpGs, call DML/DMR, count
#   Rscript 04.4-DSS-slurm.R compare    #Sample-set sensitivity: do the All DML hold up in Drop3H2, DropPC2, OutRM? (needs their summaries)
#   Rscript 04.4-DSS-slurm.R permsummary #Observed DML counts vs label permutations (needs PERM = 0 and PERM > 0 summaries)
#   Rscript 04.4-DSS-slurm.R snp-prep    #Mark CpGs with BS-SNPer SNPs (prep-bssnper/*.vcf from gannet; needs prep-raw)
#   Rscript 04.4-DSS-slurm.R snp-summary #SNP enrichment among observed and permuted DML; genotype PCA vs methylation PCA
#
# Environment variables for prep, fit, and summary (defaults in brackets):
#   LO_COV    [5]      Per-sample minimum coverage. Use 1 for no filter (DSS default)
#   HI_PERC   [99.9]   Per-sample upper coverage percentile, calculated on nuclear CpGs only. Mito CpGs are exempt,
#                      since mito coverage (~300x) is far above the nuclear cutoff (~80x). Use "none" to skip
#   PRESENCE  [cell5]  "all" = coverage in every sample; "cellK" = coverage in at least K samples of every ploidy x pH group
#   MIN_METH  [10]     Low-variation filter: keep CpGs whose mean % methylation over all samples (ignoring labels) is between
#                      MIN_METH and 100 - MIN_METH. Most oyster CpGs are ~0% or ~100% in every sample; the test has no power
#                      there, but they still count in the FDR correction. A difference between two halves of the samples is
#                      at most 2 x mean (or 2 x (100 - mean)), so MIN_METH = 10 removes only CpGs with differences of about 20%
#                      or less. The filter ignores labels, so it does not bias the tests and works the same for permutations.
#                      Mito CpGs are exempt (~2% methylation in every sample, so the filter would remove all of them).
#                      "none" skips it (runs before 2026-10-02)
#   SNP_FILTER [none]  "any" = drop CpGs where any sample has a BS-SNPer PASS SNP at the C or at the G on the other strand
#                      (either makes the site look unmethylated). Needs snp-prep
#   SAMPLES   [All]    Sample set. "All" (24 oysters, primary). Sensitivity checks chosen from the sample QC (2026-10-02):
#                      "Drop3H2" (PC1 outlier, lowest correlation with other samples), "DropPC2" (2H-1, 2H-2, 3H-5, which
#                      separate on PC2), and "OutRM" (2H-3 and 3H-2, the pair removed in 04.2/04.3, for comparison with methylKit)
#   MODEL     [interaction]  "interaction" = ~ ploidy * pH (primary); "additive" = ~ ploidy + pH (sensitivity)
#   PERM      [0]      0 = observed labels. Any other value is the seed for a label permutation (null check)
#   NCHUNK    [1]      Number of CpG blocks for the fit array. The fit is fast (~0.5 s per 7k CpGs x 24 samples), so 1 block
#                      (dispersion prior estimated from all CpGs) is the default

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
stopifnot(mode %in% c("prep-raw", "qc", "prep", "fit", "summary", "compare", "permsummary", "snp-prep", "snp-summary"))

loCov <- as.integer(Sys.getenv("LO_COV", "5"))
hiPercSetting <- Sys.getenv("HI_PERC", "99.9")
presence <- Sys.getenv("PRESENCE", "cell5")
minMethSetting <- Sys.getenv("MIN_METH", "10")
snpFilter <- Sys.getenv("SNP_FILTER", "none")
samples <- Sys.getenv("SAMPLES", "All")
model <- Sys.getenv("MODEL", "interaction")
perm <- as.integer(Sys.getenv("PERM", "0"))
nChunk <- as.integer(Sys.getenv("NCHUNK", "1"))
stopifnot(loCov >= 1,
          hiPercSetting == "none" || !is.na(as.numeric(hiPercSetting)),
          presence == "all" || grepl("^cell[0-9]+$", presence),
          snpFilter %in% c("none", "any"),
          minMethSetting == "none" || (!is.na(as.numeric(minMethSetting)) && as.numeric(minMethSetting) < 50),
          samples %in% c("All", "Drop3H2", "DropPC2", "OutRM"),
          model %in% c("interaction", "additive"))

dataDir <- "../../data/Haws"
mitoChr <- "NC_001276.1"
dropSets <- list(All = integer(0),
                 Drop3H2 = 14, #3H-2
                 DropPC2 = c(1, 2, 17), #2H-1, 2H-2, 3H-5
                 OutRM = c(3, 14)) #2H-3 and 3H-2, as in 04.2/04.3
outliers <- dropSets$OutRM

filterTag <- paste0("cov", loCov, "-hiperc", hiPercSetting, "-", presence, ifelse(minMethSetting == "none", "", paste0("-meth", minMethSetting)),
                    ifelse(snpFilter == "none", "", "-noSNP"))
settingTag <- paste0(filterTag, "-", samples)
prepDir <- paste0("prep-", settingTag)
runDir <- paste0("DSS-", settingTag, "-", model, ifelse(perm == 0, "", paste0("-perm", perm)))

# Sample metadata, as in 04.2

sampleMetadata <- read.csv(file.path(dataDir, "sample_metadata.csv")) %>%
  dplyr::select(-c(1)) %>%
  arrange(sample_number) %>%
  mutate("sampleID" = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                        "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                        "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                        "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6")) %>%
  dplyr::rename(pH = ph) %>%
  mutate(cell = paste(ploidy, pH, sep = "-"))
stopifnot(identical(sampleMetadata$sample_number, 1:24)) #File zr3644_N must match row N

keep <- setdiff(1:24, dropSets[[samples]])

# Design for one run. Effect coding (+1/2 vs -1/2) makes the ploidy and pH coefficients average over the other
# factor even with the interaction in the model, and for unequal group sizes (sample sets with samples dropped). A positive coefficient means
# higher methylation in triploids (ploidyEffect) or at low pH (pHEffect)
makeDesign <- function(metadata, perm) {
  if (perm != 0) {
    set.seed(perm)
    metadata <- metadata %>% #Shuffle ploidy within pH, then pH within the new ploidy labels. Keeps every group the same size
      group_by(pH) %>% mutate(ploidy = sample(ploidy)) %>% ungroup() %>%
      group_by(ploidy) %>% mutate(pH = sample(pH)) %>% ungroup() %>%
      mutate(cell = paste(ploidy, pH, sep = "-"))
  }
  data.frame(ploidyEffect = ifelse(metadata$ploidy == "3N", 0.5, -0.5),
             pHEffect = ifelse(metadata$pH == "low", 0.5, -0.5),
             cell = metadata$cell,
             row.names = metadata$sampleID)
}

# Mean difference in % methylation for each effect. Average of the group means, so each ploidy x pH group counts
# equally, as in the model. Samples with no coverage at a CpG are ignored
effectSizes <- function(M, Cov, design) {
  beta <- M / Cov #NaN where Cov = 0
  beta[Cov == 0] <- NA
  groupMean <- sapply(c("2N-high", "2N-low", "3N-high", "3N-low"),
                      function(g) rowMeans(beta[, design$cell == g, drop = FALSE], na.rm = TRUE))
  data.frame(meanMeth = 100 * rowMeans(groupMean),
             diffPloidy = 100 * ((groupMean[, "3N-high"] + groupMean[, "3N-low"]) - (groupMean[, "2N-high"] + groupMean[, "2N-low"])) / 2,
             diffpH = 100 * ((groupMean[, "2N-low"] + groupMean[, "3N-low"]) - (groupMean[, "2N-high"] + groupMean[, "3N-high"])) / 2,
             diffInteraction = 100 * ((groupMean[, "3N-low"] - groupMean[, "2N-low"]) - (groupMean[, "3N-high"] - groupMean[, "2N-high"])))
}

# Stage 1: read coverage files into count matrices (union of CpGs covered in any sample)

if (mode == "prep-raw") {
  dir.create("prep-raw", showWarnings = FALSE)
  covFiles <- file.path(dataDir, paste0("zr3644_", 1:24, "_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov"))
  stopifnot(all(file.exists(covFiles)))

  covData <- lapply(covFiles, function(f) {
    message("Reading ", basename(f))
    fread(f, select = c(1, 2, 5, 6), col.names = c("chr", "pos", "X", "U")) #Bismark coverage: chr, start, end, %meth, count methylated, count unmethylated
  })
  chrLevels <- sort(unique(unlist(lapply(covData, function(d) unique(d$chr)))))
  posKey <- function(d) match(d$chr, chrLevels) * 1e9 + d$pos #One number per CpG (positions are < 1e9)
  keys <- sort(unique(unlist(lapply(covData, posKey))))

  M <- matrix(0L, nrow = length(keys), ncol = 24, dimnames = list(NULL, sampleMetadata$sampleID))
  Cov <- M
  for (i in 1:24) {
    rowIndex <- match(posKey(covData[[i]]), keys)
    M[rowIndex, i] <- as.integer(covData[[i]]$X)
    Cov[rowIndex, i] <- as.integer(covData[[i]]$X + covData[[i]]$U)
    covData[i] <- list(NULL) #Free memory as we go, keeping list positions
  }
  raw <- list(chr = chrLevels[keys %/% 1e9], pos = as.integer(keys %% 1e9), M = M, Cov = Cov)
  message(length(keys), " CpGs covered in at least one sample")
  saveRDS(raw, "prep-raw/raw-counts.rds", compress = FALSE)
  file.create("prep-raw/prep-complete")
}

# Stage 2: sample QC. Exclusion rule set before looking at any DML results (plan section D3):
#   flag a sample if bisulfite conversion (from CHH) < 99%, or if mapping efficiency, deduplicated read pairs, or CpGs
#   at >= 5x are more than 3 robust SD (MAD) below the median, or duplication is more than 3 robust SD above it.
#   PCA position and correlation with other samples are reported but are not exclusion criteria by themselves,
#   since they can reflect biology. The PC-metric correlations show whether PCA outliers are explained by technical metrics

if (mode == "qc") {
  qcDir <- "sample-QC"
  dir.create(qcDir, showWarnings = FALSE)
  sampleNumber <- function(x) as.integer(str_match(x, "zr3644_([0-9]+)_")[, 2])

  trimming <- read_tsv("../Haws_01-trimgalore/multiqc_data_1/multiqc_cutadapt.txt", show_col_types = FALSE) %>%
    filter(str_detect(Sample, "_R1$")) %>%
    transmute(sample_number = sampleNumber(Sample), rawReadPairs = r_processed, percentBasesTrimmed = percent_trimmed)
  alignment <- read_tsv("../Haws_03-bismark-roslin/multiqc_data/multiqc_bismark_alignment.txt", show_col_types = FALSE) %>%
    transmute(sample_number = sampleNumber(Sample), bismarkInputPairs = total_reads, alignedPairs = aligned_reads,
              percentAligned = percent_aligned, percentAmbiguous = 100 * ambig_reads / total_reads)
  dedup <- read_tsv("../Haws_03-bismark-roslin/multiqc_data/multiqc_bismark_dedup.txt", show_col_types = FALSE) %>%
    transmute(sample_number = sampleNumber(Sample), dedupPairs = dedup_reads, percentDuplicated = dup_reads_percent)
  methExtract <- read_tsv("../Haws_03-bismark-roslin/multiqc_data/multiqc_bismark_methextract.txt", show_col_types = FALSE) %>%
    transmute(sample_number = sampleNumber(Sample),
              percentCpGMethBismark = 100 * meth_cpg / (meth_cpg + unmeth_cpg),
              percentCHGMeth = 100 * meth_chg / (meth_chg + unmeth_chg),
              percentCHHMeth = 100 * meth_chh / (meth_chh + unmeth_chh),
              conversionRate = 100 - percentCHHMeth) #Non-CpG methylation is near 0 in oysters, so CHH methylation estimates non-conversion

  raw <- readRDS("prep-raw/raw-counts.rds")
  nuclear <- raw$chr != mitoChr
  mito <- !nuclear
  coverage <- map_dfr(1:24, function(i) {
    cov <- raw$Cov[, i]
    m <- raw$M[, i]
    tibble(sample_number = i,
           CpGs1x = sum(cov > 0),
           CpGs5x = sum(cov >= 5),
           medianCovCovered = median(cov[nuclear & cov > 0]),
           meanCovCovered = mean(cov[nuclear & cov > 0]),
           hiPerc99.9Cov = unname(quantile(cov[nuclear & cov >= 2], 0.999)), #Same base set as methylKit methRead(mincov = 2)
           percentCpGMethNuclear = 100 * sum(m[nuclear]) / sum(cov[nuclear]),
           mitoCpGs5x = sum(cov[mito] >= 5),
           mitoMeanCov = mean(cov[mito]),
           percentCpGMethMito = 100 * sum(m[mito]) / sum(cov[mito]))
  })

  # PCA on nuclear CpGs with >= 5x and below the upper percentile in every sample
  hiCut <- coverage$hiPerc99.9Cov
  pcaRows <- which(nuclear & rowSums(raw$Cov >= 5) == 24 & rowSums(sweep(raw$Cov, 2, hiCut, ">")) == 0)
  beta <- raw$M[pcaRows, ] / raw$Cov[pcaRows, ]
  rm(raw); invisible(gc())
  beta <- beta[matrixStats::rowVars(beta) > 0, ]
  message(nrow(beta), " CpGs used for PCA")
  pca <- prcomp(t(beta), center = TRUE, scale. = FALSE)
  varExplained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  sampleCor <- cor(beta)
  rm(beta); invisible(gc())

  qcTable <- sampleMetadata %>%
    dplyr::select(sample_number, Library_name, sampleID, ploidy, pH) %>% #SeqID is dropped with column 1 above
    left_join(trimming, by = "sample_number") %>%
    left_join(alignment, by = "sample_number") %>%
    left_join(dedup, by = "sample_number") %>%
    left_join(methExtract, by = "sample_number") %>%
    left_join(coverage, by = "sample_number") %>%
    mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],
           meanCorOthers = (rowSums(sampleCor) - 1) / 23)

  robustZ <- function(x) {
    s <- mad(x)
    if (s == 0) return(rep(0, length(x)))
    (x - median(x)) / s
  }
  flags <- qcTable %>%
    transmute(sample_number, sampleID,
              zPercentAligned = robustZ(percentAligned),
              zDedupPairs = robustZ(dedupPairs),
              zCpGs5x = robustZ(CpGs5x),
              zPercentDuplicated = robustZ(percentDuplicated),
              zMeanCorOthers = robustZ(meanCorOthers),
              zPC1 = robustZ(PC1), zPC2 = robustZ(PC2),
              flagConversion = qcTable$conversionRate < 99,
              flagAligned = zPercentAligned < -3,
              flagDedupPairs = zDedupPairs < -3,
              flagCpGs5x = zCpGs5x < -3,
              flagDuplicated = zPercentDuplicated > 3,
              exclude = flagConversion | flagAligned | flagDedupPairs | flagCpGs5x | flagDuplicated,
              noteCorrelation = zMeanCorOthers < -3, #Reported only
              notePCA = abs(zPC1) > 3 | abs(zPC2) > 3) #Reported only

  technicalMetrics <- c("rawReadPairs", "dedupPairs", "percentAligned", "percentDuplicated", "conversionRate",
                        "CpGs5x", "medianCovCovered", "percentCpGMethNuclear", "percentBasesTrimmed")
  pcMetricCor <- expand_grid(PC = paste0("PC", 1:4), metric = technicalMetrics) %>%
    mutate(spearmanRho = map2_dbl(PC, metric, ~ cor(qcTable[[.x]], qcTable[[.y]], method = "spearman")),
           p = map2_dbl(PC, metric, ~ suppressWarnings(cor.test(qcTable[[.x]], qcTable[[.y]], method = "spearman"))$p.value),
           percentVariance = varExplained[as.integer(sub("PC", "", PC))])

  write_csv(qcTable, file.path(qcDir, "sample-QC-table.csv"))
  write_csv(flags, file.path(qcDir, "sample-QC-flags.csv"))
  write_csv(pcMetricCor, file.path(qcDir, "PC-QC-metric-correlations.csv"))
  write_csv(tibble(PC = seq_along(varExplained), percentVariance = varExplained, CpGs = nrow(pca$rotation)),
            file.path(qcDir, "PCA-variance-explained.csv"))
  write.csv(sampleCor, file.path(qcDir, "sample-correlation-matrix.csv"))

  plotData <- qcTable %>% mutate(group = paste0(ploidy, " ", pH, " pH"), outlier = sample_number %in% outliers)
  for (pcs in list(c(1, 2), c(3, 4))) {
    p <- ggplot(plotData, aes(.data[[paste0("PC", pcs[1])]], .data[[paste0("PC", pcs[2])]], color = group)) +
      geom_point(aes(shape = outlier), size = 3) +
      geom_text(aes(label = sampleID), vjust = -0.9, size = 3, show.legend = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("", "2H-3 / 3H-2")) +
      labs(x = sprintf("PC%d (%.1f%%)", pcs[1], varExplained[pcs[1]]),
           y = sprintf("PC%d (%.1f%%)", pcs[2], varExplained[pcs[2]]),
           color = NULL, shape = NULL, subtitle = paste(format(nrow(pca$rotation), big.mark = ","), "nuclear CpGs, >= 5x in all samples")) +
      theme_bw()
    ggsave(file.path(qcDir, sprintf("PCA-PC%d-PC%d.png", pcs[1], pcs[2])), p, width = 7, height = 5.5, dpi = 150)
  }
  p <- plotData %>%
    dplyr::select(sampleID, outlier, PC1, PC2, all_of(technicalMetrics)) %>%
    pivot_longer(all_of(technicalMetrics), names_to = "metric") %>%
    pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
    ggplot(aes(value, score, color = outlier)) +
    geom_point() +
    facet_grid(PC ~ metric, scales = "free_x") +
    scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "red"), labels = c("", "2H-3 / 3H-2")) +
    labs(x = NULL, y = "PC score", color = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(qcDir, "PC-vs-QC-metrics.png"), p, width = 16, height = 5, dpi = 150)

  print(as.data.frame(flags %>% dplyr::select(sampleID, starts_with("flag"), exclude, starts_with("note"))))
  print(as.data.frame(pcMetricCor %>% filter(PC %in% c("PC1", "PC2")) %>% arrange(p) %>% head(10)))
}

# Stage 3: filters for one setting

if (mode == "prep") {
  dir.create(prepDir, showWarnings = FALSE)
  raw <- readRDS("prep-raw/raw-counts.rds")
  M <- raw$M[, keep]
  Cov <- raw$Cov[, keep]
  nuclear <- raw$chr != mitoChr
  filterCounts <- tibble(step = "covered in any sample", CpGs = nrow(Cov), mitoCpGs = sum(!nuclear))

  lowCov <- Cov < loCov #Per sample: treat coverage below LO_COV as missing
  if (hiPercSetting != "none") {
    hiCut <- apply(Cov, 2, function(x) quantile(x[nuclear & x >= 2], as.numeric(hiPercSetting) / 100)) #Same base set as methylKit methRead(mincov = 2)
    highCov <- sweep(Cov, 2, hiCut, ">") & nuclear #Per sample: drop nuclear CpGs above the percentile; mito exempt
    message("Upper coverage cutoffs: ", paste(round(hiCut), collapse = ", "))
  } else {
    highCov <- FALSE
  }
  drop <- lowCov | highCov
  M[drop] <- 0L
  Cov[drop] <- 0L

  design <- makeDesign(sampleMetadata[keep, ], 0)
  if (presence == "all") {
    keepRows <- rowSums(Cov > 0) == ncol(Cov)
  } else {
    minPerCell <- as.integer(sub("cell", "", presence))
    keepRows <- Reduce(`&`, lapply(split(seq_len(ncol(Cov)), design$cell), function(cols) {
      rowSums(Cov[, cols, drop = FALSE] > 0) >= min(minPerCell, length(cols)) #A group with fewer samples (samples dropped) needs all of them
    }))
  }
  filterCounts <- bind_rows(filterCounts,
                            tibble(step = paste0("after LO_COV = ", loCov, ", HI_PERC = ", hiPercSetting, ", PRESENCE = ", presence),
                                   CpGs = sum(keepRows), mitoCpGs = sum(keepRows & !nuclear)))
  if (minMethSetting != "none") {
    minMeth <- as.numeric(minMethSetting)
    beta <- M / Cov
    beta[Cov == 0] <- NA
    meanMeth <- 100 * rowMeans(beta, na.rm = TRUE) #Mean over samples, ignoring labels
    rm(beta); invisible(gc())
    keepRows <- keepRows & ((meanMeth >= minMeth & meanMeth <= 100 - minMeth) | !nuclear) #Mito exempt (~2% methylation everywhere), as for HI_PERC
    filterCounts <- bind_rows(filterCounts,
                              tibble(step = paste0("after MIN_METH = ", minMethSetting, " (mean methylation ", minMeth, "-", 100 - minMeth, "%; mito exempt)"),
                                     CpGs = sum(keepRows), mitoCpGs = sum(keepRows & !nuclear)))
  }
  if (snpFilter == "any") {
    cpgSNP <- readRDS("prep-bssnper/CpG-SNP-flags.rds")
    stopifnot(identical(length(cpgSNP$nSamplesSNP), length(raw$pos)))
    keepRows <- keepRows & cpgSNP$nSamplesSNP == 0
    filterCounts <- bind_rows(filterCounts,
                              tibble(step = "after SNP_FILTER = any (no PASS SNP at the C or G in any sample)",
                                     CpGs = sum(keepRows), mitoCpGs = sum(keepRows & !nuclear)))
  }
  filtered <- list(chr = raw$chr[keepRows], pos = raw$pos[keepRows], M = M[keepRows, ], Cov = Cov[keepRows, ])
  print(as.data.frame(filterCounts))
  write_tsv(filterCounts, file.path(prepDir, "CpG-filter-counts.tsv"))
  saveRDS(filtered, file.path(prepDir, "filtered-counts.rds"), compress = FALSE)
  file.create(file.path(prepDir, "prep-complete"))
}

# Stage 4: fit one block of CpGs (all CpGs with the default NCHUNK = 1). DSS estimates the dispersion prior within each
# call, so if NCHUNK > 1 keep blocks large. p-values are kept; FDR is recomputed over all CpGs in the summary stage

modelTerms <- c(ploidy = "ploidyEffect", pH = "pHEffect", interaction = "ploidyEffect:pHEffect")
if (model == "additive") modelTerms <- modelTerms[c("ploidy", "pH")]

if (mode == "fit") {
  suppressPackageStartupMessages(library(DSS))
  taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  stopifnot(taskID %in% seq_len(nChunk))
  chunkDir <- file.path(runDir, "chunks")
  dir.create(chunkDir, recursive = TRUE, showWarnings = FALSE)

  filtered <- readRDS(file.path(prepDir, "filtered-counts.rds"))
  rows <- if (nChunk == 1) seq_along(filtered$pos) else
    which(cut(seq_along(filtered$pos), nChunk, labels = FALSE) == taskID) #Contiguous blocks in genome order. cut() needs at least 2 intervals
  design <- makeDesign(sampleMetadata[keep, ], perm)
  message("Task ", taskID, " of ", nChunk, ": ", length(rows), " CpGs, ", runDir)

  BSobj <- bsseq::BSseq(chr = filtered$chr[rows], pos = filtered$pos[rows],
                        M = filtered$M[rows, ], Cov = filtered$Cov[rows, ],
                        sampleNames = rownames(design))
  formula <- if (model == "interaction") ~ ploidyEffect * pHEffect else ~ ploidyEffect + pHEffect
  DMLfit <- DMLfit.multiFactor(BSobj, design = design[, c("ploidyEffect", "pHEffect")], formula = formula)
  stopifnot(all(modelTerms %in% colnames(DMLfit$X)))

  result <- data.frame(chr = filtered$chr[rows], pos = filtered$pos[rows],
                       nSamples = rowSums(filtered$Cov[rows, ] > 0),
                       effectSizes(filtered$M[rows, ], filtered$Cov[rows, ], design))
  for (term in names(modelTerms)) {
    test <- DMLtest.multiFactor(DMLfit, coef = modelTerms[[term]])
    stopifnot(identical(as.character(test$chr), result$chr), identical(as.integer(test$pos), result$pos))
    result[[paste0("stat.", term)]] <- test$stat
    result[[paste0("p.", term)]] <- test$pvals
  }
  saveRDS(result, file.path(chunkDir, sprintf("chunk-%03d.rds", taskID)), compress = FALSE)
}

# Stage 5: combine blocks, FDR over all CpGs, DML and DMR calls, counts

if (mode == "summary") {
  suppressPackageStartupMessages(library(DSS))
  chunkFiles <- list.files(file.path(runDir, "chunks"), pattern = "^chunk-[0-9]+\\.rds$", full.names = TRUE)
  stopifnot(length(chunkFiles) == nChunk)
  results <- rbindlist(lapply(chunkFiles, readRDS))
  setorder(results, chr, pos)
  dir.create(file.path(runDir, "rds"), showWarnings = FALSE)

  diffColumn <- c(ploidy = "diffPloidy", pH = "diffpH", interaction = "diffInteraction")
  counts <- list()
  for (term in names(modelTerms)) {
    p <- results[[paste0("p.", term)]]
    results[[paste0("fdr.", term)]] <- p.adjust(p, method = "BH")

    png(file.path(runDir, paste0("pvalue-histogram-", term, ".png")), width = 700, height = 500)
    hist(p, breaks = 100, main = paste0(term, " (", runDir, ")"), xlab = "p-value", col = "grey70", border = NA)
    dev.off()

    for (fdrCut in c(0.05, 0.01)) for (minDiff in c(0, 10, 25, 50)) {
      isDML <- results[[paste0("fdr.", term)]] < fdrCut & abs(results[[diffColumn[[term]]]]) >= minDiff
      counts[[length(counts) + 1]] <- tibble(term = term, fdr = fdrCut, minDiff = minDiff,
                                             CpGsTested = nrow(results), DML = sum(isDML, na.rm = TRUE),
                                             mitoCpGsTested = sum(results$chr == mitoChr),
                                             mitoDML = sum(isDML & results$chr == mitoChr, na.rm = TRUE),
                                             hypermethylated = sum(isDML & results[[diffColumn[[term]]]] > 0, na.rm = TRUE))
    }

    if (term != "interaction") { #Interaction is a nuisance term: counts only, no DML/DMR files
      outCols <- c("chr", "pos", "nSamples", "meanMeth", diffColumn[[term]], paste0(c("stat.", "p.", "fdr."), term))
      for (fdrCut in c(0.05, 0.01)) for (minDiff in c(0, 25, 50)) {
        dml <- results[results[[paste0("fdr.", term)]] < fdrCut & abs(results[[diffColumn[[term]]]]) >= minDiff, ..outCols]
        dml <- cbind(dml[, .(chr, start = pos, end = pos)], dml[, -c("chr", "pos")])
        write_tsv(dml, file.path(runDir, sprintf("DML-%s-fdr%s-diff%d.bed", term, fdrCut, minDiff)))
      }
      dmlInput <- data.frame(chr = results$chr, pos = results$pos, stat = results[[paste0("stat.", term)]],
                             pvals = p, fdrs = results[[paste0("fdr.", term)]])
      class(dmlInput) <- c("DMLtest.multiFactor", "data.frame") #callDMR checks the class to use p-values only (no delta)
      dmr <- tryCatch(callDMR(dmlInput, p.threshold = 1e-5, minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5),
                      error = function(e) { message("callDMR failed for ", term, ": ", conditionMessage(e)); NA })
      nDMR <- if (is.null(dmr)) 0L else if (identical(dmr, NA)) NA_integer_ else nrow(dmr) #callDMR returns NULL when no DMR pass
      if (identical(dmr, NA)) dmr <- NULL
      if (!is.null(dmr)) write_tsv(dmr, file.path(runDir, paste0("DMR-", term, "-p1e-5.tsv")))
      counts[[length(counts) + 1]] <- tibble(term = term, fdr = NA, minDiff = NA, CpGsTested = nrow(results), DMR = nDMR)
    }
  }
  saveRDS(results, file.path(runDir, "rds", "all-CpG-results.rds"), compress = FALSE) #Every CpG tested, for other thresholds and the concordance step

  counts <- bind_rows(counts) %>% mutate(setting = settingTag, model = model, perm = perm, .before = 1)
  write_csv(counts, file.path(runDir, "DSS-counts-long.csv"))

  if ("interaction" %in% names(modelTerms)) { #How many main-effect DML also have a strong interaction (diagnostic)
    overlap <- map_dfr(c("ploidy", "pH"), function(term) {
      isDML <- results[[paste0("fdr.", term)]] < 0.05
      tibble(term = term, DMLfdr0.05 = sum(isDML), withInteractionFdr0.05 = sum(isDML & results$fdr.interaction < 0.05))
    })
    write_csv(overlap, file.path(runDir, "main-effect-DML-with-interaction.csv"))
  }

  summaryDML <- counts %>% #Same layout as the methylKit DML-summary-table.csv (FDR < 0.01; 25% and 50% difference)
    filter(term %in% c("ploidy", "pH"), fdr == 0.01, minDiff %in% c(25, 50)) %>%
    transmute(samples = .env$samples, presence = .env$presence,
              column = paste0(term, "-DML-", minDiff), DML) %>%
    pivot_wider(names_from = column, values_from = DML) %>%
    dplyr::select(samples, presence, any_of(c("ploidy-DML-25", "ploidy-DML-50", "pH-DML-25", "pH-DML-50")))
  write.csv(summaryDML, file.path(runDir, "DML-summary-table.csv"), row.names = FALSE, quote = FALSE)
  print(as.data.frame(counts))
}

# Stage 6: sample-set sensitivity. For each DML set from All (SAMPLES ignored here; other settings as set), count how
# many CpGs were tested in each sensitivity run, are DML at the same thresholds, or have p < 0.05 with the same sign

if (mode == "compare") {
  runFor <- function(s) paste0("DSS-", filterTag, "-", s, "-", model)
  readResults <- function(s) readRDS(file.path(runFor(s), "rds", "all-CpG-results.rds"))
  primary <- readResults("All")
  diffColumn <- c(ploidy = "diffPloidy", pH = "diffpH")
  thresholds <- tibble(fdr = c(0.05, 0.01, 0.01), minDiff = c(0, 0, 25))

  comparison <- map_dfr(c("Drop3H2", "DropPC2", "OutRM"), function(s) {
    other <- readResults(s)
    map_dfr(names(diffColumn), function(term) {
      matched <- match(paste(primary$chr, primary$pos), paste(other$chr, other$pos))
      pmap_dfr(thresholds, function(fdr, minDiff) {
        isDML <- function(r, rows = TRUE) r[[paste0("fdr.", term)]][rows] < fdr & abs(r[[diffColumn[[term]]]][rows]) >= minDiff
        primaryDML <- which(isDML(primary))
        m <- matched[primaryDML]
        tested <- !is.na(m)
        sameSign <- sign(primary[[paste0("stat.", term)]][primaryDML]) == sign(other[[paste0("stat.", term)]][m])
        tibble(sampleSet = s, term = term, fdr = fdr, minDiff = minDiff,
               DMLAll = length(primaryDML),
               DMLSampleSet = sum(isDML(other), na.rm = TRUE),
               AllDMLTested = sum(tested),
               AllDMLStillDML = sum(isDML(other, m[tested]), na.rm = TRUE),
               AllDMLp0.05SameSign = sum(other[[paste0("p.", term)]][m[tested]] < 0.05 & sameSign[tested], na.rm = TRUE))
      })
    })
  })
  outFile <- paste0("sample-set-sensitivity-", filterTag, "-", model, ".csv")
  write_csv(comparison, outFile)
  print(as.data.frame(comparison))
}

# Stage 7: permutation null. Labels shuffled within strata (PERM = seed), so any DML found are false positives.
# Mean permuted count / observed count estimates the fraction of observed DML that are false

if (mode == "permsummary") {
  observedDir <- paste0("DSS-", settingTag, "-", model)
  permDirs <- list.files(".", pattern = paste0("^", gsub("\\.", "\\\\.", observedDir), "-perm[0-9]+$"))
  permDirs <- permDirs[file.exists(file.path(permDirs, "DSS-counts-long.csv"))]
  message(length(permDirs), " permutations with results")
  readCounts <- function(d) read_csv(file.path(d, "DSS-counts-long.csv"), show_col_types = FALSE) %>% filter(!is.na(fdr))
  observed <- readCounts(observedDir) %>% dplyr::select(term, fdr, minDiff, observed = DML)
  permuted <- map_dfr(permDirs, readCounts) %>% dplyr::select(perm, term, fdr, minDiff, DML)
  write_csv(permuted, paste0("permutation-counts-", settingTag, "-", model, ".csv"))

  permSummary <- permuted %>%
    group_by(term, fdr, minDiff) %>%
    summarize(nPerm = n(), permMean = mean(DML), permMedian = median(DML), permMax = max(DML), .groups = "drop") %>%
    left_join(observed, by = c("term", "fdr", "minDiff")) %>%
    left_join(permuted %>% left_join(observed, by = c("term", "fdr", "minDiff")) %>%
                group_by(term, fdr, minDiff) %>% summarize(nPermAtLeastObserved = sum(DML >= observed), .groups = "drop"),
              by = c("term", "fdr", "minDiff")) %>%
    mutate(empiricalFDR = ifelse(observed > 0, pmin(1, permMean / observed), NA),
           empiricalP = (nPermAtLeastObserved + 1) / (nPerm + 1)) %>% #Chance of a permutation giving at least as many DML
    arrange(factor(term, levels = c("ploidy", "pH", "interaction")), desc(fdr), minDiff)
  write_csv(permSummary, paste0("permutation-summary-", settingTag, "-", model, ".csv"))
  print(as.data.frame(permSummary))
}

# Stage 8: BS-SNPer SNPs at CpGs. Uses the per-sample VCFs from the earlier BS-SNPer run (code/Haws/05-BS-SNPer.ipynb;
# same bismark-2 BAMs and Roslin + mito genome as the coverage files; --mincover 5, other settings default), downloaded to
# prep-bssnper/. Only PASS calls are used (Low calls are mostly 1-2 reads). A CpG is marked for a sample if that sample has
# a SNP at the C (pos) or at the G (pos + 1): either makes the site look unmethylated in bisulfite data

vcfFile <- function(i) file.path("prep-bssnper", paste0("zr3644_", i, "_R1_val_1_val_1_val_1_bismark_bt2_pe.SNP-results.vcf"))

if (mode == "snp-prep") {
  stopifnot(all(file.exists(vcfFile(1:24))))
  raw <- readRDS("prep-raw/raw-counts.rds")
  chrLevels <- sort(unique(raw$chr))
  cpgKey <- match(raw$chr, chrLevels) * 1e9 + raw$pos
  snpMatrix <- matrix(FALSE, nrow = length(cpgKey), ncol = 24, dimnames = list(NULL, sampleMetadata$sampleID))
  altFreq <- matrix(0, nrow = length(cpgKey), ncol = 24, dimnames = list(NULL, sampleMetadata$sampleID)) #Alt allele frequency at the C or G (larger of the two)
  snpNonBS <- snpMatrix #Same, using only SNPs bisulfite conversion cannot mimic: not C>T at the C, not G>A at the G
  altFreqNonBS <- altFreq
  perSample <- list()
  for (i in 1:24) {
    message("Reading ", basename(vcfFile(i)))
    vcf <- fread(cmd = paste("grep -v '^##'", shQuote(vcfFile(i))), sep = "\t", header = TRUE,
                 select = c(1, 2, 4, 5, 7, 10), col.names = c("chr", "pos", "ref", "alt", "filter", "sample"))
    perSample[[i]] <- tibble(sample_number = i, SNPsAll = nrow(vcf), SNPsPASS = sum(vcf$filter == "PASS"),
                             CTorGA = sum(vcf$filter == "PASS" & ((vcf$ref == "C" & vcf$alt == "T") | (vcf$ref == "G" & vcf$alt == "A"))),
                             heterozygous = sum(vcf$filter == "PASS" & startsWith(vcf$sample, "0/1")))
    vcf <- vcf[filter == "PASS" & chr %in% chrLevels]
    vcf[, alfr := as.numeric(sub(".*,", "", sub(".*:", "", sample)))] #Last FORMAT field ALFR = "ref,alt"
    key <- match(vcf$chr, chrLevels) * 1e9 + vcf$pos
    for (offset in c(0, 1)) { #SNP at the C (pos) or the G (pos + 1)
      rowIndex <- match(key - offset, cpgKey)
      hit <- !is.na(rowIndex)
      snpMatrix[rowIndex[hit], i] <- TRUE
      altFreq[rowIndex[hit], i] <- pmax(altFreq[rowIndex[hit], i], vcf$alfr[hit])
      bsLike <- if (offset == 0) vcf$ref == "C" & vcf$alt == "T" else vcf$ref == "G" & vcf$alt == "A"
      hitNonBS <- hit & !bsLike
      snpNonBS[rowIndex[hitNonBS], i] <- TRUE
      altFreqNonBS[rowIndex[hitNonBS], i] <- pmax(altFreqNonBS[rowIndex[hitNonBS], i], vcf$alfr[hitNonBS])
    }
  }
  nSamplesSNP <- rowSums(snpMatrix)
  saveRDS(list(nSamplesSNP = nSamplesSNP), "prep-bssnper/CpG-SNP-flags.rds", compress = FALSE)
  snpRows <- which(nSamplesSNP > 0)
  saveRDS(list(chr = raw$chr[snpRows], pos = raw$pos[snpRows], snp = snpMatrix[snpRows, ], altFreq = altFreq[snpRows, ],
               snpNonBS = snpNonBS[snpRows, ], altFreqNonBS = altFreqNonBS[snpRows, ], Cov = raw$Cov[snpRows, ]),
          "prep-bssnper/CpG-SNP-genotypes.rds", compress = FALSE)

  perSample <- bind_rows(perSample) %>%
    mutate(CpGsWithSNP = colSums(snpMatrix), sampleID = sampleMetadata$sampleID, .after = sample_number)
  dir.create("SNP-check", showWarnings = FALSE)
  write_csv(perSample, "SNP-check/SNPs-per-sample.csv")
  write_csv(tibble(nSamplesWithSNP = 0:24, CpGs = tabulate(nSamplesSNP + 1, nbins = 25)), "SNP-check/CpGs-by-number-of-samples-with-SNP.csv")
  print(as.data.frame(perSample))
  message(sum(nSamplesSNP > 0), " of ", length(nSamplesSNP), " covered CpGs have a PASS SNP at the C or G in at least one sample")
}

if (mode == "snp-summary") {
  dir.create("SNP-check", showWarnings = FALSE)
  raw <- readRDS("prep-raw/raw-counts.rds")
  rawKey <- paste(raw$chr, raw$pos)
  nSamplesSNP <- readRDS("prep-bssnper/CpG-SNP-flags.rds")$nSamplesSNP
  rm(raw); invisible(gc())

  # 1. Are DML more often at SNP CpGs than the CpGs tested? Observed and each permutation
  runDirs <- list.files(".", pattern = "^DSS-.*-interaction(-perm[0-9]+)?$")
  runDirs <- runDirs[file.exists(file.path(runDirs, "rds", "all-CpG-results.rds"))]
  enrichment <- map_dfr(runDirs, function(d) {
    results <- readRDS(file.path(d, "rds", "all-CpG-results.rds"))
    snp <- nSamplesSNP[match(paste(results$chr, results$pos), rawKey)] > 0
    map_dfr(c("ploidy", "pH", "interaction"), function(term) {
      fdr <- results[[paste0("fdr.", term)]]
      p <- results[[paste0("p.", term)]]
      tibble(run = d, setting = sub("-interaction(-perm[0-9]+)?$", "", sub("^DSS-", "", d)),
             perm = as.integer(ifelse(grepl("-perm", d), sub(".*-perm", "", d), "0")), term = term,
             CpGsTested = length(snp), fractionSNPTested = mean(snp),
             DMLfdr0.05 = sum(fdr < 0.05), fractionSNPDML = mean(snp[fdr < 0.05]),
             p0.001 = sum(p < 0.001), fractionSNPp0.001 = mean(snp[p < 0.001]))
    })
  })
  write_csv(enrichment, "SNP-check/SNP-enrichment-among-DML.csv")
  print(as.data.frame(enrichment %>% group_by(setting, term, observed = perm == 0) %>%
                        summarize(runs = n(), fractionSNPTested = mean(fractionSNPTested), DMLfdr0.05 = mean(DMLfdr0.05),
                                  fractionSNPDML = mean(fractionSNPDML, na.rm = TRUE), fractionSNPp0.001 = mean(fractionSNPp0.001, na.rm = TRUE),
                                  .groups = "drop")))

  # 2. Genetic structure: PCA of alt allele frequency at SNP CpGs covered >= 10x in all samples and variable across samples.
  # A SNP absent at a well-covered site is treated as reference. Compared with the methylation PCA from the sample QC.
  # Run twice: all SNPs, and only SNPs bisulfite conversion cannot mimic ("nonBS": not C>T at the C, not G>A at the G).
  # At C>T / G>A sites the genotype call can carry some methylation signal, so the all-SNP comparison may be partly
  # circular; the nonBS set is the check
  geno <- readRDS("prep-bssnper/CpG-SNP-genotypes.rds")
  qc <- read_csv("sample-QC/sample-QC-table.csv", show_col_types = FALSE)
  methCor <- as.matrix(read.csv("sample-QC/sample-correlation-matrix.csv", row.names = 1, check.names = FALSE))
  genotypeStructure <- function(snp, altFreq, label) {
    suffix <- ifelse(label == "all", "", paste0("-", label))
    useRows <- rowSums(geno$Cov >= 10) == 24 & rowSums(snp) >= 2 & rowSums(snp) <= 22
    G <- altFreq[useRows, ]
    message(label, ": ", nrow(G), " SNP CpGs used for the genotype PCA")
    genoPCA <- prcomp(t(G), center = TRUE, scale. = FALSE)
    genoVar <- 100 * genoPCA$sdev^2 / sum(genoPCA$sdev^2)
    genoCor <- cor(G)
    genoScores <- tibble(sampleID = colnames(G), genoPC1 = genoPCA$x[, 1], genoPC2 = genoPCA$x[, 2], genoPC3 = genoPCA$x[, 3],
                         meanGenoCorOthers = (rowSums(genoCor) - 1) / 23) %>%
      left_join(qc %>% dplyr::select(sampleID, ploidy, pH, methPC1 = PC1, methPC2 = PC2, dedupPairs), by = "sampleID")
    write_csv(genoScores, paste0("SNP-check/genotype-PCA-scores", suffix, ".csv"))
    write.csv(genoCor, paste0("SNP-check/genotype-correlation-matrix", suffix, ".csv"))
    pcCor <- expand_grid(genoPC = paste0("genoPC", 1:3), methPC = c("methPC1", "methPC2")) %>%
      mutate(spearmanRho = map2_dbl(genoPC, methPC, ~ cor(genoScores[[.x]], genoScores[[.y]], method = "spearman")),
             p = map2_dbl(genoPC, methPC, ~ suppressWarnings(cor.test(genoScores[[.x]], genoScores[[.y]], method = "spearman"))$p.value),
             genoPercentVariance = genoVar[as.integer(sub("genoPC", "", genoPC))])
    write_csv(pcCor, paste0("SNP-check/genotype-vs-methylation-PC-correlations", suffix, ".csv"))
    pairs <- upper.tri(genoCor) #Is genotype similarity related to methylation similarity across sample pairs?
    pairRho <- cor(genoCor[pairs], methCor[colnames(G), colnames(G)][pairs], method = "spearman")
    message(label, ": Spearman correlation of genotype vs methylation similarity across sample pairs: ", round(pairRho, 3))
    writeLines(paste0("SNPCpGs\t", nrow(G), "\npairs\t", sum(pairs), "\nspearmanRho\t", pairRho),
               paste0("SNP-check/genotype-vs-methylation-similarity", suffix, ".tsv"))
    plotData <- genoScores %>% mutate(group = paste0(ploidy, " ", pH, " pH"),
                                      highlight = sampleID %in% c("3H-2", "2H-1", "2H-2", "3H-5"))
    p <- ggplot(plotData, aes(genoPC1, genoPC2, color = group)) +
      geom_point(aes(shape = highlight), size = 3) +
      geom_text(aes(label = sampleID), vjust = -0.9, size = 3, show.legend = FALSE) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("", "3H-2 / PC2 group")) +
      labs(x = sprintf("genotype PC1 (%.1f%%)", genoVar[1]), y = sprintf("genotype PC2 (%.1f%%)", genoVar[2]), color = NULL, shape = NULL,
           subtitle = paste0(format(nrow(G), big.mark = ","), " SNP CpGs, >= 10x in all samples",
                             ifelse(label == "nonBS", " (no C>T at C / G>A at G)", ""))) +
      theme_bw()
    ggsave(paste0("SNP-check/genotype-PCA", suffix, ".png"), p, width = 7, height = 5.5, dpi = 150)
    print(as.data.frame(genoScores))
    print(as.data.frame(pcCor))
  }
  genotypeStructure(geno$snp, geno$altFreq, "all")
  genotypeStructure(geno$snpNonBS, geno$altFreqNonBS, "nonBS")
}

sessionInfo()
