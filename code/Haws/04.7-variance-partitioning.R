#!/usr/bin/env Rscript
# How much of the methylation variation between oysters is explained by genetic background, ploidy, and pH?
# Plan: 04.7-variance-partitioning-plan.md. Run from analyses/Haws_04.7-variance-partitioning/ by
# 04.7-variance-partitioning-submit.sh. One job; all inputs come from the 04.4 DSS and 04.6 methylKit outputs:
#   ../Haws_04.4-DSS/prep-cov5-hiperc99.9-cell5-meth10-noSNP-All/filtered-counts.rds  CpGs (>= 5x, mean methylation 10-90%,
#                                                                                      no BS-SNPer SNP in any oyster)
#   ../Haws_04.4-DSS/prep-bssnper/CpG-SNP-genotypes.rds                                SNP genotypes at CpGs
#   ../Haws_04.6-methylKit-DMR/count-hiperc99.9-gene-cb3/region-counts.rds             reads per gene body
#
# For each sample set (SAMPLE_SETS) and feature set (CpGs, genes):
#   1. Genotype PCs from the bisulfite-safe SNPs (no C>T at the C, no G>A at the G), recomputed within the sample set
#   2. Whole-methylome tests on Euclidean distance between oysters (methylation proportions):
#      PERMANOVA (vegan::adonis2, marginal terms), variance partitioning (vegan::varpart, adjusted R2) of genotype PCs vs
#      treatment (ploidy * pH), and Mantel / partial Mantel tests of genotype distance vs methylation distance
#   3. Per-feature linear models of methylation proportion: R2 of treatment (ploidy * pH, 3 df), genotype (K PCs), and both.
#      Unique fractions = R2(both) - R2(other). The nulls: the same statistics with treatment labels shuffled (same scheme
#      as 04.4, so ploidy within pH, then pH within ploidy), and with genotype PCs shuffled between oysters
#
# Environment variables (defaults in brackets):
#   SAMPLE_SETS [All DropGeno4]  Space-separated. DropGeno4 drops the four genetic outliers (2H-1, 2H-2, 3H-2, 3H-5)
#   GENO_K      [3]              Genotype PCs in the models. Also run with 1 and 5 as a sensitivity check
#   NPERM       [9999]           Permutations for PERMANOVA and Mantel tests
#   NNULL       [100]            Label / genotype shuffles for the per-feature nulls

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(vegan)
})
options(warn = 1)

sampleSets <- strsplit(Sys.getenv("SAMPLE_SETS", "All DropGeno4"), " ")[[1]]
genoK <- as.integer(Sys.getenv("GENO_K", "3"))
nPerm <- as.integer(Sys.getenv("NPERM", "9999"))
nNull <- as.integer(Sys.getenv("NNULL", "100"))
kSettings <- sort(unique(c(1, genoK, 5)))
dropSets <- list(All = integer(0), DropGeno4 = c(1, 2, 14, 17))
stopifnot(all(sampleSets %in% names(dropSets)), genoK >= 1)

dssDir <- "../Haws_04.4-DSS"
mkDir <- "../Haws_04.6-methylKit-DMR"
dataDir <- "../../data/Haws"
mitoChr <- "NC_001276.1"
for (d in c("tables", "figures")) dir.create(d, showWarnings = FALSE)

sampleMetadata <- read.csv(file.path(dataDir, "sample_metadata.csv")) %>%
  dplyr::select(-c(1)) %>%
  arrange(sample_number) %>%
  mutate(sampleID = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                      "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                      "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                      "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6")) %>%
  dplyr::rename(pH = ph) %>%
  mutate(ploidy = factor(ploidy, levels = c("2N", "3N")), pH = factor(pH, levels = c("high", "low")))
stopifnot(identical(sampleMetadata$sample_number, 1:24))

# Label shuffle from makeDesign() in 04.4-DSS-slurm.R: ploidy within pH, then pH within the new ploidy labels
permuteLabels <- function(metadata, seed) {
  set.seed(seed)
  metadata %>%
    group_by(pH) %>% mutate(ploidy = sample(ploidy)) %>% ungroup() %>%
    group_by(ploidy) %>% mutate(pH = sample(pH)) %>% ungroup()
}

# Treatment design with effect coding (+1/2 vs -1/2), as in 04.4-DSS-slurm.R: ploidy and pH are each averaged over the
# other factor, and every term (including the interaction, a nuisance term) can be tested given the others
treatmentDesign <- function(meta) {
  ploidyEffect <- ifelse(meta$ploidy == "3N", 0.5, -0.5)
  pHEffect <- ifelse(meta$pH == "low", 0.5, -0.5)
  data.frame(ploidy = ploidyEffect, pH = pHEffect, ploidyXpH = ploidyEffect * pHEffect)
}

# Methylation proportions, features x oysters. Only features with data in every oyster of the set, nuclear only
cpgCounts <- readRDS(file.path(dssDir, "prep-cov5-hiperc99.9-cell5-meth10-noSNP-All", "filtered-counts.rds"))
cpgBeta <- cpgCounts$M / cpgCounts$Cov
cpgBeta[cpgCounts$Cov == 0] <- NA
cpgBeta <- cpgBeta[cpgCounts$chr != mitoChr, ]
colnames(cpgBeta) <- sampleMetadata$sampleID
rm(cpgCounts); invisible(gc())

geneCounts <- readRDS(file.path(mkDir, "count-hiperc99.9-gene-cb3", "region-counts.rds"))
geneData <- lapply(geneCounts, function(x) { d <- methylKit::getData(x); data.table(region = paste(d$chr, d$start, d$end), chr = d$chr, cov = d$coverage, beta = d$numCs / d$coverage) })
stopifnot(identical(vapply(geneCounts, function(x) x@sample.id, ""), sampleMetadata$sampleID))
geneKeys <- Reduce(union, lapply(geneData, `[[`, "region"))
geneBeta <- sapply(geneData, function(d) { b <- d$beta; b[d$cov < 10] <- NA; b[match(geneKeys, d$region)] }) #>= 10 reads, as in 04.6
dimnames(geneBeta) <- list(geneKeys, sampleMetadata$sampleID) #Row names are "chr start end"
geneBeta <- geneBeta[!startsWith(geneKeys, paste0(mitoChr, " ")), ]
rm(geneCounts, geneData); invisible(gc())
features <- list(CpG = cpgBeta, gene = geneBeta)

# Genotype: alt allele frequency at bisulfite-safe SNP CpGs, as in the 04.4 snp-summary stage, but recomputed per sample set
geno <- readRDS(file.path(dssDir, "prep-bssnper", "CpG-SNP-genotypes.rds"))
genotypeFor <- function(cols) {
  snp <- geno$snpNonBS[, cols]
  useRows <- rowSums(geno$Cov[, cols] >= 10) == length(cols) & rowSums(snp) >= 2 & rowSums(snp) <= length(cols) - 2
  G <- geno$altFreqNonBS[useRows, cols]
  pca <- prcomp(t(G), center = TRUE, scale. = FALSE)
  list(scores = pca$x, varExplained = 100 * pca$sdev^2 / sum(pca$sdev^2), nSNP = nrow(G),
       dist = dist(t(G))) #Euclidean genotype distance between oysters
}

# Per-feature R2 for a design matrix (with intercept), computed for all features at once. Y is oysters x features
r2 <- function(Y, X) {
  fitted <- X %*% solve(crossprod(X), crossprod(X, Y))
  tss <- colSums(sweep(Y, 2, colMeans(Y))^2)
  1 - colSums((Y - fitted)^2) / tss
}
partitionFeatures <- function(Y, meta, genoScores, k) {
  treat <- cbind(1, as.matrix(treatmentDesign(meta)))
  genoX <- cbind(1, scale(genoScores[, 1:k, drop = FALSE]))
  both <- cbind(treat, genoX[, -1, drop = FALSE])
  rT <- r2(Y, treat); rG <- r2(Y, genoX); rTG <- r2(Y, both)
  n <- nrow(Y)
  fUniqueG <- ((rTG - rT) / k) / ((1 - rTG) / (n - ncol(both)))
  fUniqueT <- ((rTG - rG) / 3) / ((1 - rTG) / (n - ncol(both)))
  tibble(R2treatment = rT, R2genotype = rG, R2both = rTG,
         uniqueTreatment = rTG - rG, uniqueGenotype = rTG - rT, shared = rT + rG - rTG,
         pUniqueGenotype = pf(fUniqueG, k, n - ncol(both), lower.tail = FALSE),
         pUniqueTreatment = pf(fUniqueT, 3, n - ncol(both), lower.tail = FALSE))
}
summarizeFeatures <- function(p) {
  p %>% summarize(features = n(),
                  across(c(R2treatment, R2genotype, uniqueTreatment, uniqueGenotype, shared), median, .names = "median_{.col}"),
                  fracGenotypeP0.01 = mean(pUniqueGenotype < 0.01), fracTreatmentP0.01 = mean(pUniqueTreatment < 0.01))
}

wholeResults <- list(); featureSummaries <- list(); nullSummaries <- list(); genotypeInfo <- list()
for (s in sampleSets) {
  keep <- setdiff(1:24, dropSets[[s]])
  meta <- sampleMetadata[keep, ]
  g <- genotypeFor(meta$sampleID)
  genotypeInfo[[s]] <- tibble(sampleSet = s, SNPCpGs = g$nSNP, PC = 1:5, percentVariance = g$varExplained[1:5])
  write_csv(as_tibble(g$scores[, 1:5], rownames = "sampleID") %>% left_join(meta %>% dplyr::select(sampleID, ploidy, pH), by = "sampleID"),
            file.path("tables", paste0("genotype-PCs-", s, ".csv")))

  for (fs in names(features)) {
    Y <- features[[fs]][, meta$sampleID]
    Y <- t(Y[rowSums(is.na(Y)) == 0 & matrixStats::rowVars(Y, na.rm = TRUE) > 0, ]) #Oysters x features
    message(s, " ", fs, ": ", ncol(Y), " features, ", nrow(Y), " oysters")
    methDist <- dist(Y)

    # Whole-methylome tests
    for (k in kSettings) {
      genoDF <- as.data.frame(scale(g$scores[, 1:k, drop = FALSE]))
      env <- bind_cols(treatmentDesign(meta), genoDF)
      set.seed(1)
      perma <- adonis2(as.formula(paste("methDist ~ ploidy + pH + ploidyXpH +", paste(names(genoDF), collapse = " + "))),
                       data = env, by = "margin", permutations = nPerm)
      vp <- varpart(methDist, as.matrix(genoDF), as.matrix(treatmentDesign(meta)))
      fr <- vp$part$indfract
      wholeResults[[length(wholeResults) + 1]] <- bind_rows(
        tibble(test = "PERMANOVA (marginal)", term = rownames(perma), R2 = perma$R2, F = perma$F, p = perma$`Pr(>F)`),
        tibble(test = "varpart (adjusted R2)", term = c("genotype only [a]", "shared [b]", "treatment only [c]", "residual [d]"),
               R2 = fr$Adj.R.squared, F = NA, p = NA)) %>%
        mutate(sampleSet = s, features = fs, nFeatures = ncol(Y), genoK = k, .before = 1)
    }
    design <- dist(model.matrix(~ ploidy + pH, meta)[, -1]) #Treatment difference between oysters
    set.seed(1)
    mt <- mantel(g$dist, methDist, method = "spearman", permutations = nPerm)
    pm <- mantel.partial(g$dist, methDist, design, method = "spearman", permutations = nPerm)
    wholeResults[[length(wholeResults) + 1]] <- tibble(sampleSet = s, features = fs, nFeatures = ncol(Y), genoK = NA,
      test = c("Mantel", "partial Mantel (| treatment)"), term = "genotype distance vs methylation distance",
      R2 = c(mt$statistic, pm$statistic), F = NA, p = c(mt$signif, pm$signif))

    # Per-feature partitioning, observed and nulls
    obs <- partitionFeatures(Y, meta, g$scores, genoK)
    featureSummaries[[length(featureSummaries) + 1]] <- summarizeFeatures(obs) %>% mutate(sampleSet = s, features = fs, genoK = genoK, .before = 1)
    if (fs == "gene") write_csv(bind_cols(tibble(region = colnames(Y)), obs), file.path("tables", paste0("per-gene-partition-", s, ".csv")))
    nulls <- bind_rows(
      map_dfr(seq_len(nNull), ~ summarizeFeatures(partitionFeatures(Y, permuteLabels(meta, .x), g$scores, genoK)) %>% mutate(null = "treatment labels shuffled", seed = .x)),
      map_dfr(seq_len(nNull), function(i) { set.seed(10000 + i); summarizeFeatures(partitionFeatures(Y, meta, g$scores[sample(nrow(Y)), , drop = FALSE], genoK)) %>% mutate(null = "genotype shuffled", seed = i) }))
    nullSummaries[[length(nullSummaries) + 1]] <- nulls %>% mutate(sampleSet = s, features = fs, .before = 1)

    plotData <- obs %>% dplyr::select(uniqueTreatment, uniqueGenotype) %>% pivot_longer(everything(), names_to = "fraction", values_to = "R2")
    p <- ggplot(plotData, aes(R2, fill = fraction)) +
      geom_density(alpha = 0.5, color = NA) +
      geom_vline(data = nulls %>% group_by(null) %>% summarize(value = ifelse(null[1] == "genotype shuffled", median(median_uniqueGenotype), median(median_uniqueTreatment))),
                 aes(xintercept = value, linetype = null)) +
      labs(x = "Per-feature R2 (unique to each part)", y = "Density", fill = NULL, linetype = "Null median",
           subtitle = sprintf("%s, %s: %s features; genotype = %d PCs; treatment = ploidy * pH", s, fs, format(ncol(Y), big.mark = ","), genoK)) +
      theme_bw()
    ggsave(file.path("figures", sprintf("per-feature-R2-%s-%s.png", fs, s)), p, width = 7, height = 4.5, dpi = 150)
    rm(Y, methDist); invisible(gc())
  }
}

wholeResults <- bind_rows(wholeResults)
featureSummaries <- bind_rows(featureSummaries)
nullSummaries <- bind_rows(nullSummaries)
write_csv(bind_rows(genotypeInfo), file.path("tables", "genotype-PCA-variance.csv"))
write_csv(wholeResults, file.path("tables", "whole-methylome-tests.csv"))
write_csv(featureSummaries, file.path("tables", "per-feature-summary.csv"))
write_csv(nullSummaries, file.path("tables", "per-feature-nulls.csv"))
write_csv(nullSummaries %>% group_by(sampleSet, features, null) %>%
            summarize(across(where(is.numeric) & !seed, list(median = median, p95 = ~ quantile(.x, 0.95))), .groups = "drop"),
          file.path("tables", "per-feature-null-summary.csv"))

p <- wholeResults %>% filter(test == "varpart (adjusted R2)", term != "residual [d]") %>%
  mutate(term = factor(term, levels = c("treatment only [c]", "shared [b]", "genotype only [a]"))) %>%
  ggplot(aes(paste0(genoK, " PCs"), pmax(R2, 0), fill = term)) +
  geom_col() + facet_grid(features ~ sampleSet) +
  labs(x = "Genotype PCs in the model", y = "Adjusted R2 of methylation distance (negative shown as 0)", fill = NULL) +
  theme_bw()
ggsave(file.path("figures", "varpart-whole-methylome.png"), p, width = 7, height = 5.5, dpi = 150)

print(as.data.frame(wholeResults))
print(as.data.frame(featureSummaries))
print(as.data.frame(nullSummaries %>% group_by(sampleSet, features, null) %>%
                      summarize(across(c(median_uniqueTreatment, median_uniqueGenotype, fracGenotypeP0.01, fracTreatmentP0.01), median), .groups = "drop")))
sessionInfo()
