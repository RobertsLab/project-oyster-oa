#!/usr/bin/env Rscript
# Is methylation associated with the genotype of nearby SNPs (cis), beyond genome-wide relatedness?
# Plan: 04.8-cis-genotype-plan.md. Run from analyses/Haws_04.8-cis-genotype/ by 04.8-cis-genotype-submit.sh. One job;
# inputs come from 04.4 DSS and 04.6 methylKit:
#   ../Haws_04.4-DSS/prep-bssnper/CpG-SNP-genotypes.rds                                SNPs at CpGs (BS-SNPer PASS calls)
#   ../Haws_04.4-DSS/prep-cov5-hiperc99.9-cell5-meth10-noSNP-All/filtered-counts.rds  CpGs with no SNP in any oyster
#   ../Haws_04.6-methylKit-DMR/count-hiperc99.9-gene-cb3/region-counts.rds             reads per gene body
#
# SNP genotypes: bisulfite-safe SNPs only (no C>T at the C, no G>A at the G), at CpGs with >= 10 reads in every oyster of the
# set, so an oyster without a SNP call is reference (as for the 04.4 genotype PCA). Genotype = alt allele frequency.
# Kept if at least MIN_CARRIERS oysters carry the SNP and at least MIN_CARRIERS don't.
#
# Tests (each SNP-feature pair, residualized on covariates; t-test of the partial correlation):
#   CpG pairs: SNP vs each non-SNP CpG within MAX_DIST bp (methylation proportion, data in every oyster), binned by distance
#   gene pairs: SNP inside a gene body vs that gene's methylation (reads summed over the gene, >= 10 reads in every oyster)
# Covariates: "treatment" (effect-coded ploidy, pH, ploidy x pH, as in 04.4) and "treatment + genotype PCs" (genome-wide PCs
# 1..GENO_K from the same SNPs, recomputed per sample set, as in 04.7). The second asks whether local genotype matters beyond
# genome-wide relatedness.
# Null ("trans"): the same SNPs paired with features drawn at random from other chromosomes, NNULL times. Keeps each SNP's
# genotype and the genome-wide structure, and breaks only the physical link. Enrichment = observed / trans-null share of
# pairs at p < 0.001, and pi1 = 1 - pi0 (Storey, lambda = 0.5)
#
# Environment variables (defaults in brackets):
#   SAMPLE_SETS   [All DropGeno4]
#   GENO_K        [3]
#   MAX_DIST      [50000]
#   MIN_CARRIERS  [3]
#   NNULL         [20]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})
options(warn = 1)

sampleSets <- strsplit(Sys.getenv("SAMPLE_SETS", "All DropGeno4"), " ")[[1]]
genoK <- as.integer(Sys.getenv("GENO_K", "3"))
maxDist <- as.integer(Sys.getenv("MAX_DIST", "50000"))
minCarriers <- as.integer(Sys.getenv("MIN_CARRIERS", "3"))
nNull <- as.integer(Sys.getenv("NNULL", "20"))
dropSets <- list(All = integer(0), DropGeno4 = c(1, 2, 14, 17))
stopifnot(all(sampleSets %in% names(dropSets)))
distBreaks <- unique(c(c(0, 250, 500, 1000, 2000, 5000, 10000, 20000)[c(0, 250, 500, 1000, 2000, 5000, 10000, 20000) < maxDist], maxDist))

dssDir <- "../Haws_04.4-DSS"
mkDir <- "../Haws_04.6-methylKit-DMR"
dataDir <- "../../data/Haws"
mitoChr <- "NC_001276.1"
for (d in c("tables", "figures", "rds")) dir.create(d, showWarnings = FALSE)

sampleMetadata <- read.csv(file.path(dataDir, "sample_metadata.csv")) %>%
  dplyr::select(-c(1)) %>%
  arrange(sample_number) %>%
  mutate(sampleID = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                      "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                      "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                      "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6")) %>%
  dplyr::rename(pH = ph)
stopifnot(identical(sampleMetadata$sample_number, 1:24))
treatmentDesign <- function(meta) { #Effect coding, as in 04.4-DSS-slurm.R and 04.7
  p <- ifelse(meta$ploidy == "3N", 0.5, -0.5); h <- ifelse(meta$pH == "low", 0.5, -0.5)
  cbind(ploidy = p, pH = h, ploidyXpH = p * h)
}

# Inputs
geno <- readRDS(file.path(dssDir, "prep-bssnper", "CpG-SNP-genotypes.rds"))
cpg <- readRDS(file.path(dssDir, "prep-cov5-hiperc99.9-cell5-meth10-noSNP-All", "filtered-counts.rds"))
cpgBeta <- cpg$M / cpg$Cov
cpgBeta[cpg$Cov == 0] <- NA
colnames(cpgBeta) <- sampleMetadata$sampleID
cpgInfo <- data.table(chr = cpg$chr, pos = cpg$pos)
rm(cpg); invisible(gc())

geneCounts <- readRDS(file.path(mkDir, "count-hiperc99.9-gene-cb3", "region-counts.rds"))
stopifnot(identical(vapply(geneCounts, function(x) x@sample.id, ""), sampleMetadata$sampleID))
geneData <- lapply(geneCounts, function(x) { d <- methylKit::getData(x); data.table(region = paste(d$chr, d$start, d$end), cov = d$coverage, beta = d$numCs / d$coverage) })
geneKeys <- Reduce(union, lapply(geneData, `[[`, "region"))
geneBeta <- sapply(geneData, function(d) { b <- d$beta; b[d$cov < 10] <- NA; b[match(geneKeys, d$region)] })
colnames(geneBeta) <- sampleMetadata$sampleID
geneInfo <- as.data.table(str_split_fixed(geneKeys, " ", 3)) %>% setnames(c("chr", "start", "end")) %>%
  mutate(start = as.integer(start), end = as.integer(end))
rm(geneCounts, geneData); invisible(gc())

# Residualize the rows of X (features x oysters) on covariates C (oysters x p, with intercept)
residualize <- function(X, C) X - X %*% C %*% solve(crossprod(C), t(C))
# Partial-correlation test for pairs (rows i of A, rows j of B), both already residualized
pairTest <- function(A, B, i, j, df) {
  num <- rowSums(A[i, , drop = FALSE] * B[j, , drop = FALSE])
  r <- num / sqrt(rowSums(A[i, , drop = FALSE]^2) * rowSums(B[j, , drop = FALSE]^2))
  t <- r * sqrt(df / pmax(1 - r^2, 1e-12))
  list(r = r, p = 2 * pt(-abs(t), df))
}
pi1 <- function(p) 1 - min(1, mean(p > 0.5) / 0.5)
summarizeP <- function(p) c(pairs = length(p), fracP001 = mean(p < 0.001), fracP01 = mean(p < 0.01), pi1 = pi1(p))

cisResults <- list(); geneResults <- list(); snpCounts <- list()
for (s in sampleSets) {
  keep <- setdiff(1:24, dropSets[[s]])
  meta <- sampleMetadata[keep, ]
  ids <- meta$sampleID
  n <- length(ids)

  # SNPs: bisulfite-safe, callable in every oyster, enough carriers and non-carriers
  carriers <- rowSums(geno$snpNonBS[, ids])
  use <- rowSums(geno$Cov[, ids] >= 10) == n & carriers >= minCarriers & carriers <= n - minCarriers & geno$chr != mitoChr
  G <- geno$altFreqNonBS[use, ids]
  snpInfo <- data.table(chr = geno$chr[use], pos = geno$pos[use], carriers = carriers[use])
  pcs <- prcomp(t(G), center = TRUE, scale. = FALSE)$x[, 1:genoK, drop = FALSE]
  message(s, ": ", nrow(G), " SNPs")

  # Features with data in every oyster of the set
  cpgUse <- which(rowSums(is.na(cpgBeta[, ids])) == 0 & cpgInfo$chr != mitoChr)
  Y <- cpgBeta[cpgUse, ids]; yInfo <- cpgInfo[cpgUse]
  geneUse <- which(rowSums(is.na(geneBeta[, ids])) == 0 & geneInfo$chr != mitoChr)
  Z <- geneBeta[geneUse, ids]; zInfo <- geneInfo[geneUse]

  # Physical pairs. CpG pairs within maxDist (the noSNP CpG set already excludes every CpG with a SNP)
  yInfo[, idx := .I]; snpInfo[, idx := .I]
  cpgPairs <- rbindlist(lapply(split(snpInfo, snpInfo$chr), function(sn) {
    cp <- yInfo[chr == sn$chr[1]]
    if (nrow(cp) == 0) return(NULL)
    lo <- findInterval(sn$pos - maxDist - 1, cp$pos) + 1
    hi <- findInterval(sn$pos + maxDist, cp$pos)
    ok <- hi >= lo
    if (!any(ok)) return(NULL)
    data.table(snp = rep(sn$idx[ok], hi[ok] - lo[ok] + 1),
               feature = cp$idx[unlist(mapply(seq, lo[ok], hi[ok], SIMPLIFY = FALSE))])
  }))
  cpgPairs[, distance := abs(yInfo$pos[feature] - snpInfo$pos[snp])]
  cpgPairs <- cpgPairs[distance > 1] #Skip the G of the SNP's own CpG
  cpgPairs[, bin := cut(distance, distBreaks, include.lowest = TRUE, dig.lab = 6)]
  zInfo[, idx := .I]
  genePairs <- rbindlist(lapply(split(snpInfo, snpInfo$chr), function(sn) {
    gn <- zInfo[chr == sn$chr[1]]
    if (nrow(gn) == 0) return(NULL)
    hits <- lapply(sn$pos, function(p) gn$idx[gn$start <= p & gn$end >= p])
    data.table(snp = rep(sn$idx, lengths(hits)), feature = unlist(hits))
  }))
  snpCounts[[s]] <- tibble(sampleSet = s, SNPs = nrow(G), CpGs = nrow(Y), genes = nrow(Z),
                           SNPCpGPairs = nrow(cpgPairs), SNPsWithCpGPair = uniqueN(cpgPairs$snp),
                           SNPGenePairs = nrow(genePairs), SNPsInGenes = uniqueN(genePairs$snp))
  message(s, ": ", nrow(cpgPairs), " SNP-CpG pairs, ", nrow(genePairs), " SNP-gene pairs")

  for (covSet in c("treatment", "treatment + genotype PCs")) {
    C <- cbind(1, treatmentDesign(meta), if (covSet == "treatment") NULL else scale(pcs))
    df <- n - ncol(C) - 1
    Gr <- residualize(G, C); Yr <- residualize(Y, C); Zr <- residualize(Z, C)

    obs <- pairTest(Gr, Yr, cpgPairs$snp, cpgPairs$feature, df)
    cpgPairs[, `:=`(r = obs$r, p = obs$p)]
    if (covSet == "treatment + genotype PCs")
      saveRDS(cbind(cpgPairs, snpChr = snpInfo$chr[cpgPairs$snp], snpPos = snpInfo$pos[cpgPairs$snp], cpgPos = yInfo$pos[cpgPairs$feature]),
              file.path("rds", paste0("cis-CpG-pairs-", s, ".rds")))
    gobs <- pairTest(Gr, Zr, genePairs$snp, genePairs$feature, df)

    # Trans null: replace each pair's feature with a random feature on another chromosome
    yChr <- yInfo$chr; zChr <- zInfo$chr
    drawOther <- function(featureChr, snpChr) {
      out <- sample(length(featureChr), length(snpChr), replace = TRUE)
      bad <- featureChr[out] == snpChr
      while (any(bad)) { out[bad] <- sample(length(featureChr), sum(bad), replace = TRUE); bad <- featureChr[out] == snpChr }
      out
    }
    nullCpG <- list(); nullGene <- list()
    for (k in seq_len(nNull)) {
      set.seed(k)
      nf <- drawOther(yChr, snpInfo$chr[cpgPairs$snp])
      nullCpG[[k]] <- data.table(bin = cpgPairs$bin, p = pairTest(Gr, Yr, cpgPairs$snp, nf, df)$p, null = k)
      nf <- drawOther(zChr, snpInfo$chr[genePairs$snp])
      nullGene[[k]] <- data.table(p = pairTest(Gr, Zr, genePairs$snp, nf, df)$p, null = k)
    }
    nullCpG <- rbindlist(nullCpG); nullGene <- rbindlist(nullGene)

    binObs <- cpgPairs[, as.list(summarizeP(p)), by = bin]
    binNull <- nullCpG[, as.list(summarizeP(p)), by = .(bin, null)][
      , .(nullFracP001 = mean(fracP001), nullFracP001_95 = quantile(fracP001, 0.95), nullPi1 = mean(pi1), nullPi1_95 = quantile(pi1, 0.95)), by = bin]
    allObs <- as.list(summarizeP(cpgPairs$p))
    allNull <- nullCpG[, as.list(summarizeP(p)), by = null][, .(nullFracP001 = mean(fracP001), nullFracP001_95 = quantile(fracP001, 0.95), nullPi1 = mean(pi1), nullPi1_95 = quantile(pi1, 0.95))]
    cisResults[[length(cisResults) + 1]] <- bind_rows(
      as_tibble(merge(binObs, binNull, by = "bin")) %>% mutate(bin = as.character(bin)),
      bind_cols(tibble(bin = paste0("all (<= ", maxDist, ")")), as_tibble(allObs), as_tibble(allNull))) %>%
      mutate(sampleSet = s, covariates = covSet, df = df, enrichmentP001 = fracP001 / nullFracP001, .before = 1)

    gNull <- nullGene[, as.list(summarizeP(p)), by = null]
    geneResults[[length(geneResults) + 1]] <- bind_cols(
      tibble(sampleSet = s, covariates = covSet, df = df), as_tibble(as.list(summarizeP(gobs$p))),
      tibble(nullFracP001 = mean(gNull$fracP001), nullFracP001_95 = quantile(gNull$fracP001, 0.95, names = FALSE),
             nullPi1 = mean(gNull$pi1), nullPi1_95 = quantile(gNull$pi1, 0.95, names = FALSE))) %>%
      mutate(enrichmentP001 = fracP001 / nullFracP001)
    if (covSet == "treatment + genotype PCs") {
      write_csv(tibble(snpChr = snpInfo$chr[genePairs$snp], snpPos = snpInfo$pos[genePairs$snp],
                       gene = paste(zInfo$chr[genePairs$feature], zInfo$start[genePairs$feature], zInfo$end[genePairs$feature]),
                       carriers = snpInfo$carriers[genePairs$snp], r = gobs$r, p = gobs$p) %>% arrange(p),
                file.path("tables", paste0("cis-gene-pairs-", s, ".csv")))
      qq <- bind_rows(tibble(type = "cis (observed)", p = sort(cpgPairs$p)),
                      tibble(type = "trans null (pooled)", p = sort(sample(nullCpG$p, min(nrow(nullCpG), nrow(cpgPairs)))))) %>%
        group_by(type) %>% mutate(expected = -log10(ppoints(n())), observed = -log10(p)) %>% ungroup()
      pq <- ggplot(qq, aes(expected, observed, color = type)) + geom_point(size = 0.4) + geom_abline(linetype = 2) +
        labs(x = "Expected -log10 p", y = "Observed -log10 p", color = NULL,
             subtitle = sprintf("%s: SNP-CpG pairs within %s bp, adjusted for treatment + %d genotype PCs", s, format(maxDist, big.mark = ","), genoK)) +
        theme_bw()
      ggsave(file.path("figures", paste0("qq-cis-CpG-", s, ".png")), pq, width = 6.5, height = 5, dpi = 150)
    }
  }
}

cisResults <- bind_rows(cisResults); geneResults <- bind_rows(geneResults)
write_csv(bind_rows(snpCounts), file.path("tables", "SNP-and-pair-counts.csv"))
write_csv(cisResults, file.path("tables", "cis-CpG-by-distance.csv"))
write_csv(geneResults, file.path("tables", "cis-gene.csv"))

p <- cisResults %>% filter(!startsWith(bin, "all")) %>%
  mutate(bin = factor(bin, levels = unique(bin))) %>%
  ggplot(aes(bin, enrichmentP001, color = covariates, group = covariates)) +
  geom_hline(yintercept = 1, linetype = 2) + geom_line() + geom_point() +
  facet_wrap(~ sampleSet) +
  labs(x = "Distance from SNP (bp)", y = "Share of pairs at p < 0.001, cis / trans null", color = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path("figures", "cis-enrichment-by-distance.png"), p, width = 8, height = 4.5, dpi = 150)

print(as.data.frame(bind_rows(snpCounts)))
print(as.data.frame(cisResults))
print(as.data.frame(geneResults))
sessionInfo()
