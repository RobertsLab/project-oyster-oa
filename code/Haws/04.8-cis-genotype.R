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
# Nulls: the same SNPs paired with features drawn at random, NNULL times, either
#   "other chromosome": from other chromosomes. Keeps each SNP's genotype and the genome-wide structure, and breaks the
#   physical link, but also breaks any relatedness that runs along a chromosome (shared haplotype stretches in families)
#   "same chromosome, far": from the same chromosome, more than FAR_DIST bp away. Keeps chromosome-scale relatedness, so
#   enrichment over this null is the part that is local
# Enrichment = observed / null share of pairs at p < 0.001, and pi1 = 1 - pi0 (Storey, lambda = 0.5)
#
# LD decay: r2 between the genotypes of every pair of these SNPs on the same chromosome, binned by distance, against random
# pairs on different chromosomes. Raw genotypes and genotypes with the genome-wide PCs removed (as in the cis test). Shows
# whether shared haplotypes (linkage) extend far enough to explain enrichment that stays flat out to 50 kb
#
# Methylation correlation decay (treatment + genotype PCs removed): correlation between the methylation of pairs of CpGs by
# distance (20,000 random anchor CpGs, each with every CpG within 100 kb), against CpGs on the same chromosome > FAR_DIST
# away and on other chromosomes. Shows whether methylation varies in regional blocks.
# Conditioning test: for SNPs with CpGs within 1 kb, SNP-CpG pairs 1-50 kb apart are retested with the mean methylation of
# the CpGs within 1 kb of the SNP as an extra covariate, and compared with the same-chromosome-far null treated the same
# way. If the flat cis enrichment beyond 1 kb is a local effect spreading through a methylation domain, conditioning on
# the local methylation should remove most of it
#
# Environment variables (defaults in brackets):
#   SAMPLE_SETS   [All DropGeno4]
#   GENO_K        [3]
#   MAX_DIST      [50000]
#   MIN_CARRIERS  [3]
#   NNULL         [20]
#   FAR_DIST      [1000000]

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
farDist <- as.integer(Sys.getenv("FAR_DIST", "1000000"))
nullTypes <- c("other chromosome", "same chromosome, far")
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
summarizeP <- function(p) { p <- p[!is.na(p)]; c(pairs = length(p), fracP001 = mean(p < 0.001), fracP01 = mean(p < 0.01), pi1 = pi1(p)) } #NA if a residualized genotype or feature is constant

cisResults <- list(); geneResults <- list(); snpCounts <- list(); ldResults <- list(); corDecay <- list(); conditioned <- list()
corBreaks <- c(0, 250, 500, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5)
unitRows <- function(X) X / sqrt(rowSums(X^2)) #Residualized rows have mean 0, so dot products of unit rows are correlations
ldBreaks <- c(0, 1e3, 5e3, 1e4, 2e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, Inf)
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

  # LD decay
  for (adj in c("none", "genotype PCs")) {
    Gx <- if (adj == "none") G else residualize(G, cbind(1, scale(pcs)))
    Gs <- t(scale(t(Gx))) / sqrt(n - 1) #Rows with unit length, so a cross product is a correlation
    chrSNPs <- split(seq_len(nrow(Gs)), geno$chr[use])
    within <- rbindlist(lapply(chrSNPs[lengths(chrSNPs) >= 2], function(i) {
      R <- tcrossprod(Gs[i, , drop = FALSE])
      pos <- geno$pos[use][i]
      ut <- which(upper.tri(R), arr.ind = TRUE)
      data.table(distance = abs(pos[ut[, 1]] - pos[ut[, 2]]), r2 = R[ut]^2)
    }))
    within[, bin := cut(distance, ldBreaks, include.lowest = TRUE, dig.lab = 10)]
    set.seed(1)
    a <- sample(nrow(Gs), 4e5, replace = TRUE); b <- sample(nrow(Gs), 4e5, replace = TRUE)
    other <- geno$chr[use][a] != geno$chr[use][b]
    a <- head(a[other], 2e5); b <- head(b[other], 2e5)
    otherR2 <- rowSums(Gs[a, , drop = FALSE] * Gs[b, , drop = FALSE])^2
    ldResults[[length(ldResults) + 1]] <- bind_rows(
      as_tibble(within[, .(pairs = .N, meanR2 = mean(r2), medianR2 = median(r2), fracR2above0.5 = mean(r2 > 0.5)), by = bin][order(bin)]) %>%
        mutate(bin = as.character(bin)),
      tibble(bin = "other chromosome", pairs = length(otherR2), meanR2 = mean(otherR2), medianR2 = median(otherR2), fracR2above0.5 = mean(otherR2 > 0.5))) %>%
      mutate(sampleSet = s, genotypeAdjustment = adj, oysters = n, .before = 1)
    rm(within); invisible(gc())
  }

  # Features with data in every oyster of the set
  cpgUse <- which(rowSums(is.na(cpgBeta[, ids])) == 0 & cpgInfo$chr != mitoChr)
  cpgUse <- cpgUse[matrixStats::rowVars(cpgBeta[cpgUse, ids]) > 0] #Constant features have no correlation
  Y <- cpgBeta[cpgUse, ids]; yInfo <- cpgInfo[cpgUse]
  geneUse <- which(rowSums(is.na(geneBeta[, ids])) == 0 & geneInfo$chr != mitoChr)
  geneUse <- geneUse[matrixStats::rowVars(geneBeta[geneUse, ids]) > 0]
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

    # Nulls: replace each pair's feature with a random feature on another chromosome, or on the same chromosome far away
    drawOther <- function(featureChr, snpChr) {
      out <- sample(length(featureChr), length(snpChr), replace = TRUE)
      bad <- featureChr[out] == snpChr
      while (any(bad)) { out[bad] <- sample(length(featureChr), sum(bad), replace = TRUE); bad <- featureChr[out] == snpChr }
      out
    }
    drawSameFar <- function(featureChr, featurePos, snpChr, snpPos) {
      out <- integer(length(snpChr))
      byChr <- split(seq_along(featureChr), featureChr)
      for (ch in unique(snpChr)) {
        i <- which(snpChr == ch); cand <- byChr[[ch]]
        if (is.null(cand)) { out[i] <- NA; next }
        draw <- cand[sample.int(length(cand), length(i), replace = TRUE)]
        bad <- abs(featurePos[draw] - snpPos[i]) <= farDist
        tries <- 0
        while (any(bad) && tries < 100) {
          draw[bad] <- cand[sample.int(length(cand), sum(bad), replace = TRUE)]
          bad <- abs(featurePos[draw] - snpPos[i]) <= farDist; tries <- tries + 1
        }
        draw[bad] <- NA #No feature far enough away on this chromosome
        out[i] <- draw
      }
      out
    }
    geneMid <- (zInfo$start + zInfo$end) / 2
    pairNull <- function(A, B, snp, feature, df) { ok <- !is.na(feature); p <- rep(NA_real_, length(snp)); p[ok] <- pairTest(A, B, snp[ok], feature[ok], df)$p; p }
    nullCpG <- list(); nullGene <- list()
    for (nt in nullTypes) for (k in seq_len(nNull)) {
      set.seed(k + ifelse(nt == nullTypes[1], 0, 1000))
      nf <- if (nt == nullTypes[1]) drawOther(yInfo$chr, snpInfo$chr[cpgPairs$snp]) else
        drawSameFar(yInfo$chr, yInfo$pos, snpInfo$chr[cpgPairs$snp], snpInfo$pos[cpgPairs$snp])
      nullCpG[[length(nullCpG) + 1]] <- data.table(nullType = nt, bin = cpgPairs$bin, p = pairNull(Gr, Yr, cpgPairs$snp, nf, df), null = k)
      nf <- if (nt == nullTypes[1]) drawOther(zInfo$chr, snpInfo$chr[genePairs$snp]) else
        drawSameFar(zInfo$chr, geneMid, snpInfo$chr[genePairs$snp], snpInfo$pos[genePairs$snp])
      nullGene[[length(nullGene) + 1]] <- data.table(nullType = nt, p = pairNull(Gr, Zr, genePairs$snp, nf, df), null = k)
    }
    nullCpG <- rbindlist(nullCpG); nullGene <- rbindlist(nullGene)

    binObs <- cpgPairs[, as.list(summarizeP(p)), by = bin]
    nullStats <- function(d) d[, .(nullFracP001 = mean(fracP001), nullFracP001_95 = quantile(fracP001, 0.95, names = FALSE),
                                   nullPi1 = mean(pi1), nullPi1_95 = quantile(pi1, 0.95, names = FALSE)), by = nullType]
    binNull <- nullCpG[, as.list(summarizeP(p)), by = .(nullType, bin, null)][
      , .(nullFracP001 = mean(fracP001), nullFracP001_95 = quantile(fracP001, 0.95, names = FALSE),
          nullPi1 = mean(pi1), nullPi1_95 = quantile(pi1, 0.95, names = FALSE)), by = .(nullType, bin)]
    allObs <- as.list(summarizeP(cpgPairs$p))
    allNull <- nullStats(nullCpG[, as.list(summarizeP(p)), by = .(nullType, null)])
    cisResults[[length(cisResults) + 1]] <- bind_rows(
      as_tibble(merge(binObs, binNull, by = "bin")) %>% mutate(bin = as.character(bin)),
      as_tibble(allNull) %>% mutate(bin = paste0("all (<= ", maxDist, ")"), !!!allObs)) %>%
      mutate(sampleSet = s, covariates = covSet, df = df, enrichmentP001 = fracP001 / nullFracP001, .before = 1)

    gNull <- nullStats(nullGene[, as.list(summarizeP(p)), by = .(nullType, null)])
    geneResults[[length(geneResults) + 1]] <- as_tibble(gNull) %>%
      mutate(sampleSet = s, covariates = covSet, df = df, !!!as.list(summarizeP(gobs$p)), .before = 1) %>%
      mutate(enrichmentP001 = fracP001 / nullFracP001)
    if (covSet == "treatment + genotype PCs") {
      # Methylation correlation decay
      Yu <- unitRows(Yr)
      set.seed(2)
      anchors <- sort(sample(nrow(Yu), min(20000, nrow(Yu))))
      near <- rbindlist(lapply(split(anchors, yInfo$chr[anchors]), function(a) {
        cp <- yInfo[chr == yInfo$chr[a[1]]]
        lo <- findInterval(yInfo$pos[a] - 1e5 - 1, cp$pos) + 1
        hi <- findInterval(yInfo$pos[a] + 1e5, cp$pos)
        data.table(anchor = rep(a, hi - lo + 1), other = cp$idx[unlist(mapply(seq, lo, hi, SIMPLIFY = FALSE))])
      }))[anchor != other]
      near[, distance := abs(yInfo$pos[anchor] - yInfo$pos[other])]
      near[, r := rowSums(Yu[anchor, , drop = FALSE] * Yu[other, , drop = FALSE])]
      near[, bin := cut(distance, corBreaks, include.lowest = TRUE, dig.lab = 10)]
      set.seed(3)
      a <- sample(anchors, 2e5, replace = TRUE)
      far <- drawSameFar(yInfo$chr, yInfo$pos, yInfo$chr[a], yInfo$pos[a]); okFar <- !is.na(far)
      oth <- drawOther(yInfo$chr, yInfo$chr[a])
      rFar <- rowSums(Yu[a[okFar], , drop = FALSE] * Yu[far[okFar], , drop = FALSE])
      rOth <- rowSums(Yu[a, , drop = FALSE] * Yu[oth, , drop = FALSE])
      corDecay[[length(corDecay) + 1]] <- bind_rows(
        as_tibble(near[, .(pairs = .N, meanR = mean(r), meanR2 = mean(r^2)), by = bin][order(bin)]) %>% mutate(bin = as.character(bin)),
        tibble(bin = paste0("same chromosome, > ", format(farDist, big.mark = ",", scientific = FALSE)), pairs = length(rFar), meanR = mean(rFar), meanR2 = mean(rFar^2)),
        tibble(bin = "other chromosome", pairs = length(rOth), meanR = mean(rOth), meanR2 = mean(rOth^2))) %>%
        mutate(sampleSet = s, oysters = n, .before = 1)
      rm(near); invisible(gc())

      # Conditioning on methylation within 1 kb of the SNP
      local <- cpgPairs[distance <= 1000]
      Mloc <- rowsum(Yr[local$feature, , drop = FALSE], local$snp) / as.vector(table(local$snp)[as.character(sort(unique(local$snp)))])
      locSNP <- as.integer(rownames(Mloc))
      Mu <- unitRows(Mloc)
      Gu <- unitRows(Gr)
      testPairs <- cpgPairs[distance > 1000 & distance <= 5e4 & snp %in% locSNP]
      condTest <- function(snp, feature) {
        m <- match(snp, locSNP)
        rgy <- rowSums(Gu[snp, , drop = FALSE] * Yu[feature, , drop = FALSE])
        rgm <- rowSums(Gu[snp, , drop = FALSE] * Mu[m, , drop = FALSE])
        rym <- rowSums(Yu[feature, , drop = FALSE] * Mu[m, , drop = FALSE])
        rc <- (rgy - rgm * rym) / sqrt(pmax((1 - rgm^2) * (1 - rym^2), 1e-12))
        pt2 <- function(r, d) 2 * pt(-abs(r * sqrt(d / pmax(1 - r^2, 1e-12))), d)
        tibble(unconditioned = pt2(rgy, df), conditioned = pt2(rc, df - 1))
      }
      obsCond <- condTest(testPairs$snp, testPairs$feature)
      nullCond <- map_dfr(seq_len(nNull), function(k) {
        set.seed(5000 + k)
        nf <- drawSameFar(yInfo$chr, yInfo$pos, snpInfo$chr[testPairs$snp], snpInfo$pos[testPairs$snp]); ok <- !is.na(nf)
        condTest(testPairs$snp[ok], nf[ok]) %>% summarize(across(everything(), ~ mean(.x < 0.001))) %>% mutate(null = k)
      })
      conditioned[[length(conditioned) + 1]] <- tibble(
        sampleSet = s, SNPsWithLocalCpG = length(locSNP), pairs1to50kb = nrow(testPairs),
        model = c("unconditioned", "conditioned on methylation within 1 kb"),
        fracP001 = c(mean(obsCond$unconditioned < 0.001), mean(obsCond$conditioned < 0.001)),
        nullFracP001 = c(mean(nullCond$unconditioned), mean(nullCond$conditioned)),
        nullFracP001_95 = c(quantile(nullCond$unconditioned, 0.95, names = FALSE), quantile(nullCond$conditioned, 0.95, names = FALSE))) %>%
        mutate(enrichmentP001 = fracP001 / nullFracP001)
    }

    if (covSet == "treatment + genotype PCs") {
      write_csv(tibble(snpChr = snpInfo$chr[genePairs$snp], snpPos = snpInfo$pos[genePairs$snp],
                       gene = paste(zInfo$chr[genePairs$feature], zInfo$start[genePairs$feature], zInfo$end[genePairs$feature]),
                       carriers = snpInfo$carriers[genePairs$snp], r = gobs$r, p = gobs$p) %>% arrange(p),
                file.path("tables", paste0("cis-gene-pairs-", s, ".csv")))
      qq <- bind_rows(tibble(type = "cis (observed)", p = sort(cpgPairs$p)),
                      map_dfr(nullTypes, function(nt) { np <- nullCpG[nullType == nt & !is.na(p), p]
                        tibble(type = paste0("null: ", nt), p = sort(sample(np, min(length(np), nrow(cpgPairs))))) })) %>%
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
corDecay <- bind_rows(corDecay); conditioned <- bind_rows(conditioned)
write_csv(corDecay, file.path("tables", "methylation-correlation-decay.csv"))
write_csv(conditioned, file.path("tables", "cis-conditioned-on-local-methylation.csv"))
pc <- corDecay %>% filter(!str_detect(bin, "chromosome")) %>%
  mutate(bin = factor(bin, levels = unique(bin))) %>%
  ggplot(aes(bin, meanR, color = sampleSet, group = sampleSet)) +
  geom_line() + geom_point() +
  geom_hline(data = corDecay %>% filter(str_detect(bin, "same chromosome")), aes(yintercept = meanR, color = sampleSet), linetype = 2) +
  geom_hline(data = corDecay %>% filter(bin == "other chromosome"), aes(yintercept = meanR, color = sampleSet), linetype = 3) +
  labs(x = "Distance between CpGs (bp)", y = "Mean correlation of methylation", color = NULL,
       subtitle = "Treatment + genotype PCs removed\nDashed: same chromosome > 1 Mb; dotted: other chromosome") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path("figures", "methylation-correlation-decay.png"), pc, width = 7, height = 4.5, dpi = 150)
ldResults <- bind_rows(ldResults)
write_csv(ldResults, file.path("tables", "LD-decay.csv"))
pl <- ldResults %>% filter(bin != "other chromosome") %>%
  mutate(bin = factor(bin, levels = unique(bin))) %>%
  ggplot(aes(bin, meanR2, color = genotypeAdjustment, group = genotypeAdjustment)) +
  geom_line() + geom_point() +
  geom_hline(data = ldResults %>% filter(bin == "other chromosome"), aes(yintercept = meanR2, color = genotypeAdjustment), linetype = 2) +
  facet_wrap(~ sampleSet) +
  labs(x = "Distance between SNPs (bp)", y = "Mean r2 (dashed: other chromosome)", color = "Genotype PCs removed") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path("figures", "LD-decay.png"), pl, width = 8, height = 4.5, dpi = 150)

p <- cisResults %>% filter(!startsWith(bin, "all")) %>%
  mutate(bin = factor(bin, levels = unique(bin))) %>%
  ggplot(aes(bin, enrichmentP001, color = covariates, linetype = nullType, group = interaction(covariates, nullType))) +
  geom_hline(yintercept = 1, linetype = 3) + geom_line() + geom_point() +
  facet_wrap(~ sampleSet) +
  labs(x = "Distance from SNP (bp)", y = "p < 0.001 share, observed / null", color = NULL, linetype = "Null") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.box = "vertical")
ggsave(file.path("figures", "cis-enrichment-by-distance.png"), p, width = 8, height = 5, dpi = 150)

print(as.data.frame(bind_rows(snpCounts)))
print(as.data.frame(cisResults))
print(as.data.frame(geneResults))
print(as.data.frame(ldResults))
print(as.data.frame(corDecay))
print(as.data.frame(conditioned))
sessionInfo()
