#!/usr/bin/env Rscript
# methylKit differentially methylated regions (DMRs) for the Haws ploidy x pH experiment. Plan: 04.6-methylKit-DMR-plan.md.
# Stages follow 04.3-methylKit-slurm.R. 04.6-methylKit-DMR-slurm-submit.sh runs them from analyses/Haws_04.6-methylKit-DMR/:
#   Rscript 04.6-methylKit-DMR-slurm.R prep       #Read the 24 coverage files (>= 1x) and remove high-coverage CpGs, mito exempt
#   Rscript 04.6-methylKit-DMR-slurm.R regions    #Annotation region sets (genes, exons, introns, flanks, lncRNA, TEs, mito)
#   Rscript 04.6-methylKit-DMR-slurm.R segments   #Group-blind segments from methylation pooled over all 24 oysters (needs prep)
#   Rscript 04.6-methylKit-DMR-slurm.R count      #Reads per region for one REGIONS / COV_BASES setting (needs prep, and regions or segments)
#   Rscript 04.6-methylKit-DMR-slurm.R dmr        #One test per array task: ploidy or pH x min.per.group All/10/8, from SLURM_ARRAY_TASK_ID
#   Rscript 04.6-methylKit-DMR-slurm.R summary    #Combine counts from every run, and compare observed counts with permutations
#
# Environment variables (defaults in brackets):
#   HI_PERC         [99.9]      Per-sample upper coverage percentile, calculated on nuclear CpGs only. Mito CpGs are exempt,
#                               since mito coverage (~300x) is far above the nuclear cutoff (~80x). Use "none" to skip
#   REGIONS         [tile1000]  tile<W> (non-overlapping W bp windows), an annotation set (gene, exonUTR, intron, upstream,
#                               downstream, lncRNA, TE, TE-DNA, TE-RC, TE-LINE, TE-LTR, repeat-unknown, repeat-simple, mito), or segments
#   COV_BASES       [3]         Minimum CpGs with coverage in a region, per sample
#   LO_COUNT        [10]        Minimum reads per region (summed over its CpGs), per sample
#   SAMPLES         [All]       "All" (24 oysters) or "OutRM" (drop 2H-3 and 3H-2)
#   OVERDISPERSION  [MN]        "MN" (primary), "shrinkMN", or "none" (negative control only)
#   PERM            [0]         0 = observed labels. Any other value is the seed for a label permutation. Uses the same shuffle
#                               as 04.4-DSS-slurm.R, so the same seed gives the same permuted labels in both analyses
#
# DMR threshold (fixed in the plan before any results): methylKit q < 0.05 and |meth.diff| >= 10%.
# Counts at q < 0.01 and |meth.diff| >= 15/25% are written as sensitivity checks only.
# meth.diff is treatment minus control: positive = higher methylation in triploids (ploidy) or at low pH (pH)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(methylKit)
})
options(warn = 1) #Print warnings as they happen, so they show up in the Slurm logs

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
stopifnot(mode %in% c("prep", "regions", "segments", "count", "dmr", "summary"))

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) #Use all cores allocated by Slurm
hiPercSetting <- Sys.getenv("HI_PERC", "99.9")
regionSet <- Sys.getenv("REGIONS", "tile1000")
covBases <- as.integer(Sys.getenv("COV_BASES", "3"))
loCount <- as.integer(Sys.getenv("LO_COUNT", "10"))
samples <- Sys.getenv("SAMPLES", "All")
overdispersion <- Sys.getenv("OVERDISPERSION", "MN")
perm <- as.integer(Sys.getenv("PERM", "0"))

annotationSets <- c("gene", "exonUTR", "intron", "upstream", "downstream", "lncRNA",
                    "TE", "TE-DNA", "TE-RC", "TE-LINE", "TE-LTR", "repeat-unknown", "repeat-simple", "mito")
stopifnot(hiPercSetting == "none" || !is.na(as.numeric(hiPercSetting)),
          grepl("^tile[0-9]+$", regionSet) || regionSet %in% c(annotationSets, "segments"),
          covBases >= 1, loCount >= 1,
          samples %in% c("All", "OutRM"),
          overdispersion %in% c("none", "MN", "shrinkMN"))

dataDir <- "../../data/Haws"
featureDir <- "../../genome-feature-files"
rmOutFile <- "downloads/GCF_902806645.1_cgigas_uk_roslin_v1_rm.out.gz" #Downloaded by the submit script
mitoChr <- "NC_001276.1"
outliers <- c(3, 14) #Samples 2H-3 and 3H-2
presenceSettings <- c("All", "10", "8")
qCuts <- c(0.05, 0.01)
diffCuts <- c(10, 15, 25)
primaryQ <- 0.05
primaryDiff <- 10

prepDir <- paste0("prep-hiperc", hiPercSetting)
segmentDir <- paste0("segments-hiperc", hiPercSetting)
countDir <- paste0("count-hiperc", hiPercSetting, "-", regionSet, "-cb", covBases)
runTag <- paste0(regionSet, "-cb", covBases, "-lo", loCount, "-hiperc", hiPercSetting, "-", samples, "-", overdispersion)
runDir <- paste0("DMR-", runTag, ifelse(perm == 0, "", paste0("-perm", perm)))

tasks <- expand.grid(min.per.group = presenceSettings, test = c("ploidy", "pH"),
                     stringsAsFactors = FALSE)[, c("test", "min.per.group")] #Tasks 1-3 are ploidy, 4-6 are pH

# Sample metadata, as in 04.2

sampleMetadata <- read.csv(file.path(dataDir, "sample_metadata.csv")) %>%
  dplyr::select(-c(1)) %>%
  arrange(sample_number) %>%
  mutate("sampleID" = c("2H-1", "2H-2", "2H-3", "2H-4", "2H-5", "2H-6",
                        "2L-1", "2L-2", "2L-3", "2L-4", "2L-5", "2L-6",
                        "3H-1", "3H-2", "3H-3", "3H-4", "3H-5", "3H-6",
                        "3L-1", "3L-2", "3L-3", "3L-4", "3L-5", "3L-6")) %>%
  dplyr::rename(pH = ph)
stopifnot(identical(sampleMetadata$sample_number, 1:24)) #File zr3644_N must match row N

keep <- if (samples == "OutRM") setdiff(1:24, outliers) else 1:24

# Label permutation, copied from makeDesign() in 04.4-DSS-slurm.R so a seed gives the same labels in both analyses.
# Shuffles ploidy within pH, then pH within the new ploidy labels. Keeps every group the same size
permuteLabels <- function(metadata, perm) {
  if (perm == 0) return(metadata)
  set.seed(perm)
  metadata %>%
    group_by(pH) %>% mutate(ploidy = sample(ploidy)) %>% ungroup() %>%
    group_by(ploidy) %>% mutate(pH = sample(pH)) %>% ungroup()
}

# Read a GFF from genome-feature-files as GRanges (1-based, as in the GFF), with the ID attribute as the name
readGFF <- function(file) {
  gff <- fread(file.path(featureDir, file), header = FALSE, sep = "\t", quote = "",
               select = c(1, 4, 5, 7, 9), col.names = c("chr", "start", "end", "strand", "attributes"))
  GRanges(gff$chr, IRanges(gff$start, gff$end), strand = gff$strand,
          name = str_match(gff$attributes, "ID=([^;]+)")[, 2])
}

# Merge overlapping intervals so a CpG is only counted in one region of a set. Strand is ignored, since the
# coverage files are already merged across strands
mergeRegions <- function(gr) reduce(gr, ignore.strand = TRUE)

# Stage 1: read coverage files and remove high-coverage CpGs (likely PCR duplicates or repeats) in each sample.
# Low coverage is filtered later, on summed region counts, so every CpG with at least 1 read is kept here

if (mode == "prep") {
  dir.create(prepDir, showWarnings = FALSE)
  covFiles <- file.path(dataDir, paste0("zr3644_", 1:24, "_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov"))
  stopifnot(all(file.exists(covFiles)))

  raw <- methylKit::methRead(as.list(covFiles),
                             sample.id = as.list(sampleMetadata$sampleID),
                             assembly = "cgigas_uk_roslin_v1",
                             treatment = ifelse(sampleMetadata$ploidy == "3N", 1, 0),
                             pipeline = "bismarkCoverage",
                             header = FALSE, #Bismark coverage files have no header. The default (TRUE) drops the first CpG, which is the first mito CpG
                             mincov = 1)

  # Same rule as filterByCoverage(hi.perc), which drops CpGs above the percentile, except that the percentile is
  # calculated on nuclear CpGs only and mito CpGs are never removed
  filterStats <- list()
  for (i in seq_along(raw)) {
    covData <- methylKit::getData(raw[[i]])
    isMito <- covData$chr == mitoChr
    cutoff <- if (hiPercSetting == "none") Inf else quantile(covData$coverage[!isMito], as.numeric(hiPercSetting) / 100, names = FALSE)
    keepRows <- which(covData$coverage <= cutoff | isMito)
    filterStats[[i]] <- tibble(sampleID = raw[[i]]@sample.id,
                               hiCutoff = cutoff,
                               CpGs = nrow(covData),
                               CpGsRemoved = nrow(covData) - length(keepRows),
                               mitoCpGs = sum(isMito),
                               mitoMeanCoverage = mean(covData$coverage[isMito]),
                               nuclearMeanCoverage = mean(covData$coverage[!isMito]))
    raw[[i]] <- methylKit::select(raw[[i]], keepRows)
    message(raw[[i]]@sample.id, ": cutoff ", cutoff, "x, removed ", filterStats[[i]]$CpGsRemoved, " CpGs, kept ", sum(isMito), " mito CpGs")
  }
  rm(covData)

  write_tsv(bind_rows(filterStats), file.path(prepDir, "high-coverage-filter.tsv"))
  saveRDS(raw, file.path(prepDir, "filtered-raw.rds"), compress = FALSE)
  file.create(file.path(prepDir, "prep-complete")) #Lets the submit script skip prep on later runs
}

# Stage 2: annotation region sets. Genes keep one region per gene, so results can be used for gene-level enrichment.
# All other sets are merged (see mergeRegions)

if (mode == "regions") {
  dir.create("regions", showWarnings = FALSE)
  stopifnot(file.exists(rmOutFile))

  genes <- readGFF("cgigas_uk_roslin_v1_gene.gff")
  exonUTR <- mergeRegions(readGFF("cgigas_uk_roslin_v1_exonUTR.gff"))

  # TEs from the NCBI RepeatMasker output. The lab bed (cgigas_uk_roslin_v1_rm.te.bed) is this file without the
  # class column, so it also has simple repeats, and its starts are not converted to 0-based
  rmLines <- readLines(gzfile(rmOutFile))[-(1:3)] #Three header lines
  rmLines <- rmLines[nzchar(trimws(rmLines))]
  rmFields <- strsplit(trimws(rmLines), "[[:space:]]+")
  rm(rmLines)
  repeats <- tibble(chr = map_chr(rmFields, 5),
                    start = as.integer(map_chr(rmFields, 6)), #1-based, inclusive
                    end = as.integer(map_chr(rmFields, 7)),
                    strand = ifelse(map_chr(rmFields, 9) == "C", "-", "+"),
                    repeatName = map_chr(rmFields, 10),
                    classFamily = map_chr(rmFields, 11)) %>%
    mutate(class = str_remove(classFamily, "/.*"))
  rm(rmFields)
  teClasses <- c("DNA", "RC", "LINE", "LTR", "SINE")
  repeatGR <- function(classes) {
    r <- filter(repeats, class %in% classes)
    GRanges(r$chr, IRanges(r$start, r$end))
  }
  te <- filter(repeats, class %in% teClasses)
  write_tsv(transmute(te, chr, start = start - 1L, end, name = repeatName, classFamily, strand) %>% arrange(chr, start),
            "regions/cgigas_uk_roslin_v1_TE.bed", col_names = FALSE) #0-based BED, TEs only, for 07 and other downstream steps

  # Mito genome as one region. It isn't in the genome-feature-files sequence lengths, so take its length from the
  # last CpG with coverage in sample 1
  mitoEnd <- fread(cmd = paste("grep", mitoChr, shQuote(file.path(dataDir, "zr3644_1_R1_val_1_val_1_val_1_bismark_bt2_pe..CpG_report.merged_CpG_evidence.cov"))),
                   select = 2)[[1]] %>% max()

  regionList <- list(gene = genes,
                     exonUTR = exonUTR,
                     intron = GenomicRanges::setdiff(mergeRegions(genes), exonUTR),
                     upstream = mergeRegions(readGFF("cgigas_uk_roslin_v1_upstream.gff")),
                     downstream = mergeRegions(readGFF("cgigas_uk_roslin_v1_downstream.gff")),
                     lncRNA = mergeRegions(readGFF("cgigas_uk_roslin_v1_lncRNA.gff")),
                     TE = mergeRegions(repeatGR(teClasses)),
                     "TE-DNA" = mergeRegions(repeatGR("DNA")),
                     "TE-RC" = mergeRegions(repeatGR("RC")),
                     "TE-LINE" = mergeRegions(repeatGR("LINE")),
                     "TE-LTR" = mergeRegions(repeatGR("LTR")),
                     "repeat-unknown" = mergeRegions(repeatGR("Unknown")),
                     "repeat-simple" = mergeRegions(repeatGR(c("Simple_repeat", "Low_complexity", "Satellite"))),
                     mito = GRanges(mitoChr, IRanges(1, mitoEnd)))
  stopifnot(setequal(names(regionList), annotationSets))

  regionSummary <- imap_dfr(regionList, ~ tibble(regions = .y, intervals = length(.x), Mb = sum(as.numeric(width(.x))) / 1e6,
                                                 medianWidth = median(width(.x))))
  print(as.data.frame(regionSummary))
  write_tsv(regionSummary, "regions/region-set-summary.tsv")
  write_tsv(count(repeats, class, name = "intervals"), "regions/RepeatMasker-class-counts.tsv")
  saveRDS(regionList, "regions/annotation-regions.rds")
  file.create("regions/regions-complete")
}

# Stage 3: group-blind segments (plan C). Methylation is summed over all 24 oysters at each CpG, ignoring treatment,
# and methSeg splits the genome into segments of similar methylation. Segmenting on treatment differences instead
# would make the later tests on those segments look more significant than they are

if (mode == "segments") {
  dir.create(segmentDir, showWarnings = FALSE)
  raw <- readRDS(file.path(prepDir, "filtered-raw.rds"))
  pooled <- rbindlist(lapply(raw, function(x) as.data.table(methylKit::getData(x))[, .(chr, start, end, numCs, numTs)]))
  rm(raw)
  pooled <- pooled[, .(numCs = sum(numCs), numTs = sum(numTs)), by = .(chr, start, end)]
  pooled[, coverage := numCs + numTs]
  pooled <- pooled[coverage >= 10][order(chr, start)] #At least 10 reads over all oysters
  pooledRaw <- new("methylRaw",
                   data.frame(chr = pooled$chr, start = pooled$start, end = pooled$end, strand = "*",
                              coverage = pooled$coverage, numCs = pooled$numCs, numTs = pooled$numTs),
                   sample.id = "pooled", assembly = "cgigas_uk_roslin_v1", context = "CpG", resolution = "base")
  message(nrow(pooled), " CpGs with >= 10 pooled reads")
  rm(pooled)

  segments <- methylKit::methSeg(pooledRaw, diagnostic.plot = FALSE, initialize.on.subset = 0.1, minSeg = 10) #minSeg: at least 10 CpGs per segment (passed to fastseg)
  segments <- segments[width(segments) > 1]
  write_tsv(as.data.frame(segments) %>% dplyr::select(chr = seqnames, start, end, num.mark, seg.mean, seg.group),
            file.path(segmentDir, "segments.tsv"))
  saveRDS(segments, file.path(segmentDir, "segments.rds"))
  file.create(file.path(segmentDir, "segments-complete"))
}

# Stage 4: reads per region in each sample. Tiles are counted from the CpG-level data, so they aren't affected by
# the per-CpG normalization in 04.3. Coverage filtering and normalization are done per test, in the dmr stage

if (mode == "count") {
  dir.create(countDir, showWarnings = FALSE)
  raw <- readRDS(file.path(prepDir, "filtered-raw.rds"))

  if (grepl("^tile", regionSet)) {
    windowSize <- as.integer(sub("tile", "", regionSet))
    regionCounts <- methylKit::tileMethylCounts(raw, win.size = windowSize, step.size = windowSize, #Non-overlapping, so each CpG is tested once
                                                cov.bases = covBases, mc.cores = nCores)
  } else {
    regions <- if (regionSet == "segments") readRDS(file.path(segmentDir, "segments.rds")) else readRDS("regions/annotation-regions.rds")[[regionSet]]
    regionCounts <- methylKit::regionCounts(raw, regions = regions, cov.bases = covBases, strand.aware = FALSE)
  }

  write_tsv(tibble(sampleID = map_chr(regionCounts, ~ .x@sample.id),
                   regions = map_int(regionCounts, ~ nrow(methylKit::getData(.x))),
                   medianReads = map_dbl(regionCounts, ~ median(methylKit::getData(.x)$coverage))),
            file.path(countDir, "regions-per-sample.tsv"))
  saveRDS(regionCounts, file.path(countDir, "region-counts.rds"), compress = FALSE)
  file.create(file.path(countDir, "count-complete"))
}

# Stage 5: one test per array task

if (mode == "dmr") {
  taskID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  stopifnot(taskID %in% seq_len(nrow(tasks)))
  test <- tasks$test[taskID]
  presence <- tasks$min.per.group[taskID]
  message("Task ", taskID, ": ", test, ", min.per.group = ", presence, ", ", runDir, ", cores = ", nCores)
  for (d in file.path(runDir, c("counts", "rds", "pvalue-histograms"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  metadata <- permuteLabels(sampleMetadata[keep, ], perm)
  treatment <- if (test == "ploidy") as.integer(metadata$ploidy == "3N") else as.integer(metadata$pH == "low")
  covariates <- if (test == "ploidy") data.frame(pH = metadata$pH) else data.frame(ploidy = metadata$ploidy) #Other factor is the covariate

  regionCounts <- readRDS(file.path(countDir, "region-counts.rds")) %>%
    methylKit::reorganize(sample.ids = metadata$sampleID, treatment = treatment) %>%
    methylKit::filterByCoverage(lo.count = loCount) %>%
    methylKit::normalizeCoverage() #Median normalization over the samples in this test

  methylationInformation <- if (presence == "All") {
    methylKit::unite(regionCounts, destrand = FALSE, mc.cores = 2) #Regions with data in all samples
  } else {
    methylKit::unite(regionCounts, destrand = FALSE, mc.cores = 2, min.per.group = as.integer(presence))
  }
  rm(regionCounts)

  differentialMethylationStats <- methylKit::calculateDiffMeth(methylationInformation,
                                                               covariates = covariates,
                                                               overdispersion = overdispersion,
                                                               mc.cores = nCores) #Default SLIM q-values, as in 04.3

  results <- methylKit::getData(differentialMethylationStats) %>%
    mutate(qvalueBH = p.adjust(pvalue, "BH")) #BH q-values for comparison with DSS
  if (regionSet == "gene") { #Add gene IDs. Genes with the same coordinates get all their IDs
    geneIDs <- as.data.frame(readRDS("regions/annotation-regions.rds")$gene) %>%
      group_by(chr = as.character(seqnames), start, end) %>%
      summarize(regionID = paste(name, collapse = ","), .groups = "drop")
    results <- left_join(results, geneIDs, by = c("chr", "start", "end"))
  }

  isMito <- results$chr == mitoChr
  counts <- expand_grid(qvalue = qCuts, difference = diffCuts) %>%
    mutate(DMR = map2_int(qvalue, difference, ~ sum(results$qvalue < .x & abs(results$meth.diff) >= .y)),
           hypermethylated = map2_int(qvalue, difference, ~ sum(results$qvalue < .x & results$meth.diff >= .y)),
           hypomethylated = map2_int(qvalue, difference, ~ sum(results$qvalue < .x & results$meth.diff <= -.y)),
           DMR.mito = map2_int(qvalue, difference, ~ sum(results$qvalue < .x & abs(results$meth.diff) >= .y & isMito)),
           DMR.BH = map2_int(qvalue, difference, ~ sum(results$qvalueBH < .x & abs(results$meth.diff) >= .y)),
           primary = qvalue == primaryQ & difference == primaryDiff)
  counts <- tibble(method = "methylKit", regions = regionSet, cov.bases = covBases, lo.count = loCount,
                   hi.perc = hiPercSetting, samples = samples, overdispersion = overdispersion, perm = perm,
                   test = test, min.per.group = presence,
                   regions.tested = nrow(results), regions.tested.mito = sum(isMito)) %>%
    bind_cols(counts)
  print(as.data.frame(counts))

  taskTag <- paste0(test, "-", presence)
  write_tsv(counts, file.path(runDir, "counts", paste0(taskTag, ".tsv")))
  pBins <- cut(results$pvalue, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)
  write_tsv(as_tibble(table(bin = pBins)), file.path(runDir, "pvalue-histograms", paste0(taskTag, ".tsv"))) #Calibration check (plan E1)

  if (perm == 0) { #Only keep full results and BED files for the observed labels
    saveRDS(results, file.path(runDir, "rds", paste0("diffMeth-", taskTag, ".rds")), compress = FALSE)
    dmr <- results %>%
      filter(qvalue < primaryQ, abs(meth.diff) >= primaryDiff) %>%
      transmute("#chr" = chr, start = start - 1L, end, across(any_of("regionID")), meth.diff, pvalue, qvalue, qvalueBH) %>% #0-based BED. Header starts with # so bedtools skips it
      arrange(`#chr`, start)
    write_tsv(dmr, file.path(runDir, paste0("DMR-", taskTag, "-q", primaryQ, "-diff", primaryDiff, ".bed")))
  }
}

# Stage 6: combine counts from every run in this directory (observed and permuted), and compare each observed
# count with its permutation distribution (plan E2)

if (mode == "summary") {
  countFiles <- list.files(".", pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE) %>%
    str_subset("^\\./DMR-[^/]+/counts/")
  stopifnot(length(countFiles) > 0)
  counts <- map_dfr(countFiles, read_tsv, col_types = cols(.default = "c", perm = "i", cov.bases = "i", lo.count = "i",
                                                           regions.tested = "i", regions.tested.mito = "i",
                                                           qvalue = "d", difference = "d", DMR = "i", hypermethylated = "i",
                                                           hypomethylated = "i", DMR.mito = "i", DMR.BH = "i", primary = "l")) %>%
    arrange(regions, cov.bases, lo.count, hi.perc, samples, overdispersion, test, min.per.group, qvalue, difference, perm)
  write_csv(counts, "DMR-counts-all.csv")

  settingColumns <- c("method", "regions", "cov.bases", "lo.count", "hi.perc", "samples", "overdispersion",
                      "test", "min.per.group", "qvalue", "difference", "primary")
  observed <- counts %>% filter(perm == 0) %>% dplyr::select(all_of(settingColumns), regions.tested, DMR)
  permuted <- counts %>%
    filter(perm != 0) %>%
    group_by(across(all_of(settingColumns))) %>%
    summarize(permutations = n(), perm.median = median(DMR), perm.95 = quantile(DMR, 0.95, names = FALSE),
              perm.DMR = list(DMR), .groups = "drop")
  permSummary <- left_join(observed, permuted, by = settingColumns) %>%
    mutate(empirical.p = map2_dbl(perm.DMR, DMR, ~ if (is.null(.x)) NA_real_ else (1 + sum(.x >= .y)) / (1 + length(.x))),
           above.perm.95 = DMR > perm.95) %>%
    dplyr::select(-perm.DMR)
  write_csv(permSummary, "DMR-permutation-summary.csv")
  print(as.data.frame(filter(permSummary, primary)))
}

sessionInfo()
