#!/usr/bin/env Rscript
# methylKit differentially methylated regions (DMRs) for the Haws ploidy x pH experiment. Plan: 04.6-methylKit-DMR-plan.md.
# Stages follow 04.3-methylKit-slurm.R. 04.6-methylKit-DMR-slurm-submit.sh runs them from analyses/Haws_04.6-methylKit-DMR/:
#   Rscript 04.6-methylKit-DMR-slurm.R prep       #Read the 24 coverage files (>= 1x) and remove high-coverage CpGs, mito exempt
#   Rscript 04.6-methylKit-DMR-slurm.R regions    #Annotation region sets (genes, exons, introns, flanks, lncRNA, TEs, mito)
#   Rscript 04.6-methylKit-DMR-slurm.R segments   #Group-blind segments from methylation pooled over all 24 oysters (needs prep)
#   Rscript 04.6-methylKit-DMR-slurm.R count      #Reads per region for one REGIONS / COV_BASES setting (needs prep, and regions or segments)
#   Rscript 04.6-methylKit-DMR-slurm.R dmr        #One test per array task: ploidy or pH x min.per.group All/10/8, from SLURM_ARRAY_TASK_ID
#   Rscript 04.6-methylKit-DMR-slurm.R summary    #Combine counts from every run, and compare observed counts with permutations
#   Rscript 04.6-methylKit-DMR-slurm.R enrich     #GO gene-set enrichment (GSEA) of gene-level results, tested against the label permutations
#   Rscript 04.6-methylKit-DMR-slurm.R simsummary #Region-level power simulation: detection of spiked-in effects (needs SIM_DELTA > 0 runs and permutations)
#
# Environment variables (defaults in brackets):
#   HI_PERC         [99.9]      Per-sample upper coverage percentile, calculated on nuclear CpGs only. Mito CpGs are exempt,
#                               since mito coverage (~300x) is far above the nuclear cutoff (~80x). Use "none" to skip
#   REGIONS         [tile1000]  tile<W> (non-overlapping W bp windows), an annotation set (gene, exonUTR, intron, upstream,
#                               downstream, lncRNA, TE, TE-DNA, TE-RC, TE-LINE, TE-LTR, repeat-unknown, repeat-simple, mito), or segments
#   COV_BASES       [3]         Minimum CpGs with coverage in a region, per sample
#   LO_COUNT        [10]        Minimum reads per region (summed over its CpGs), per sample
#   SAMPLES         [All]       "All" (24 oysters), or a sensitivity set (same sets as 04.4 D4, plus DropGeno4):
#                               "Drop3H2" (drop 3H-2), "DropPC2" (drop 2H-1, 2H-2, 3H-5), "DropGeno4" (drop all four genetic
#                               outliers from the 04.4 SNP check: 2H-1, 2H-2, 3H-2, 3H-5), or "OutRM" (drop 2H-3 and 3H-2, as in 04.2/04.3)
#   OVERDISPERSION  [MN]        "MN" (primary), "shrinkMN", or "none" (negative control only)
#   PERM            [0]         0 = observed labels. Any other value is the seed for a label permutation. Uses the same shuffle
#                               as 04.4-DSS-slurm.R, so the same seed gives the same permuted labels in both analyses
#   PERM_SCHEME     [free]      "free" = shuffle all oysters (above). "fixGeno4" = the four genetic outliers keep their real
#                               labels and only the other oysters are shuffled. All four outliers are high-pH oysters, so
#                               this keeps the real link between genotype and pH in every permutation
#   SIM_DELTA       [0]         Region-level power simulation (plan: 04.6-methylKit-DMR-power-simulation.md). > 0 = before the
#                               coverage filter, spike a difference of SIM_DELTA percentage points into SIM_N regions for ploidy and
#                               another SIM_N for pH, using the real labels and the real region counts (binomial thinning of the
#                               summed reads, coverage unchanged; same method as 04.4-DSS-slurm.R). Up to 30
#   SIM_N           [1000]      Number of spiked regions per term
#   SIM_REP         [1]         Replicate; the seed for which regions are spiked and their direction (same regions for every SIM_DELTA)
#
# DMR threshold (fixed in the plan before any results): methylKit q < 0.05 and |meth.diff| >= 10%.
# Counts at q < 0.01 and |meth.diff| >= 15/25% are written as sensitivity checks only. difference = 0 is a q-only tier (any
# effect size), added after the first results to check with permutations whether regions with small, precisely
# measured differences pass q < 0.05 more often than by chance. It is exploratory and does not change the DMR definition.
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
stopifnot(mode %in% c("prep", "regions", "segments", "count", "dmr", "summary", "enrich", "simsummary"))

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) #Use all cores allocated by Slurm
hiPercSetting <- Sys.getenv("HI_PERC", "99.9")
regionSet <- Sys.getenv("REGIONS", "tile1000")
covBases <- as.integer(Sys.getenv("COV_BASES", "3"))
loCount <- as.integer(Sys.getenv("LO_COUNT", "10"))
samples <- Sys.getenv("SAMPLES", "All")
overdispersion <- Sys.getenv("OVERDISPERSION", "MN")
perm <- as.integer(Sys.getenv("PERM", "0"))
permScheme <- Sys.getenv("PERM_SCHEME", "free")
simDelta <- as.integer(Sys.getenv("SIM_DELTA", "0"))
simN <- as.integer(Sys.getenv("SIM_N", "1000"))
simRep <- as.integer(Sys.getenv("SIM_REP", "1"))

annotationSets <- c("gene", "exonUTR", "intron", "upstream", "downstream", "lncRNA",
                    "TE", "TE-DNA", "TE-RC", "TE-LINE", "TE-LTR", "repeat-unknown", "repeat-simple", "mito")
stopifnot(hiPercSetting == "none" || !is.na(as.numeric(hiPercSetting)),
          grepl("^tile[0-9]+$", regionSet) || regionSet %in% c(annotationSets, "segments"),
          covBases >= 1, loCount >= 1,
          samples %in% c("All", "Drop3H2", "DropPC2", "DropGeno4", "OutRM"),
          permScheme %in% c("free", "fixGeno4"),
          permScheme == "free" || samples == "All", #fixGeno4 is for the all-24 analysis, where all four outliers are present
          overdispersion %in% c("none", "MN", "shrinkMN"),
          simDelta >= 0, simDelta <= 30, simN >= 1, simRep >= 1,
          simDelta == 0 || perm == 0) #Spike-ins use the real labels

dataDir <- "../../data/Haws"
featureDir <- "../../genome-feature-files"
rmOutFile <- "downloads/GCF_902806645.1_cgigas_uk_roslin_v1_rm.out.gz" #Downloaded by the submit script
gafFile <- "downloads/GCF_902806645.1_cgigas_uk_roslin_v1_gene_ontology.gaf.gz" #NCBI GO annotation, downloaded by the submit script
mitoChr <- "NC_001276.1"
dropSets <- list(All = integer(0),
                 Drop3H2 = 14,               #3H-2
                 DropPC2 = c(1, 2, 17),      #2H-1, 2H-2, 3H-5
                 DropGeno4 = c(1, 2, 14, 17), #2H-1, 2H-2, 3H-2, 3H-5: genetic outliers in the 04.4 SNP check
                 OutRM = c(3, 14))           #2H-3, 3H-2: the pair removed in 04.2/04.3
geneticOutliers <- c("2H-1", "2H-2", "3H-2", "3H-5")
presenceSettings <- c("All", "10", "8")
qCuts <- c(0.05, 0.01)
diffCuts <- c(0, 10, 15, 25) #0 = q-only tier
primaryQ <- 0.05
primaryDiff <- 10

prepDir <- paste0("prep-hiperc", hiPercSetting)
segmentDir <- paste0("segments-hiperc", hiPercSetting)
countDir <- paste0("count-hiperc", hiPercSetting, "-", regionSet, "-cb", covBases)
runTag <- paste0(regionSet, "-cb", covBases, "-lo", loCount, "-hiperc", hiPercSetting, "-", samples, "-", overdispersion)
permTag <- ifelse(permScheme == "free", "", paste0("-", permScheme)) #Only used for permuted runs
runDir <- paste0("DMR-", runTag, ifelse(perm == 0, "", paste0(permTag, "-perm", perm)),
                 ifelse(simDelta == 0, "", paste0("-simD", simDelta, "-n", simN, "-rep", simRep)))

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

keep <- setdiff(1:24, dropSets[[samples]])

# Label permutation, copied from makeDesign() in 04.4-DSS-slurm.R so a seed gives the same labels in both analyses.
# Shuffles ploidy within pH, then pH within the new ploidy labels. Keeps every group the same size
permuteLabels <- function(metadata, perm, scheme = "free") {
  if (perm == 0) return(metadata)
  set.seed(perm)
  shuffle <- function(m) m %>%
    group_by(pH) %>% mutate(ploidy = sample(ploidy)) %>% ungroup() %>%
    group_by(ploidy) %>% mutate(pH = sample(pH)) %>% ungroup()
  if (scheme == "free") return(shuffle(metadata))
  fixed <- metadata$sampleID %in% geneticOutliers #fixGeno4: outliers keep their labels; the rest are shuffled the same way
  metadata[!fixed, ] <- shuffle(metadata[!fixed, ])
  metadata
}

# Region-level power simulation. Same binomial thinning as spikeIn() in 04.4-DSS-slurm.R, applied to each region's
# reads summed over its CpGs: for a spiked region with pooled methylation m (all samples, ignoring labels), oysters on the
# "up" side have each unmethylated read switched to methylated with probability (delta / 2) / (1 - m), and oysters on the
# "down" side have each methylated read switched with probability (delta / 2) / m. Coverage is unchanged. Regions are
# chosen with the SIM_REP seed from nuclear regions with at least LO_COUNT reads in every sample (so they pass the
# coverage filter and are tested with min.per.group = All) and pooled methylation of 15-85% (so both probabilities are at
# most 1 for delta up to 30). The other seed depends only on SIM_REP and SIM_DELTA, so every array task gets the same data
simPool <- c(0.15, 0.85)

setColumn <- function(x, column, value) { #Replace one column of a methylRaw (an S4 object that extends data.frame)
  x@.Data[[match(column, x@names)]] <- value
  x
}

spikeRegions <- function(regionCounts, metadata, delta, nSpike, rep, minReads) {
  data <- lapply(regionCounts, methylKit::getData)
  keys <- lapply(data, function(d) paste(d$chr, d$start, d$end))
  pool <- Reduce(intersect, lapply(seq_along(data), function(j) keys[[j]][data[[j]]$coverage >= minReads & data[[j]]$chr != mitoChr]))
  rows <- lapply(keys, function(k) match(pool, k))
  m <- Reduce(`+`, lapply(seq_along(data), function(j) data[[j]]$numCs[rows[[j]]])) /
       Reduce(`+`, lapply(seq_along(data), function(j) data[[j]]$coverage[rows[[j]]]))
  inRange <- m >= simPool[1] & m <= simPool[2]
  pool <- pool[inRange]
  m <- m[inRange]
  stopifnot(length(pool) >= 2 * nSpike)
  set.seed(rep)
  chosen <- sample(length(pool), 2 * nSpike)
  spikeM <- m[chosen]
  spikes <- tibble(key = pool[chosen], meanMeth = 100 * spikeM,
                   term = rep(c("ploidy", "pH"), each = nSpike),
                   direction = sample(c(-1L, 1L), 2 * nSpike, replace = TRUE)) %>% #+1 = higher in triploids / at low pH
    separate(key, c("chr", "start", "end"), sep = " ", convert = TRUE, remove = FALSE)
  groupSign <- list(ploidy = ifelse(metadata$ploidy == "3N", 1L, -1L), pH = ifelse(metadata$pH == "low", 1L, -1L))
  set.seed(1000 * rep + delta)
  for (j in seq_along(data)) {
    r <- match(spikes$key, keys[[j]])
    up <- spikes$direction * ifelse(spikes$term == "ploidy", groupSign$ploidy[j], groupSign$pH[j]) > 0
    numCs <- data[[j]]$numCs
    numTs <- data[[j]]$numTs
    gained <- rbinom(length(r), numTs[r], (delta / 200) / (1 - spikeM))
    lost <- rbinom(length(r), numCs[r], (delta / 200) / spikeM)
    numCs[r] <- numCs[r] + ifelse(up, gained, -lost)
    regionCounts[[j]] <- setColumn(setColumn(regionCounts[[j]], "numCs", numCs), "numTs", data[[j]]$coverage - numCs)
  }
  list(regionCounts = regionCounts, spikes = dplyr::select(spikes, -key))
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
  smallestGroup <- min(table(if (test == "ploidy") sampleMetadata$ploidy[keep] else sampleMetadata$pH[keep]))
  if (presence != "All" && as.integer(presence) > smallestGroup) { #e.g. min.per.group 10 with only 8 high-pH oysters in DropGeno4
    message("Skipping: min.per.group = ", presence, " but the smallest ", test, " group has ", smallestGroup, " oysters")
    quit(save = "no", status = 0)
  }
  for (d in file.path(runDir, c("counts", "rds", "pvalue-histograms"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  metadata <- permuteLabels(sampleMetadata[keep, ], perm, permScheme)
  treatment <- if (test == "ploidy") as.integer(metadata$ploidy == "3N") else as.integer(metadata$pH == "low")
  covariates <- if (test == "ploidy") data.frame(pH = metadata$pH) else data.frame(ploidy = metadata$ploidy) #Other factor is the covariate

  regionCounts <- readRDS(file.path(countDir, "region-counts.rds")) %>%
    methylKit::reorganize(sample.ids = metadata$sampleID, treatment = treatment)
  if (simDelta > 0) {
    stopifnot(presence == "All")
    sim <- spikeRegions(regionCounts, metadata, simDelta, simN, simRep, loCount)
    regionCounts <- sim$regionCounts
    write_csv(sim$spikes, file.path(runDir, "spiked-regions.csv"))
    message("Spiked ", simDelta, "% into ", simN, " regions per term (replicate ", simRep, ")")
  }
  regionCounts <- regionCounts %>%
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
                   perm.scheme = ifelse(perm == 0, "observed", permScheme),
                   test = test, min.per.group = presence,
                   regions.tested = nrow(results), regions.tested.mito = sum(isMito)) %>%
    bind_cols(counts)
  print(as.data.frame(counts))

  taskTag <- paste0(test, "-", presence)
  write_tsv(counts, file.path(runDir, "counts", paste0(taskTag, ".tsv")))
  pBins <- cut(results$pvalue, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)
  write_tsv(as_tibble(table(bin = pBins)), file.path(runDir, "pvalue-histograms", paste0(taskTag, ".tsv"))) #Calibration check (plan E1)

  if (simDelta > 0) { #Power simulation: per-region results only
    saveRDS(dplyr::select(results, chr, start, end, meth.diff, pvalue, qvalue),
            file.path(runDir, "rds", paste0("diffMeth-", taskTag, ".rds")), compress = FALSE)
  } else if (perm == 0) { #Only keep full results and BED files for the observed labels
    saveRDS(results, file.path(runDir, "rds", paste0("diffMeth-", taskTag, ".rds")), compress = FALSE)
    dmr <- results %>%
      filter(qvalue < primaryQ, abs(meth.diff) >= primaryDiff) %>%
      transmute("#chr" = chr, start = start - 1L, end, across(any_of("regionID")), meth.diff, pvalue, qvalue, qvalueBH) %>% #0-based BED. Header starts with # so bedtools skips it
      arrange(`#chr`, start)
    write_tsv(dmr, file.path(runDir, paste0("DMR-", taskTag, "-q", primaryQ, "-diff", primaryDiff, ".bed")))
  } else if (!grepl("^tile", regionSet)) { #Per-region results for permuted labels, for enrichment tests against the same permutations. Tiles are too large to keep
    saveRDS(dplyr::select(results, chr, start, end, any_of("regionID"), meth.diff, pvalue),
            file.path(runDir, "rds", paste0("diffMeth-", taskTag, ".rds")))
  } else { #Tiles: p-values only, as the null for the power simulation's empirical FDR
    saveRDS(sort(results$pvalue), file.path(runDir, "rds", paste0("pvalues-", taskTag, ".rds")))
  }
}

# Stage 6: combine counts from every run in this directory (observed and permuted), and compare each observed
# count with its permutation distribution (plan E2)

if (mode == "summary") {
  countFiles <- list.files(".", pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE) %>%
    str_subset("^\\./DMR-[^/]+/counts/") %>%
    str_subset("-simD[0-9]+-n[0-9]+-rep[0-9]+/", negate = TRUE) #Power-simulation runs are summarized by simsummary
  stopifnot(length(countFiles) > 0)
  counts <- map_dfr(countFiles, read_tsv, col_types = cols(.default = "c", perm = "i", cov.bases = "i", lo.count = "i",
                                                           regions.tested = "i", regions.tested.mito = "i",
                                                           qvalue = "d", difference = "d", DMR = "i", hypermethylated = "i",
                                                           hypomethylated = "i", DMR.mito = "i", DMR.BH = "i", primary = "l")) %>%
    mutate(perm.scheme = if ("perm.scheme" %in% names(.)) perm.scheme else NA_character_,
           perm.scheme = coalesce(perm.scheme, ifelse(perm == 0, "observed", "free"))) %>% #Runs from before PERM_SCHEME existed
    arrange(regions, cov.bases, lo.count, hi.perc, samples, overdispersion, test, min.per.group, qvalue, difference, perm.scheme, perm)
  write_csv(counts, "DMR-counts-all.csv")

  settingColumns <- c("method", "regions", "cov.bases", "lo.count", "hi.perc", "samples", "overdispersion",
                      "test", "min.per.group", "qvalue", "difference", "primary")
  observed <- counts %>% filter(perm == 0) %>% dplyr::select(all_of(settingColumns), regions.tested, DMR)
  permuted <- counts %>% #One comparison per permutation scheme
    filter(perm != 0) %>%
    group_by(across(all_of(c(settingColumns, "perm.scheme")))) %>%
    summarize(permutations = n(), perm.median = median(DMR), perm.95 = quantile(DMR, 0.95, names = FALSE),
              perm.DMR = list(DMR), .groups = "drop")
  permSummary <- left_join(observed, permuted, by = settingColumns) %>%
    mutate(empirical.p = map2_dbl(perm.DMR, DMR, ~ if (is.null(.x)) NA_real_ else (1 + sum(.x >= .y)) / (1 + length(.x))),
           above.perm.95 = DMR > perm.95) %>%
    dplyr::select(-perm.DMR)
  write_csv(permSummary, "DMR-permutation-summary.csv")
  print(as.data.frame(filter(permSummary, primary)))
}

# Stage 7: GO gene-set enrichment of the gene-level results (plan F), for ploidy and pH with coverage in all samples.
# Genes are ranked by signed -log10 p (sign of meth.diff: positive = higher methylation in triploids or at low pH),
# and each GO set gets a GSEA enrichment score (weighted Kolmogorov-Smirnov running sum, weight = |statistic|).
# Significance comes from the label permutations, not from methylKit p-values, which are anti-conservative here:
#   NES = ES / mean ES of the same sign over the permutations of the same set (GSEA phenotype-permutation normalization)
#   nominal p = share of same-sign permutation ES at least as extreme
#   FDR = share of permutation NES at least as extreme, pooled over all sets, divided by the same share for observed NES
# Each permutation is also scored against the others, to show how many sets the procedure flags without any real signal

if (mode == "enrich") {
  stopifnot(regionSet == "gene", file.exists(gafFile))
  suppressPackageStartupMessages({
    library(GO.db)
    library(AnnotationDbi)
  })
  enrichDir <- paste0("enrichment-", runTag, permTag)
  dir.create(enrichDir, showWarnings = FALSE)
  minSize <- 10
  maxSize <- 500

  # GO annotation: direct terms from the GAF plus all their ancestors (GO.db), so each gene is in every term above its own
  gaf <- fread(gafFile, header = FALSE, sep = "\t", quote = "", skip = "NCBIGene", select = c(3, 4, 5, 9),
               col.names = c("symbol", "qualifier", "GO", "aspect"))
  gaf <- unique(gaf[!grepl("^NOT", qualifier), .(symbol, GO, aspect)])
  ancestorMaps <- list(P = as.list(GOBPANCESTOR), F = as.list(GOMFANCESTOR), C = as.list(GOCCANCESTOR))
  ontologyName <- c(P = "BP", F = "MF", C = "CC")
  annotation <- rbindlist(lapply(names(ancestorMaps), function(a) {
    direct <- gaf[aspect == a]
    known <- direct$GO %in% names(ancestorMaps[[a]]) #Drops terms obsolete in this GO.db release
    message(ontologyName[a], ": ", sum(!known), " of ", nrow(direct), " annotations use GO IDs not in GO.db, dropped")
    direct <- direct[known]
    expanded <- data.table(symbol = rep(direct$symbol, lengths(ancestorMaps[[a]][direct$GO]) + 1L),
                           GO = unlist(Map(c, direct$GO, ancestorMaps[[a]][direct$GO]), use.names = FALSE))
    unique(expanded[GO != "all"])[, ontology := ontologyName[a]]
  }))

  # Signed -log10 p for each gene. Genes with the same coordinates share a region, so each of their IDs gets its result
  readGeneStats <- function(file) {
    d <- as.data.table(readRDS(file))[, .(regionID, meth.diff, pvalue)]
    d <- d[, .(symbol = sub("^gene-", "", unlist(strsplit(regionID, ",")))), by = .(regionID, meth.diff, pvalue)]
    d[, stat := sign(meth.diff) * -log10(pmax(pvalue, .Machine$double.xmin))]
    d[, .(symbol, stat)]
  }

  # GSEA enrichment score for one ranked statistic (sorted decreasing) and the sorted positions of one set's genes
  enrichmentScore <- function(rankedStat, positions) {
    n <- length(rankedStat)
    k <- length(positions)
    hitWeights <- abs(rankedStat[positions])
    hitSum <- cumsum(hitWeights) / sum(hitWeights)
    missSum <- (positions - seq_len(k)) / (n - k) #Misses before each hit
    top <- max(hitSum - missSum) #Running sum just after each hit
    bottom <- min(0, c(0, hitSum[-k]) - missSum) #Running sum just before each hit
    if (top >= -bottom) top else bottom
  }

  for (test in c("ploidy", "pH")) {
    taskTag <- paste0(test, "-All")
    observedFile <- file.path(paste0("DMR-", runTag), "rds", paste0("diffMeth-", taskTag, ".rds"))
    permFiles <- Sys.glob(file.path(paste0("DMR-", runTag, permTag, "-perm*"), "rds", paste0("diffMeth-", taskTag, ".rds")))
    permIDs <- as.integer(str_match(permFiles, "-perm([0-9]+)/")[, 2])
    permFiles <- permFiles[order(permIDs)]
    message(test, ": ", length(permFiles), " permutations")
    stopifnot(file.exists(observedFile), length(permFiles) >= 20)

    observed <- readGeneStats(observedFile)
    universe <- unique(observed$symbol)
    statList <- c(list(observed), lapply(permFiles, readGeneStats))
    statMatrix <- sapply(statList, function(d) { #Genes x rankings; column 1 is the observed labels
      stopifnot(setequal(unique(d$symbol), universe)) #Same genes in every run (presence "All" doesn't depend on labels)
      d$stat[match(universe, d$symbol)]
    })

    # GO sets within the tested genes, kept if they have minSize-maxSize genes. Sets with identical genes are tested once
    sets <- annotation[symbol %in% universe]
    sets <- sets[, .(genes = list(sort(unique(symbol)))), by = .(GO, ontology)]
    sets[, size := lengths(genes)]
    sets <- sets[size >= minSize & size <= maxSize]
    sets[, geneKey := map_chr(genes, paste, collapse = ",")]
    sets <- sets[, .(GO = GO[1], ontology = ontology[1], genes = genes[1], size = size[1], sameGenesAs = paste(GO[-1], collapse = ",")), by = geneKey][, geneKey := NULL]
    message(test, ": ", nrow(sets), " GO sets with ", minSize, "-", maxSize, " genes, from ", length(universe), " ranked genes (",
            sum(universe %in% annotation$symbol), " with GO terms)")

    # ES for every set and ranking
    setIndex <- lapply(sets$genes, match, table = universe)
    es <- matrix(NA_real_, nrow(sets), ncol(statMatrix))
    for (j in seq_len(ncol(statMatrix))) {
      ord <- order(statMatrix[, j], decreasing = TRUE)
      rankPos <- integer(length(universe))
      rankPos[ord] <- seq_along(ord)
      ranked <- statMatrix[ord, j]
      es[, j] <- vapply(setIndex, function(i) enrichmentScore(ranked, sort(rankPos[i])), numeric(1))
    }

    # Normalize, p, and FDR for one ranking (column) against a set of null columns
    scoreRanking <- function(j, nullCols) {
      obsES <- es[, j]
      nullES <- es[, nullCols, drop = FALSE]
      posMean <- rowMeans(ifelse(nullES >= 0, nullES, NA), na.rm = TRUE)
      negMean <- -rowMeans(ifelse(nullES < 0, nullES, NA), na.rm = TRUE)
      scale <- ifelse(obsES >= 0, posMean, negMean)
      nes <- obsES / scale
      nullNES <- ifelse(nullES >= 0, nullES / posMean, nullES / negMean)
      nominalP <- vapply(seq_along(obsES), function(i) {
        sameSign <- if (obsES[i] >= 0) nullES[i, nullES[i, ] >= 0] else nullES[i, nullES[i, ] < 0]
        if (length(sameSign) == 0) NA_real_ else (1 + sum(abs(sameSign) >= abs(obsES[i]))) / (1 + length(sameSign))
      }, numeric(1))
      pooledNull <- as.vector(nullNES)
      pooledNull <- pooledNull[is.finite(pooledNull)]
      fdr <- vapply(nes, function(x) {
        if (!is.finite(x)) return(NA_real_)
        if (x >= 0) {
          (mean(pooledNull[pooledNull >= 0] >= x)) / mean(nes[nes >= 0 & is.finite(nes)] >= x)
        } else {
          (mean(pooledNull[pooledNull < 0] <= x)) / mean(nes[nes < 0 & is.finite(nes)] <= x)
        }
      }, numeric(1))
      tibble(ES = obsES, NES = nes, nominalP = nominalP, FDR = pmin(fdr, 1))
    }

    nRankings <- ncol(es)
    result <- bind_cols(as_tibble(sets[, .(GO, ontology, size, sameGenesAs)]), scoreRanking(1, 2:nRankings)) %>%
      mutate(term = suppressMessages(AnnotationDbi::Term(GO)),
             meanStat = map_dbl(setIndex, ~ mean(statMatrix[.x, 1]))) #Mean signed -log10 p of the set's genes
    observedDiff <- as.data.table(readRDS(observedFile))[, .(symbol = sub("^gene-", "", unlist(strsplit(regionID, ",")))), by = .(regionID, meth.diff)]
    result <- result %>%
      mutate(meanMethDiff = map_dbl(sets$genes, ~ mean(observedDiff$meth.diff[match(.x, observedDiff$symbol)])),
             direction = ifelse(NES >= 0, "higher in treatment", "lower in treatment")) %>%
      relocate(term, .after = GO) %>%
      arrange(FDR, nominalP)
    write_csv(result, file.path(enrichDir, paste0("GSEA-GO-", taskTag, ".csv")))

    # Same procedure with each permutation as the "observed" ranking and the other permutations as the null
    calibration <- map_dfr(2:nRankings, function(j) {
      s <- scoreRanking(j, setdiff(2:nRankings, j))
      tibble(ranking = paste0("perm", sort(permIDs)[j - 1]), FDR.05 = sum(s$FDR < 0.05, na.rm = TRUE),
             FDR.25 = sum(s$FDR < 0.25, na.rm = TRUE), nominalP.01 = sum(s$nominalP < 0.01, na.rm = TRUE))
    })
    observedCounts <- tibble(ranking = "observed", FDR.05 = sum(result$FDR < 0.05, na.rm = TRUE),
                             FDR.25 = sum(result$FDR < 0.25, na.rm = TRUE), nominalP.01 = sum(result$nominalP < 0.01, na.rm = TRUE))
    calibration <- bind_rows(observedCounts, calibration) %>% mutate(test = test, sets = nrow(sets), .before = 1)
    write_csv(calibration, file.path(enrichDir, paste0("GSEA-GO-", taskTag, "-calibration.csv")))
    message(test, ": observed sets at FDR < 0.05 / 0.25 = ", observedCounts$FDR.05, " / ", observedCounts$FDR.25,
            "; permutations scored the same way: median ", median(calibration$FDR.05[-1]), " / ", median(calibration$FDR.25[-1]),
            ", max ", max(calibration$FDR.05[-1]), " / ", max(calibration$FDR.25[-1]))
    print(head(as.data.frame(dplyr::select(result, GO, term, ontology, size, NES, nominalP, FDR, meanMethDiff)), 15))
  }
}

# Stage 8: region-level power simulation. For each SIM_DELTA > 0 run of this setting (REGIONS, SAMPLES = All, coverage in
# all samples), how many spiked regions are found, and would the experiment as a whole have looked different from the
# label-permutation null (free scheme, same setting)? Detection rules, as in the 04.4 DSS simulation:
#   1. DMR: methylKit q < 0.05 and |meth.diff| >= 10 (the fixed DMR definition); also q < 0.05 at any |meth.diff|
#   2. Empirical FDR < 0.05: the largest p-value cutoff at which (mean permuted count below it) / (simulated count below
#      it) is at most 0.05
#   3. Experiment-level: the run's DMR count (rule 1) is above the 95th percentile of the permuted counts; also for the
#      q-only count

if (mode == "simsummary") {
  stopifnot(samples == "All", perm == 0)
  escaped <- gsub("\\.", "\\\\.", paste0("DMR-", runTag))
  simDirs <- list.files(".", pattern = paste0("^", escaped, "-simD[0-9]+-n[0-9]+-rep[0-9]+$"))
  permDirs <- list.files(".", pattern = paste0("^", escaped, "-perm[0-9]+$"))
  terms <- c("ploidy", "pH")
  readPerm <- function(d, t) { #Sorted p-values of one permuted run, or NULL
    f <- file.path(d, "rds", paste0(c("pvalues-", "diffMeth-"), t, "-All.rds"))
    if (file.exists(f[1])) return(readRDS(f[1]))
    if (file.exists(f[2])) return(sort(readRDS(f[2])$pvalue))
    NULL
  }
  pGrid <- 10^seq(-12, -1, by = 0.02)
  permNull <- map(setNames(terms, terms), function(t) {
    ps <- compact(map(permDirs, readPerm, t = t))
    counts <- map_dfr(permDirs, function(d) {
      f <- file.path(d, "counts", paste0(t, "-All.tsv"))
      if (file.exists(f)) read_tsv(f, show_col_types = FALSE) else NULL
    })
    list(meanBelow = colMeans(do.call(rbind, map(ps, ~ findInterval(pGrid, .x, left.open = TRUE)))),
         nP = length(ps),
         dmr = counts %>% filter(qvalue == 0.05, difference == primaryDiff) %>% pull(DMR),
         qOnly = counts %>% filter(qvalue == 0.05, difference == 0) %>% pull(DMR))
  })
  message(length(simDirs), " simulation runs; permutations with p-values: ",
          paste(terms, map_int(permNull, "nP"), collapse = ", "), "; with counts: ", paste(terms, map_int(permNull, ~ length(.x$dmr)), collapse = ", "))
  stopifnot(length(simDirs) > 0, all(map_int(permNull, "nP") >= 20))

  runs <- map_dfr(simDirs, function(d) {
    spikes <- read_csv(file.path(d, "spiked-regions.csv"), show_col_types = FALSE)
    map_dfr(terms, function(t) {
      f <- file.path(d, "rds", paste0("diffMeth-", t, "-All.rds"))
      if (!file.exists(f)) return(NULL)
      r <- readRDS(f)
      s <- spikes %>% filter(term == t)
      idx <- match(paste(s$chr, s$start, s$end), paste(r$chr, r$start, r$end))
      found <- !is.na(idx)
      isSpike <- seq_len(nrow(r)) %in% idx
      dmr <- r$qvalue < primaryQ & abs(r$meth.diff) >= primaryDiff
      qHit <- r$qvalue < primaryQ
      simBelow <- findInterval(pGrid, sort(r$pvalue), left.open = TRUE)
      nul <- permNull[[t]]
      ok <- simBelow > 0 & nul$meanBelow / pmax(simBelow, 1) <= 0.05
      cutoff <- if (any(ok)) max(pGrid[ok]) else 0
      tibble(run = d,
             delta = as.integer(sub(".*-simD([0-9]+)-.*", "\\1", d)),
             nSpiked = as.integer(sub(".*-n([0-9]+)-rep.*", "\\1", d)),
             rep = as.integer(sub(".*-rep([0-9]+)$", "\\1", d)),
             term = t, regionsTested = nrow(r), spikedTested = mean(found),
             realizedDelta = mean(s$direction[found] * r$meth.diff[idx[found]]),
             powerDMR = sum(dmr[idx[found]]) / nrow(s),
             powerQ = sum(qHit[idx[found]]) / nrow(s),
             correctSign = mean(sign(r$meth.diff[idx[found]][qHit[idx[found]]]) == s$direction[found][qHit[idx[found]]]),
             falseDMR = sum(dmr & !isSpike),
             empiricalCutoff = cutoff,
             powerEmpiricalFDR = sum(r$pvalue[idx[found]] < cutoff) / nrow(s),
             totalDMR = sum(dmr), permDMR95 = quantile(nul$dmr, 0.95, names = FALSE),
             experimentDetectedDMR = sum(dmr) > quantile(nul$dmr, 0.95, names = FALSE),
             totalQ = sum(qHit), permQ95 = quantile(nul$qOnly, 0.95, names = FALSE),
             experimentDetectedQ = sum(qHit) > quantile(nul$qOnly, 0.95, names = FALSE))
    })
  })

  powerSummary <- runs %>%
    group_by(nSpiked, delta, term) %>%
    summarize(reps = n(), realizedDelta = mean(realizedDelta), spikedTested = mean(spikedTested),
              across(c(powerDMR, powerQ, powerEmpiricalFDR, correctSign, falseDMR, totalDMR, totalQ), ~ mean(.x, na.rm = TRUE)),
              fractionExperimentDetectedDMR = mean(experimentDetectedDMR), fractionExperimentDetectedQ = mean(experimentDetectedQ),
              permDMR95 = dplyr::first(permDMR95), permQ95 = dplyr::first(permQ95), .groups = "drop") %>%
    mutate(regions = regionSet, .before = 1) %>%
    arrange(nSpiked, factor(term, levels = terms), delta)
  simDir <- "power-simulation"
  dir.create(simDir, showWarnings = FALSE)
  write_csv(runs, file.path(simDir, paste0("power-runs-", runTag, ".csv")))
  write_csv(powerSummary, file.path(simDir, paste0("power-summary-", runTag, ".csv")))

  p <- runs %>%
    transmute(nSpiked, delta, term,
              `DMR (q < 0.05, |diff| >= 10%)` = powerDMR, `q < 0.05, any |diff|` = powerQ, `Empirical FDR < 0.05` = powerEmpiricalFDR,
              `Experiment-level, DMR count > permutation 95th pct` = as.numeric(experimentDetectedDMR),
              `Experiment-level, q-only count > permutation 95th pct` = as.numeric(experimentDetectedQ)) %>%
    pivot_longer(-c(nSpiked, delta, term), names_to = "rule", values_to = "power") %>%
    group_by(nSpiked, delta, term, rule) %>% summarize(power = mean(power), .groups = "drop") %>%
    ggplot(aes(delta, power, color = rule)) +
    geom_line() + geom_point() +
    facet_grid(rows = vars(paste(nSpiked, "spiked regions")), cols = vars(term)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Spiked difference (percentage points)", y = "Power (fraction of spiked regions found, or of replicates detected)",
         color = NULL, subtitle = paste0(runTag, "; null from ", permNull$ploidy$nP, " permutations")) +
    theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave(file.path(simDir, paste0("power-curves-", runTag, ".png")), p, width = 8, height = 6.5, dpi = 150)
  print(as.data.frame(powerSummary))
}

sessionInfo()
