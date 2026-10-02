# Haws ploidy × pH methylation: where things stand (2026-10-02)

Summary of the differential-methylation work in 04.2–04.6 (PRs #48–#55) and options for what comes next. Details are in [04.4-DSS-plan.md](04.4-DSS-plan.md) and [04.6-methylKit-DMR-results.md](04.6-methylKit-DMR-results.md).

## Bottom line

With 6 oysters per ploidy × pH group, **no ploidy or pH effect on methylation can be detected** at any scale: single CpGs, 250 bp or 1 kb windows, whole genes, or GO gene sets. What does stand out is that **methylation tracks each oyster's genetic background**.

## What was done

| Step | What | Status |
|---|---|---|
| 04.2 / 04.3 methylKit per-CpG | Fixed argument and outlier bugs; Slurm version comparing overdispersion settings | Merged |
| 04.4 DSS per-CpG | New pipeline: ploidy × pH model with the interaction kept as a nuisance term, mito kept in, sample QC, label permutations, SNP check, genotype covariate, runs with outlier oysters dropped | Merged |
| 04.6 methylKit regions | Plan and pipeline: windows, genes, a new TE file, GO enrichment, runs with genetic outliers dropped, a permutation that keeps the four outliers at high pH (`fixGeno4`) | Merged |

## Key findings

1. **Overdispersion decides everything.** Without it (methylKit `none`), there are thousands of hits: 730–6,896 per-CpG DML and 4,500–9,000 DMRs. With it (`MN`), there are about 0. The `none` hits are differences between individual oysters, not treatment effects.
2. **Shuffled labels give as many hits as the real ones.**
   - DSS: 187 ploidy DML from the real labels vs. a mean of 223 from 100 shuffles. pH gave 185 vs. 166.
   - methylKit regions: in all 72 comparisons the real count is inside the shuffled distribution.
   - So the reported FDRs from both tools are too optimistic for this data. Only the permutations give a fair null.
3. **All 24 samples pass technical QC** (conversion 99.5%, mapping 60–62%). 2H-3, one of the two oysters dropped in earlier analyses, isn't an outlier.
4. **Genetic background is the one clear result.**
   - Genetic similarity between pairs of oysters predicts their methylation similarity (rho = 0.86). This holds with only SNPs that bisulfite conversion can't fake, and still holds at 0.67 among the 20 non-outlier oysters.
   - Four oysters are genetically distinct (2H-1, 2H-2, 3H-5, 3H-2), and **all four are high pH**. That would happen by chance about 5% of the time, so in the full 24-oyster analysis pH is confounded with genotype.
5. **Accounting for genotype doesn't uncover a hidden effect.**
   - Dropping the outliers (three different ways) still gives counts within the permutation distribution.
   - So does the `fixGeno4` permutation.
   - The only leads were 9 pH genes and one GO set (mitotic spindle checkpoint, FDR 0.078). Neither replicated once all four outliers were removed or the genotype-matched permutation was used.

## Housekeeping for the paper

- **Earlier results are probably false positives.** The old 04-methylKit and 04-DSS DML (178 ploidy / 154 pH DML), and the 07–09 annotation and enrichment built on them, came from settings now shown to give false positives. Anything citing them needs to be flagged or retracted.
- **07's TE overlaps used the wrong file.** That file is off by one base and includes simple repeats, so those overlaps need redoing (if they're kept at all).
- **Pull the hatchery records for the four genetic outliers.** If they share a source or family, the pH arm has a design confound that should be reported.

## Options from here

**A. Write it up as a solid null with a positive genetics result (recommended).** The paper would say that at this sample size, ploidy and pH have no detectable locus-specific effect on methylation, and that genetic background explains much of the variation between oysters. To make the null convincing:

- **A power simulation:** add known Δ effects into the real data, rerun DSS/MN with permutations, and report the smallest effect that could have been detected. **Done for DSS per-CpG ([04.4-DSS-power-simulation.md](04.4-DSS-power-simulation.md)):** an effect of about ≥ 20 points at 1,000 CpGs would have been detected against the permutation null; ≤ 15 points would not. Finding the individual CpGs at permutation-calibrated FDR takes about 30–40 points. **Region level done ([04.6-methylKit-DMR-power-simulation.md](04.6-methylKit-DMR-power-simulation.md)):** gene bodies are much more sensitive. A 5-point shift across ~300 genes (1%) would have been detected every time, and at 10 points about two-thirds of those genes are found individually. For 1 kb windows the limit is about 10 points. So the real null rules out gene-body effects of about ≥ 5 points across ~1% of genes.
- **Variance partitioning:** a whole-methome PERMANOVA or similar on methylation distance, with genotype PCs, ploidy and pH as terms. Pair it with the 06 global-methylation results. **Done ([04.7-variance-partitioning-plan.md](04.7-variance-partitioning-plan.md)):** genotype (3 PCs) explains 8–12% of the methylation distance between oysters (Mantel r = 0.73–0.77), and 5–9% without the four outliers (r = 0.28–0.41). Ploidy and pH explain −0.5% to 0.4%, and per CpG or gene they are within their shuffle nulls. Genotype-associated features are 4–5× above chance (2× without the outliers).

**B. Close out the planned tests, then stop.** The exon, intron, upstream and TE region sets are built but not tested. Run each once against permutations so the region analysis is complete. Don't widen the parameter grid: about 200 settings have already been tried, and more searching mostly adds chances for false positives. Drop shrinkMN, DSS smoothing and group-blind segmentation unless a reviewer asks.

**C. Reframe around genetics (more ambitious).** **Started ([04.8-cis-genotype-plan.md](04.8-cis-genotype-plan.md)):** methylation is linked to the genotype of SNPs within 50 kb, beyond genome-wide relatedness. Over a same-chromosome, > 1 Mb null, enrichment is about 1.7× overall, 2.8× within 250 bp, and 1.9× for gene bodies, with a similar pattern without the outliers. LD decay shows that linkage between SNPs is short (clearly elevated r² only within ~1 kb), plus weak chromosome-scale relatedness. So the peak within 1 kb is consistent with local (cis) control of methylation. The flat ~1.6× out to 50 kb isn't explained by LD, and is more likely regional methylation domains. Methylation–methylation correlation decay is the next check. The rho = 0.86 result suggests SNPs that drive nearby methylation (cis-mQTL-like) or methylation that follows family. With 24 oysters this is descriptive only. It could seed a follow-up proposal or connect to other lab datasets that have genotypes.

**D. Get more power (only if the biology needs a positive answer).** More oysters per cell, or a targeted assay on candidate genes, rather than more analysis of this dataset. Any new sampling should balance families across treatments. If RNA-seq exists for these oysters, testing expression against genotype is a cheaper next step.

**Suggested order:** hatchery check → power simulation and variance partitioning (A) → B in the background → decide between C and D based on what the paper needs.
