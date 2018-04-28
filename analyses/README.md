# `analyses` Subdirectory Structure

This subdirectory hosts the output of all analyses regarding *C. gigas* and ocean acidification. Every folder houses all necessary files for a discrete analysis step in a larger pipeline. Folders follow this naming convention:

**ProjectID_AnalysisStep_DateInitiated** *(ex. DNR_SRM_20170902)*

The following pipelines are represented in this subdirectory:

## DNR

**Environmental Data**:

- [Summary Information and Visualizations](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Environmental_Data_Analysis_20161115)
- [Environmental Data and Biomarker Analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-11-15-Environmental-Data-and-Biomarker-Analyses)
  - [Quality-Controlled Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-13-Environmental-Data-Quality-Control): Visualized environmental data after corrections for tidal data
  - [Biomarker Boxplots](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-27-Biomarker-Boxplots): Comparing biomarkers across sites
  - [Biomarker Scatterplots](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-11-29-Biomarker-Scatterplots): Comparing biomarker and environmental data
  - [Growth Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-11-15-Environmental-Data-and-Biomarker-Analyses/2017-12-19-Growth-Data-Analyses): Comparing growth patterns across sites

**Protein Extraction**:

- [BCA Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_BCA_Analysis): Used to quantify the amount of protein in individual oyster samples

**Shotgun Proteomics**:

- [MZ Ratios](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MZratios_Larval_Samples_20170118): Identification of MZ ratios for shotgun mass spectrometry methods based on previous runs from *C. gigas* larval samples
- [MSConvert](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MSConvert_20170412): Used to demultiplex RAW mass spectrometry files before using in Skyline
- [Protein Expression Characterization](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170524): Exploratory analyses for shotgun proteomics data
  - [Error checking](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170524/error-checking): Calculation of Skyline error rate
  - [MSStats](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170524/2017-06-22-MSstats): Determination of differentially expressed proteins between sites and habitats

**Selected Reaction Monitoring**:

- [Transition Selection](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_TransitionSelection_20170707): Identification of protein targets for Selected Reaction Monitoring Assay based on MSStats output
- [Target analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902): Multivariate analyses to determine targets differentially expressed based on sites and eelgrass habitat. Also includes data visualization instructions and data troubleshooting.
  - [Coefficient of Variation filtering](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-24-Coefficient-of-Variation): Used a CV benchmark to remove technical replicate data. Found that a threshold of CV > 20 provides results consistent with geoduck data
  - [Integrated Dataset](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset): Used CV filtering information to remove samples with a coefficient of variation greater than 20 for more than half of the transitions used in SRM. Integrated data at the peptide level and proceeded with NMDS, ANOSIM, boxplots, ANOVA, and Tukey Honest Significant Difference Tests
    - [Peptide-level Boxplots](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2017-11-06-Boxplots): Boxplots for each peptide across site only. ANOVA p-value on boxplots, Tukey HSD p-values written out in separate dataframe in this subdirectory
    - [Eelgrass-Unvegetated Differences by Site](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-Individual-Site-Eelgrass-Differences): Used to confirm that the global trend of no difference between habitat was not driven by one site
    - [DNR Paper Figure](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-11-05-Integrated-Dataset/2018-02-15-DNR-Paper-Figure): Used to create publication-quality figure

**Irrelevant Analyses**: Outputs that were not used further in the pipeline (i.e. failed attempts and dead ends)

- PECAN: Program used to create a `.blib` file for Skyline analyses. Contents of the folder were used in unsuccessful attempts, primarily due to large computing capacity necessary for this program. PECAN was used successfully by Emma Timmins-Schiffman for this project.
  - [Attempt 1](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_PECAN_20170228)
  - [Attempt 2](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_PECAN_RUN_2_20170307)
  - [Attempt 3](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_PECAN_Run_3_20170308)
- Preliminary Skyline Analysis: Protein expression data from Skyline using a `.blib` generated from *C. gigas* seed as opposed to adult *C. gigas* samples used for proteomics.
  - [Data from Skyline](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170314)
  - [Analysis products](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Preliminary_Analyses_20170321)

- Skyline
  - [Data Normalization](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170511): Conducted using inefficient version of `brecan`
  - [MSStats and Error Checking](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170512): Conducted using inefficient version of `brecan`
  
- [SRM Data Troubleshooting](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting)
  - [Boxplots](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-09-12-Protein-Area-Boxplots): Visualized protein area at different sites and eelgrass habitats
  - [Transition Replicate Correlations](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations): Regressed technical replicates against eachother to identify transitions with poor technical replication
  - [Confidence Interval filtering](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-24-Confidence-Interval-Transitions): Tried creating a confidence interval around a regression line to weed out bad samples
  - [PRTC](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-24-PRTC-NMDS): Made technical replicate NMDS plots with PRTC data, found that each run was most consistent with those that used the same batch of PRTC


## Manchester

**Histology**:

- [Pie Charts](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/Manchester_Gonad_Histology): Visualized maturation state and sex ratio before pH exposure and after pH expsoure (low pH and ambient pH groups)
- [Binomial GLM Analyses](https://github.com/RobertsLab/project-oyster-oa/blob/master/analyses/Manchester_Gonad_Histology/2018-02-27-Histology-Classification-Analyses.R): Analyzed differences in maturation states between pre and post-treatment sampling and pH treatments

**Reproductive Output**:

- [D-Hinge Abundance](https://github.com/RobertsLab/project-oyster-oa/blob/master/analyses/Manchester_ReproductiveOutput_20180214/2018-02-14-Reproductive-Output.R): Found significant difference in D-hinge abundance between maternal treatments. Egg production component is irrelevant.

**Larval Mortality**:

- [Abundance by Parental Treatment](https://github.com/RobertsLab/project-oyster-oa/blob/master/analyses/Manchester_Larval_Mortality_20180329/2018-03-29-Larval-Mortality.R): Looked at larval counts over time by parental treatment