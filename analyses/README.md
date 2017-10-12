# `analyses` Subdirectory Structure

This subdirectory hosts the output of all analyses regarding *C. gigas* and ocean acidification. Every folder houses all necessary files for a discrete analysis step in a larger pipeline. Folders follow this naming convention:

**ProjectID_AnalysisStep_DateInitiated** *(ex. DNR_SRM_20170902)*

The following pipelines are represented in this subdirectory:

## DNR

**Environmental Data**:

- [Summary Information and Visualizations](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Environmental_Data_Analysis_20161115)

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
  - [Boxplots](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-09-12-Protein-Area-Boxplots): Visualized protein area at different sites and eelgrass habitats
  - [Transition Replicate Correlations](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902/2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations): Regressed technical replicates against eachother to identify transitions with poor technical replication

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


## Manchester

No analyses as of 2017-09-27.