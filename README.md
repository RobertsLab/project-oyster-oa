# project-oyster-oa

This repository focuses on efforts examining the impact of ocean acidification on Pacific oysters, *Crassostrea gigas*.

## Projects

Each project has a distinct ID used in all folder naming conventions to make sorting through the repository easier.

### Exploring Proteomic Variation in Pacific Oysters (ID: DNR)

In June 2016, 150 sibling *C. gigas* were outplanted at five different locations in Washington (leads: Micah Horowith (DNR) and Alex Lowe (UW Biology)). Four of these sites were located in Puget Sound — Case Inlet, Fidalgo Bay, Port Gamble Bay, Skokomish River Delta — and another in southern Washington, Willapa Bay. At each site, oysters were placed in both bare and eelgrass patches. After one month, gill tissue samples were collected and flash frozen. [Yaamini Venkataraman](yaaminiv.github.io) began proteomic analyses on samples to see if environmental variables like pH, dissolved oxygen, and temperature at each site and habitat drove differential protein expression.

Shotgun proteomic methods were used to characterize expression of 9,047 proteins present in ten ctenidia samples in January 2017. Each sample represented a unique site and habitat condition. A subset of 273 environmental response proteins were identified. Proteins within this subset were evaluated as potential biomarkers. Fifteen proteins with varying functions (environmental response, metabolism and growth) were selected as targets for a Selected Reaction Monitoring assay in July 2017. The efficacy of these proteins as early-response bioindicators will be evaluated using multivariate statistical methods.

### Pacific Oyster Ocean Acidification Trials at Chew Hatchery (ID: Manchester)

[Yaamini Venkataraman](yaaminiv.github.io) exposed adult *C. gigas* from Willapa Bay to either low and ambient pH conditions for seven weeks (February 4, 2017 to April 8, 2017). Oysters were separated into six tanks: Tanks 1-3 experienced low pH conditions, and Tanks 4-6 experienced ambient pH conditions. Water samples from each experimental tank and the two header tanks were collected twice a week and poisoned with mercuric chloride for chemistry analysis. Before and after pH exposure, ctenida, mantle and adductor tissues were sampled from twenty individuals and flash frozen. Additional ctenida tissue was placed in ethanol for downstream DNA analyses, and gonad tissue was used for histological analyses to evaluate maturation. Live weight, empty shell weight and shell lengths were also recorded. This trial will elucidate how a single stressor affects methylation in different tissues.

After exposure, oysters were pooled and redistributed into two 100 L tanks with ambient water conditions. On June 5, 2017, a subset of oysters were exposed to a one hour heat shock of 40ºC on June 5, 2017. Oysters were placed back in ambient conditions. Starting June 30, 2017, oysters began conditioning for spawning. Temperature was ramped up to 23ºC and algae dosing increased compensatorily. Towards the end of the conditioning period, gonad maturation was evaluated by sampling a few oysters from both holding tanks.

Oysters were strip spawned on July 30, 2017. Individual males were with female pools. Egg pools were determined based on life history. There were three resulting pools: low pH females, ambient pH females and heat shock females. Individual crosses were fertilized together, then distributed into two replicates. After 24 hours, hatch rate was calculated. Based on life history of both parents in a cross, d-hinge larvae were pooled together and redistributed to five replicates. Larvae that were not stocked were frozen for genetic analyses. There were five treatment pools: low pH female x low pH male, low pH female x ambient pH male, ambient pH female x low pH male, ambient pH female x ambient pH male, and heat shock female x heat shock female. Larvae stopped growing after reaching 60 microns. On August 15, 2017, 10,000 larvae from each of the five treatment groups were sampled and frozen for genetic analyses. The experiment was terminated on August 17, 2017 and all remaining larvae were collected and frozen. This trial will allow us to understand how parental life history affects larval production and differential mortality.

## Repository Contents

Each larger folder also has its own README.md file with additional information. With the exception of `analyses`, all folder contents are divided by project ID. Some files in [`.gitignore`](https://github.com/RobertsLab/project-oyster-oa/blob/master/.gitignore); highlights are outlined below.

### [analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses)

Output for multiple analyses. Project IDs are detailed in the folder names themselves.

**DNR**

- [BCA Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_BCA_Analysis)
- [Environmental Data Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Environmental_Data_Analysis_20161115)
- [MZ Ratios for Mass Spectrometry](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MZratios_Larval_Samples_20170118)
- [MSConvert](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MSConvert_20170412)
- [Preliminary DIA Analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Preliminary_Analyses_20170321)
- [MSStats Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170524)
- [SRM Transition Selection](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_TransitionSelection_20170707)
- [SRM Analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902)
  - Several subdirectories are also housed here relating to troubleshooting SRM data

**Manchester**

- [Gonad Histology](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/Manchester_Gonad_Histology)
- [Reproductive Output](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/Manchester_ReproductiveOutput_20180214)
- [Larval Mortality](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/Manchester_Larval_Mortality_20180329)

### [data](https://github.com/RobertsLab/project-oyster-oa/tree/master/data)

Raw data used for project analyses.

**[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/DNR)**

- [*C. gigas* Sample Key](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/PugetSound-DNR-samples-2016.csv)
- Environmental Data
	- [Temperature, dissolved oxygen and pH data](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-14-Environmental-Data-from-Micah.csv)
	- [Salinity data](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-25-Calculated-Salinity-Output-from-Micah.csv)
	- [Tidal Data](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-12-13-Tidal-Data-by-Site.csv)
	- [Environmental Summary Data](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Environmental%20Summary%20Data%20for%20Proteomics%20Project.csv)
- Biomarkers
	- [Biomarker data for Yaamini samples](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-21-Alex-Data-Yaamini-Samples-Only.csv)
	- [Oysters used for biomarkers by Alex Lowe](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Biomarker_proteomeCrossRef.xlsx)

**[Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester)**

- [Adult *C. gigas* Tissue Sampling Key](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae)
- [Pacific Oyster Larvae Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling)
- [Water Chemistry Information](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/Water-Chemistry-Data)

### [images](https://github.com/RobertsLab/project-oyster-oa/tree/master/images)

Relevant photos for downstream analysis.

**[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR)**

- [General Lab Notebook](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR/Lab-Notebook)
- [Skyline Cleaning Screenshots](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR/2017-08-31-Skyline-Cleaning-Screenshots)

**[Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester)**

- [February Oyster Sampling](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/2017-2-4-Oyster-Sampling)
- [*C. gigas* Gonad Histology](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-gonad-histology)
  - [Before OA experiment](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-gonad-histology/2017-02-04-Sampling)
  - [After OA experiment](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-gonad-histology/2017-04-08-Sampling)
- [*C. gigas* Larval Images](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-larvae)

### [miscellaneous](https://github.com/RobertsLab/project-oyster-oa/tree/master/miscellaneous)

Meeting notes and other background information

**[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/miscellaneous/DNR)**

**Manchester**

No contents as of 2017-09-27

### [notebooks](https://github.com/RobertsLab/project-oyster-oa/tree/master/notebooks)

Jupyter notebook that detail reproducible methods for data analysis

**[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/notebooks/DNR)**

- [Preliminary Proteomic Analyses](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-03-21-Preliminary-Proteomic-Data-Analyses.ipynb)
- [Demultiplexing Proteomic Data](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-04-12-Demultiplex-Raw-Files.ipynb)
- [Shotgun Proteomic Analyses](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-06-13-Full-Skyline-Preliminary-Analysis.ipynb)
- [DIA Analysis Pipeline](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-10-11-DIA-Skyline-Data-Pipeline.ipynb)
- [SRM Target Selection in MSStats](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-06-22-Selecting-SRM-Targets-with-MSstats-Part-2.ipynb)
- [Protein Annotations](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-07-05-Examining-Protein-Annotations.ipynb)
- [SRM Target ID in Skyline](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-07-07-SRM-Target-Identification-in-Skyline.ipynb)
- [SRM Analysis Pipeline](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-09-28-SRM-Skyline-Data-Pipeline.ipynb)

**Manchester**

No contents as of 2017-09-27

### [presentations](https://github.com/RobertsLab/project-oyster-oa/tree/master/presentations)

Final slides for various conferences

**[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/presentations/DNR)**

- [UW SAFS Graduate Student Symposium 2016](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/DNR/GSS2016_Venkataraman.pdf)
- [Pacific Coast Shellfish Growers Association 2017](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/DNR/PCSGA2017_Venkataraman.pptx)
- [UW SAFS Graduate Student Symposium 2017](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/DNR/PCSGA2017_Venkataraman.pptx)

**[Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/presentations/Manchester)**

- [National Shellfisheries Association 2018](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/Manchester/Venkataraman_NSA2018.pptx)