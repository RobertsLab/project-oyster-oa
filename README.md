# project-oyster-oa

This repository focuses on efforts examining the impact of ocean acidification on Pacific oysters, *Crassostrea gigas*.

## Projects

Each project has a distinct ID used in all folder naming conventions to make sorting through the repository easier.

### Exploring Proteomic Variation in Pacific Oysters (ID: DNR)

In June 2016, 150 sibling *C. gigas* were outplanted at five different locations in Washington (leads: Micah Horowith (DNR) and Alex Lowe (UW Biology)). Four of these sites were located in Puget Sound — Case Inlet, Fidalgo Bay, Port Gamble Bay, Skokomish River Delta — and another in southern Washington, Willapa Bay. At each site, oysters were placed in both bare and eelgrass patches. After one month, gill tissue samples were collected and flash frozen. [Yaamini Venkataraman](yaaminiv.github.io) began proteomic analyses on samples to see if environmental variables like pH, dissolved oxygen, and temperature at each site and habitat drove differential protein expression.

### Pacific Oyster Ocean Acidification Trials at Chew Hatchery (ID: Manchester)

[Yaamini Venkataraman](yaaminiv.github.io) exposed adult *C. gigas* from Willapa Bay to either low and ambient pH conditions for seven weeks (February 4, 2017 to April 8, 2017). Before and after exposure, ctnidia, mantle and adductor tissues were collected, along with gonad tissue for histology. Live weights and shell lengths were also recorded.

After exposure, a subset of oysters were exposed to a one hour heat shock of 40ºC on June 5, 2017. Oysters were then conditioned for a strip spawn at 23ºC. Oysters were spawned on July 30, 2017 to create five groups based on parental life histories. Larvae survived for 20 days without growing past 60 microns. This trial will allow us to understand how parental life history affected larval production and differential mortality.

## Repository Contents

Each larger folder also has its own README.md file with additional information. With the exception of `analyses`, all folder contents are divided by project ID. Some highlights are outlined below.

[`analyses`](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses): Output for multiple analyses. Project IDs are detailed in the folder names themselves.

DNR:

- [BCA Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_BCA_Analysis)
- [Environmental Data Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Environmental_Data_Analysis_20161115)
- [MZ Ratios for Mass Spectrometry](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MZratios_Larval_Samples_20170118)
- [MSConvert](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_MSConvert_20170412)
- [Preliminary DIA Analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Preliminary_Analyses_20170321)
- [MSStats Analysis](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_Skyline_20170524)
- [SRM Transition Selection](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_TransitionSelection_20170707)
- [SRM Analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses/DNR_SRM_20170902)

Manchester: No analyses as of 2017-09-27

[`data`](https://github.com/RobertsLab/project-oyster-oa/tree/master/data): Raw data used for project analyses.

[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/DNR):

- [*C. gigas* Sample Key](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/PugetSound-DNR-samples-2016.csv)
- [Oysters used for biomarkers by Alex Lowe](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Biomarker_proteomeCrossRef.xlsx)
- [Environmental Summary Data](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Environmental%20Summary%20Data%20for%20Proteomics%20Project.csv)

[Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester):

- [Adult *C. gigas* Tissue Sampling Key](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae)
- [Pacific Oyster Larvae Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling)
- [Water Chemistry Information](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/Water-Chemistry-Data) 

[`images`](https://github.com/RobertsLab/project-oyster-oa/tree/master/images): Relevant photos for downstream analysis

[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR):

- [General Lab Notebook](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR/Lab-Notebook)
- [Skyline Cleaning Screenshots](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/DNR/2017-08-31-Skyline-Cleaning-Screenshots)

[Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester):

- [February Oyster Sampling](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/2017-2-4-Oyster-Sampling)
- [*C. gigas* Gonad Histology](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-gonad-histology)
- [*C. gigas* Larval Images](https://github.com/RobertsLab/project-oyster-oa/tree/master/images/Manchester/Gigas-larvae)

[`miscellaneous`](https://github.com/RobertsLab/project-oyster-oa/tree/master/miscellaneous): Meeting notes and other background information

[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/miscellaneous/DNR)

Manchester: No contents as of 2017-09-27

[`notebooks`](https://github.com/RobertsLab/project-oyster-oa/tree/master/notebooks): Jupyter notebook that detail reproducible methods for data analysis

[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/notebooks/DNR):

- [Preliminary Proteomic Analyses](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-03-21-Preliminary-Proteomic-Data-Analyses.ipynb)
- [Demultiplexing Proteomic Data](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-04-12-Demultiplex-Raw-Files.ipynb)
- [Shotgun Proteomic Analyses](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-06-13-Full-Skyline-Preliminary-Analysis.ipynb)
- [SRM Target Selection in MSStats](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-06-22-Selecting-SRM-Targets-with-MSstats-Part-2.ipynb)
- [Protein Annotations](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-07-05-Examining-Protein-Annotations.ipynb)
- [SRM Target ID in Skyline](https://github.com/RobertsLab/project-oyster-oa/blob/master/notebooks/DNR/2017-07-07-SRM-Target-Identification-in-Skyline.ipynb)

Manchester: No contents as of 2017-09-27

[`presentations`](https://github.com/RobertsLab/project-oyster-oa/tree/master/presentations): Final slides for various conferences

[DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/presentations/DNR):

- [UW SAFS Graduate Student Symposium 2016](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/DNR/GSS2016_Venkataraman.pdf)
- [Pacific Coast Shellfish Growers Association 2017](https://github.com/RobertsLab/project-oyster-oa/blob/master/presentations/DNR/PCSGA2017_Venkataraman.pptx)

Manchester: No presentations as of 2017-09-27