# `data` Subdirectory Structure

Raw data used for [analyses](https://github.com/RobertsLab/project-oyster-oa/tree/master/analyses), organized by project ID. Repository contents and relevant metadata below.

## [DNR](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/DNR)

### **[PugetSound-DNR-samples-2016.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/PugetSound-DNR-samples-2016.csv)**: 
Guide to samples outplanted in the summer of 2016

- Round: Exposure period 1 (June to July 2016), exposure period 2 (July to August 2016)       
- Site: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)      
- Spp: Pacific oyster (Cg)
- Patch: Bare (B), Eelgrass (E)       
- Exclosure: There were 3 predator exclosures in each patch (mesh cages for oysters)      
- Rep: There were 5 juveniles in each exclosure. Due to mortality, we didn't always collect 5 individuals at the end of the exposure period     
- PRVial: Label for vial containing dissected ctenidia

### **[Biomarker_proteomeCrossRef.xlsx](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Biomarker_proteomeCrossRef.xlsx)**: 
From Alex Lowe (UW Biology), indicating which samples he used for his fatty acid analysis

- Sample: Randomly assigned number- Species: Pacific (*C. gigas*)- Round: Exposure period 1 (June to July 2016), exposure period 2 (July to August 2016)       
- Site: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)      
- Habitat: Bare (B), Eelgrass (E) - Replicate: Number, based on same site and eelgrass conditions- Vial: Label for vial containing fatty acid- PRVial: Label for vial containing dissected ctenidia### **[2017-11-14-Environmental-Data-from-Micah.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-14-Environmental-Data-from-Micah.csv)**:Temperature, dissolved oxygen, and pH data from probes processed by Micah Horwith (DNR)- DateTime: Date and military time when reading was taken- Date: Date portion of DateTime- Time: Time portion of DateTime
The rest of the column headers are as follows: SiteParameter (Probe ID)- Site: Combination of Location and habitat- Location: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)  - Habitat: Bare (B), Eelgrass (E)
- Parameter: pH, pHT (temperature from pH logger), do (dissolved oxygen), doT (temperature from dissolved oxygen logger)
- Probe ID: irrelevant  
### **[2017-11-25-Calculated-Salinity-Output-from-Micah.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-25-Calculated-Salinity-Output-from-Micah.csv)**:
Salinity calculated by Micah Horwith (DNR)

- Date: Day-Month-Year
- Time: Military time when reading was taken

The rest of the column headers are as follows: SiteParameter- Site: Combination of Location and habitat- Location: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)  - Habitat: Bare (B), Eelgrass (E)
- Parameter: ctS (salinity), ctT (temperature from salinity logger)

### **[2017-12-13-Tidal-Data-by-Site.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-12-13-Tidal-Data-by-Site.csv)**:
Tidal data for each outplant location

- Date: Month-Day-Year
- Time: Military time

The rest of the column headers are as follows: Location-Tide

- Location: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)

### **[Environmental Summary Data for Proteomics Project.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/Environmental%20Summary%20Data%20for%20Proteomics%20Project.csv)**: 
Average temperature and dissolved oxygen at each site and habitat over the course of the one month outplant

- Outplant #: Exposure period 1 (June to July 2016), exposure period 2 (July to August 2016)
- Site: Case Inlet (CI), Fidalgo Bay (FB), Port Gamble Bay (PG), Skokomish River Delta (SK), Willapa Bay (WB)      
- Habitat: Bare (B), Eelgrass (E) - Mean Temp: ºC, from continuous monitoring
- Mean [DO]: From continuous monitoring

### **[2017-11-21-Alex-Data-Yaamini-Samples-Only.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-11-21-Alex-Data-Yaamini-Samples-Only.csv)**:
Biomarker data from Alex Low (UW Biology) for oysters used in proteomic analyses

- Sample: Alex's sample ID
- X
- PRVial: Yaamini's sample ID
- Round: Outplant 1 or 2
- Site.x: Location for outplant
- Habitat: Bare or Eelgrass
- Spp
- Bag
- Rep.x
- FAvial: Alex's sample ID
- FinalHeight: Distance between umbow to top of shell
- topValve
- bothValves
- TissueMass: Wet weight
- shellThick: Shell thickness
- peakLoad: Raw strength measure
- Strength: Calculated shell strength
- delc
- C.perc
- delN
- N.perc
- CN.ratio: C.perc/N.perc

### **[2017-12-05-OysterGrowth.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/DNR/2017-12-05-OysterGrowth.csv)**: 
*C. gigas* growth data from Micah Horwith (DNR)

- Round: Outplant 1 or 2
- Site: Location for outplant
- Habitat: Bare or Eelgrass
- Exclosure
- Rep
- PRVial: Yaamini's sample ID
- Fshell: Final shell length
- Ishell1 through Ishell5: Initial shell lengths of five random oysters in same exclosure and replicate
- AvgIshell: Average initial shell length for exclosure and replicate

## [Manchester](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester)

### **[Water-Chemistry-Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/Water-Chemistry-Data)**: 
HOBO data from *C. gigas* tanks

### **[2017-Adult-Gigas-Tissue-Sampling](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling)**: 
Tissue sampling information [before](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling/20170204-GigasTissueSamplingInformation.csv) and [after](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling/20170408-GigasTissueSamplingInformation.xlsx) exposure to differing pH conditions

#### **[2017-02-04 Sampling](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling/20170204-GigasTissueSamplingInformation.csv)**:

- Oyster number: Used for animal ID
- Adductuor Tissue Tube ID: A-OysterNumber, tube holds adductor tissue for specific animal
- Ctenida (Gill Tissue) Tube ID: C-OysterNumber, tube holds ctenida tissue for specific animal
- Mantle Tissue Tube ID: M-OysterNumber, tube holds mantle tissue for specific animal
- Ethanol Tube ID: E-OysterNumber, tube holds ctenidia tissue for specific animal in ethanol for DNA analyses
- Gonad Histology Cassette ID: Location of gonad tissue for animal. Cassettes distinguished by number of notches on the side of the cassette (ex. 1 notch, 4 notches)
- Location in Cassette: Cassettes have four quadrants. top-left, bottom-left, top-right, bottom-right
- Length: cm, measured using calipers
- Shell and animal weight: g, measured using scale
- Empty shell weight: g, measured using scale
- Calculated animal weight: g, Empty shell weight subtracted from shell and animal weight 
- Additional notes

#### **[2017-04-08 Sampling](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-Adult-Gigas-Tissue-Sampling/20170408-GigasTissueSamplingInformation.xlsx)**:

- Number: Used for oyster-ID
- oyster-ID: Number and tank the oyster came from. Six tanks were used, three with low pH conditions (T1, T2 and T3) and three ambient pH conditions (T4, T5 and T6). (ex. Oyster 6 from Tank 1 is 6-T1)
- length: cm, measured using calipers
- shell-and-animal-weight: g, measured using scale
- empty-shell-weight: g, measured using scale
- calculated-biomass: g, empty-shell-weight subtracted from shell-and-animal-weight 
- adductor-ID: oysterID-A, tube holds adductor tissue for specific animal
- ctenida-ID: oysterID-C, tube holds ctenida tissue for specific animal
- mantle-ID: oysterID-M, tube holds mantle tissue for specific animal
- ethanol-ID: oysterID-E, tube holds ctenidia tissue for specific animal in ethanol for DNA analyses
- histology-cassette: Location of gonad tissue for animal. Cassettes distinguished by species and a number ID (ex. gigas-1)
- cassette-position: Cassettes have four quadrants. top-left, bottom-left, top-right, bottom-right
- notes

### **[2017-07-30-Pacific-Oyster-Larvae](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae)**: 
Data collected for larval experiment during August 2017

#### **[HOBO Data](https://github.com/RobertsLab/project-oyster-oa/tree/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/HOBO-Data)**: 
Temperature data from HOBO loggers in larval buckets

#### **[2017-07-29-Spawning-Calculations.xlsx](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-17-Live-Larvae-Counts.jpg)**: 
Information used to generate crosses used for *C. gigas* spawning on 2017-07-30

**Sex**:

- Oyster-ID: Tank from original pH exposure (low pH; 1-3, ambient pH: 4-6) (ex. 3-4)
- Sex: Male (M) and Female (F)

**Pools**:

- Oyster-ID: Tank from original pH exposure (low pH; 1-3, ambient pH: 4-6) (ex. 3-4)
- Weight Contributed: g of gametes used, dry weight
Notes
- Hydrate Start: When 23ºC was added to gametes
- Hydrate End: When gametes were mixed together, 45 minute hydration goal

**Crosses**:

- Bucket number: 1-60
- Female Pool: Life history of female used for cross. Low, Ambient or Heat Shock
- Sample Egg Density: Number of eggs counted 250 µL samples
- Female Egg Count: Total eggs in a female pool
- Container Volume: mL, bucket volume
- Stocking Density: eggs/mL, goal per bucket
- Amount Low Female Added: mL, volume of pool used for spawn
- Ambient Male ID: Male used for cross
- Amount Male Added: mL, volume of sperm used for spawn
- Fertilization: Yes or No
- Location: Bucket placement, 1 (Outside) or 0 (Inside)
- Notes

**Hatching**:

- Tripour Number: Bucket IDs for two replicates with identical cross (ex. Bucket 1 and 15 were the same cross, so the Tripour Number is 1-15)
- Amount Female Added: mL, volume of pool used for spawn
- Number of Eggs Added: Eggs used for spawn
- D-Hinge Count: Number of d-hinge counted
- Volume: mL, Volume of sample used for d-hinge count
- Average D-Hinge Count: Average of three counts
- Tripour Volume: µL, total volume of tripour with larvae
- D-Hinge in Tripour: Number of d-hinge larvae produced by each cross
- Average Hatch Rate: D-Hinge in tripour from total number of eggs used

**Restocking**:

- Treatment Cross: Life history of male and female used in cross
- Tripours Numbers: Tripour IDs with parental life histories designated in Treatment Cross
- Color: Tape color used to identify life histories
- Total D-Hinge: Total d-hinge produced for each treatment cross
- Larvae Needed per Bucket: Target number of larvae for experiment
- Tripour Volume: mL, volume of tripour with larvae
- Volume Needed Per Bucket: mL, amount of tripour needed for target larvae
- Final Bucket Numbers: Bucket IDs for treatment cross
- Maximum Total Larvae Per Replicate: Highest number of larvae that could be added to five replicates

**Algae**:

- strains: Algae used
- count: from hemocytometer
- average: average algal cell count
- Cells/mL: number of cells in hemocytometer
- Cells Needed: target number of cells to feed larvae
- Volume Needed: mL, target volume of strains per bucket
- Volume Added: mL, volume algae added

#### **[2017-07-30-Feeding.xlsx](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-07-30-Feeding.xlsx)**:

**Pre-Feeding Concentrations**: 
Concentration of algae in buckets prior to feeding

- Date
- Bucket: Bucket ID from [`2017-07-29-Spawning-Calculations.xlsx`](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-17-Live-Larvae-Counts.jpg)
- Count: Algal cells counted in hemocytometer
- Average: Average algal cell count
- Cells/mL: Cells in hemocytometer
- Cells Present: Cells in bucket
- Cells Needed: Target algal cell count for feeding
- Cells To Be Added: Difference between cells needed and present
- Notes

**Algae Added**: 

- Date
- Strains: Algae used
- Count: Algal cells counted in hemocytometer
- Average Cell Count
- Cells/mL: Cells in hemocytometer
- Cells Needed: Target algal cell count for feeding
- Volume Needed: mL, Target volume of strains per bucket
- Volume Added: mL, Volume algae added
- Cells Added: Number of algal cells added to bucket based on Volume Added
- Algal Concentration: Cells Added/Volume Added (mL)
- Notes

#### **[2017-07-31-Temperature.xlsx](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-07-31-Temperature.xlsx)**:

- Date
- Bucket: Bucket ID from [`2017-07-29-Spawning-Calculations.xlsx`](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-17-Live-Larvae-Counts.jpg)
- Temperature: Upon arrival, used temperature gun
- Tote: Tote number where bucket is located
- Notes

#### **[2017-08-02-Larvae-Counts.xlsx](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-02-Larvae-Counts.xlsx)**:

- Date Sampled
- Date Counted
- Bucket: Bucket ID from [`2017-07-29-Spawning-Calculations.xlsx`](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-17-Live-Larvae-Counts.jpg)
- Screen Size: microns, Mesh size used to catch larvae
- Tripour Volume: mL, Volume of container larvae were in
- Sample Volume: µL, Volume sampled for counts
- Count Alive
- Count Dead
- Number Alive: Total in bucket
- Total Live Larvae in Treatment: Number of larvae alive in each parental life history group
- Larvae Needed: mL, used to sample 10,000 larvae per treatment group for epigenetic analyses
- Number Sampled per Treatment: Larvae sampled for epigenetic analyses
- Remaining Live Larvae Per Treatment: Larvae left after sampling
- Notes

#### **[2017-08-17-Live-Larvae-Counts.jpg](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/2017-07-30-Pacific-Oyster-Larvae/2017-08-17-Live-Larvae-Counts.jpg)**:
 
Number of live larvae remaining based on parental life history over the duration of the experiment (August 2017)

#### **[2018-04-26-Total-Alkalinity-per-Tank.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/Water-Chemistry-Data/2018-04-26-Total-Alkalinity-per-Tank.csv)**:

Information taken directly from [Sam's notebook entry](http://onsnetwork.org/kubu4/2018/04/24/total-alkalinity-calculations-yaaminis-ocean-chemistry-samples/).

- Tank: Experimental tank number from Manchester set-up
- Date: Discrete sampling date
- Total Alkalinity: µmol/kg, calculated by Sam using T5 Excellence Titrator (Metter Toledo)

#### **[2018-04-26-Average-Total-Alkalinity.csv](https://github.com/RobertsLab/project-oyster-oa/blob/master/data/Manchester/Water-Chemistry-Data/2018-04-26-Average-Total-Alkalinity.csv)**:

- Treatment: Experimental (low pH) or Control (ambient pH)
- Date: Discrete sampling date
- averageAlkalinity: µmol/kg, averages from three replicates
- standardError: µmol/kg, calculated from three replicates