#Step 1: Randomly choose 10 peptide files to download and create a histogram from

peptideFiles <- c("11A", "11", "12A", "12", "16A", "16", "19A", "19", "1A", "1", "20A", "20", "24A", "24", "27A", "27", "28A", "28", "32A", "32", "35A", "35", "36A", "36", "3A", "3", "40A", "40", "43A", "43", "44A", "44", "48A", "48", "4", "51A", "51", "52A", "52", "56A", "56", "8A", "8") #Create a set with peptide sample numbers to randomly draw from

set.seed(124) #Set seed to ensure reproducibility
sample(peptideFiles, 10) #Pick 10 sample peptide files. The files chosen were: "12"  "28"  "35"  "27"  "1A"  "20"  "56"  "8A"  "48A" "1". However, because I can't have both 1 and 1A, I chose the following files instead: 1A, 8A, 12, 20, 27, 28, 32, 35, 48A, and 56.

#Step 2: Import datasets

sample1A <- read.table(file = "interact-20161205_sample_1A.pep.xls", header = TRUE)
sample8A <- read.table(file = "interact-20161205_sample_8A.pep.xls", header = TRUE)
sample12 <- read.table(file = "interact-20161205_sample_12.pep.xls", header = TRUE)
sample20 <- read.table(file = "interact-20161205_sample_20.pep.xls", header = TRUE)
sample27 <- read.table(file = "interact-20161205_sample_27.pep.xls", header = TRUE)
sample28 <- read.table(file = "interact-20161205_sample_28.pep.xls", header = TRUE)
sample32 <- read.table(file = "interact-20161205_sample_32.pep.xls", header = TRUE)
sample35 <- read.table(file = "interact-20161205_sample_35.pep.xls", header = TRUE)
sample48A <- read.table(file = "interact-20161205_sample_48A.pep.xls", header = TRUE)
sample56 <- read.table(file = "interact-20161205_sample_56.pep.xls", header = TRUE)

#Step 3: Use rbind to merge all datasets together

allsamples <- rbind(sample1A, sample8A, sample12, sample20, sample27, sample28, sample32, sample35, sample48A, sample56)

#Step 4: Create a histogram using frequency of MZ ratios and export as a png

png("MZratios-larval-samples-histogram", width = 800, height = 800)
hist(allsamples$MZratio, freq = TRUE, main= "M/Z ratio frequency for larval Crassostrea gigas samples", col="darkturquoise", xlab= "M/Z ratios", ylab = "Frequency")
dev.off()
