#Before going through this script, I went through "2017-09-10-NMDS-ANOSIM-for-Cluster-Analysis." In this script, I'll depict normalized protein area across samples as bar charts

#### DATA MANIPULATION ####
SRMDataNMDSAveragedCorrected #From 2017-09-06-NMDS-for-Technical-Replication. Average normalized area data.
boxplotData <- data.frame(t(SRMDataNMDSAveragedCorrected)) #Transpose the data
boxplotData$Sample.Number <- rownames(boxplotData) #Save rownames as a new column
head(boxplotData) #Confirm changes

biologicalReplicates <- read.csv("2017-09-06-Biological-Replicate-Information.csv", na.strings = "N/A") #Import site and eelgrass condition information (i.e. biological replicate information)
head(biologicalReplicates) #Confirm import
rownames(biologicalReplicates) <- biologicalReplicates$Sample.Number #Set sample number as row names
head(biologicalReplicates) #Confirm changes
biologicalReplicates <- biologicalReplicates[-c(50,100),] #Remove blanks
biologicalReplicates$Site <- factor(biologicalReplicates$Site) #Remove 0 as a factor
biologicalReplicates$Eelgrass.Condition <- factor(biologicalReplicates$Eelgrass.Condition) #Remove 0 as a factor
str(biologicalReplicates) #Confirm factor reset

boxplotData <- merge(x = biologicalReplicates, y = boxplotData, by = "Sample.Number") #Merge together
head(boxplotData) #Confirm merge
rownames(boxplotData) <- boxplotData$Sample.Number #Set sample number as row names
boxplotData <- boxplotData[-1] #Remove Sample.Number column
head(boxplotData) #Confirm changes

#### MAKE BOXPLOTS JUST BASED ON SITES ####

nTransitions <- (length(boxplotData)) #The number of columns in the dataframe. The first 2 columns are Site and Eelgrass.Condition
boxplotFilenames <- colnames(boxplotData)
for(i in 3:nTransitions) { #For all of my columns with transition IDs
  fileName <- boxplotFilenames[i] #Set the file name
  jpeg(filename = fileName, width = 1000, height = 1000) #Save using set file name
  boxplot(boxplotData[,i] ~ boxplotData$Site, xlab = "Sites", ylab = "Abundance") #Create the boxplot
  dev.off() #Close file
}

boxplot(boxplotData[,3] ~ boxplotData$Site) #Make b


#### SUBSET DATA BASED ON SITES ####

CaseInlet <- subset(x = boxplotData, subset = boxplotData$Site == "CI")
FidalgoBay <- subset(x = boxplotData, subset = boxplotData$Site == "FB")
PortGamble <- subset(x = boxplotData, subset = boxplotData$Site == "PG")
SkokomishRiver <- subset(x = boxplotData, subset = boxplotData$Site == "SK")
WillapaBay <- subset(x = boxplotData, subset = boxplotData$Site == "WB")


#### BAR CHARTS WITHOUT AVERAGING ACROSS SITES AND CONDITIONS ####

#### BAR CHARTS AFTER AVERAGING ACROSS UNIQUE SITE + CONDITION COMBINATIONS ####

#### BAR CHARTS AFTER AVERAGING BY SITE ONLY ####

#Subset relevant data
CaseInletBarplot <- subset(x = SRMDataBarplots, subset = SRMDataBarplots$Site == "CI")
FidalgoBayBarplot <- subset(x = SRMDataBarplots, subset = SRMDataBarplots$Site == "FB")
PortGambleBarplot <- subset(x = SRMDataBarplots, subset = SRMDataBarplots$Site == "PG")
SkokomishRiverBarplot <- subset(x = SRMDataBarplots, subset = SRMDataBarplots$Site == "SK")
WillapaBayBarplot <- subset(x = SRMDataBarplots, subset = SRMDataBarplots$Site == "WB")

#Sort by protein name
attach(CaseInletBarplot)
CaseInletBarplot <- CaseInletBarplot[order(Protein.Peptide.Transition),] #Reorder by protein name
detach(CaseInletBarplot)
head(CaseInletBarplot) #Confirm reordering

#Average all values for one Protein.Peptide.Transition
mean(CaseInletBarplot$Normalized.Area[1:92])


#ONLY FOR PRACTICE TALK WILL MAKE THIS BETTER WHEN I GET MORE SLEEP.

caseCatalase <- subset(CaseInletBarplot, CaseInletBarplot$Protein.Name == "CHOYP_CATA.1.3|m.11120") #Subset catalase
fidalgoCatalase <- subset(FidalgoBayBarplot, FidalgoBayBarplot$Protein.Name == "CHOYP_CATA.1.3|m.11120")
portCatalase <- subset(PortGambleBarplot, PortGambleBarplot$Protein.Name == "CHOYP_CATA.1.3|m.11120")
skokomishCatalase <- subset(SkokomishRiverBarplot, SkokomishRiverBarplot$Protein.Name == "CHOYP_CATA.1.3|m.11120")
willapaCatalase <- subset(WillapaBayBarplot, WillapaBayBarplot$Protein.Name == "CHOYP_CATA.1.3|m.11120")

caseCatalaseP1T1 <- mean(caseCatalase$Normalized.Area[1:18]) #Average catalase values for the first peptide and first transition
fidalgoCatalaseP1T1 <- mean(fidalgoCatalase$Normalized.Area[1:18])
portCatalaseP1T1 <- mean(portCatalase$Normalized.Area[1:18])
skokomishCatalaseP1T1 <- mean(skokomishCatalase$Normalized.Area[1:18])
willapaCatalaseP1T1 <- mean(willapaCatalase$Normalized.Area[1:18])

catalaseP1T1 <- data.frame(Site = c("Case Inlet", "Fidalgo Bay", "Port Gamble Bay", "Skokomish River", "Willapa Bay"),
                           Area = c(caseCatalaseP1T1, fidalgoCatalaseP1T1, portCatalaseP1T1, skokomishCatalaseP1T1, willapaCatalaseP1T1))

barplot(catalaseP1T1$Area, xlab = "Sites", ylab = "Abundance", main = "Catalase Abundance across Sites", col = c("red", "blue", "black", "green", "magenta"), legend = TRUE)
legend("topleft", cex = .5, pch = c(rep(x = 16, times = 5)), legend=c('Case Inlet', "Fidalgo Bay", "Port Gamble", "Skokomish River", "Willapa Bay"), col=c('red', 'blue', 'black', 'green', 'magenta'))
