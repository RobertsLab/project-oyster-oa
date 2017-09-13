#Before going through this script, I went through "2017-09-10-NMDS-ANOSIM-for-Cluster-Analysis." In this script, I'll depict normalized protein area across samples at bar charts

#### BASELINE DATA MANIPULATION ####

SRMDataBarplots <- read.csv("2017-09-07-Master-SRM-Data-BiologicalReplicates-NoBlanks-NoPivot.csv", na.strings = "N/A") #Read in data
head(SRMDataBarplots)
SRMDataBarplots <- SRMDataBarplots[,-c(1, 3, 6, 8)] #Remove extraneous columns
SRMDataBarplots$Area <- as.numeric(SRMDataBarplots$Area)
SRMDataBarplots$TIC <- as.numeric(SRMDataBarplots$TIC)
SRMDataBarplots$Normalized.Area <- (SRMDataBarplots$Area)/(SRMDataBarplots$TIC) #Normalize area
SRMDataBarplots$Protein.Peptide.Transition <- paste(SRMDataBarplots$Protein.Name, SRMDataBarplots$Peptide.Sequence, SRMDataBarplots$Fragment.Ion) #Combine protein IDs
head(SRMDataBarplots)
SRMDataBarplots <- SRMDataBarplots[,-c(3:6)] #Remove extraneous columns, but keep protein name
head(SRMDataBarplots)
attach(SRMDataBarplots)
SRMDataBarplots <- SRMDataBarplots[order(Protein.Name),] #Reorder by protein name
head(SRMDataBarplots)
detach(SRMDataBarplots)

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
