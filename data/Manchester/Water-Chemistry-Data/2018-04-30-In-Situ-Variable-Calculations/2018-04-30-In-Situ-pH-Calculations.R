#In this script, I will calculate in situ carbonate chemistry parameters for the three days with total alkalinity measurements.

#### SET WORKING DIRECTORY ####
setwd(getwd()) #Working directory should be the "2018-04-30-In-Situ-Variable-Calculations" directory

#### LOAD DEPENDENCIES ####
library(seacarb)

#### IMPORT DATAFRAME ####
seacarbInputs <- 

  
  
  
  
#Calculate CO2 parameters using seacarb
carb.ouptput <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$Corrected.TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb
carb.ouptput$ALK <- carb.ouptput$ALK*1000000 #convert to µmol kg-1
carb.ouptput$CO2 <- carb.ouptput$CO2*1000000 #convert to µmol kg-1
carb.ouptput$HCO3 <- carb.ouptput$HCO3*1000000 #convert to µmol kg-1
carb.ouptput$CO3 <- carb.ouptput$CO3*1000000 #convert to µmol kg-1
carb.ouptput$DIC <- carb.ouptput$DIC*1000000 #convert to µmol kg-1
carb.ouptput <- carb.ouptput[,-c(1,4,5,8,10:13,19)] #subset variables of interest
carb.ouptput <- cbind(SW.chem$Date,  SW.chem$Tank,  SW.chem$Treatment, SW.chem$Period1,SW.chem$Period2, SW.chem$Period3, carb.ouptput) #combine the sample information with the seacarb output
colnames(carb.ouptput) <- c("Date",  "Tank",  "Treatment",	"Period1", "Period2", "Period3",	"Salinity",	"Temperature", "pH",	"CO2",	"pCO2","HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat") #Rename columns to describe contents
write.table(carb.ouptput, "~/MyProjects/HI_Pdam_Parental/RAnalysis/Output/Seawater_chemistry_table_Output_All.csv", sep=",", row.names = FALSE) #save data

#### EXPORT DATA ####
#write.csv(averageAlkalinity, "2018-04-26-Average-Total-Alkalinity.csv") #Export dataframe
