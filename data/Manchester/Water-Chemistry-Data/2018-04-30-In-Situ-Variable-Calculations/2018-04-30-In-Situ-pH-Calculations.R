#In this script, I will calculate in situ carbonate chemistry parameters for the three days with total alkalinity measurements.

#### SET WORKING DIRECTORY ####
setwd(getwd()) #Working directory should be the "2018-04-30-In-Situ-Variable-Calculations" directory

#### LOAD DEPENDENCIES ####
library(seacarb)

#### IMPORT DATAFRAME ####
seacarbInputs <- read.csv("2018-04-30-Seacarb-Inputs.csv", header = TRUE)
head(seacarbInputs) #Confirm import

#### CALCULATE PARAMETERS USING SEACARB ####
carb.output <- carb(flag = 8, var1 = seacarbInputs$pH, var2 = seacarbInputs$TotalAlkalinity/1000000, S = seacarbInputs$Salinity, T = seacarbInputs$Temperature, P = 0, Pt = 0, Sit = 0, pHscale = "T", kf = "pf", k1k2 = "1", ks = "d") #Calculate seawater chemistry paramters
carb.ouptput$ALK <- carb.ouptput$ALK*1000000 #convert to µmol kg-1
carb.ouptput$CO2 <- carb.ouptput$CO2*1000000 #convert to µmol kg-1
carb.ouptput$HCO3 <- carb.ouptput$HCO3*1000000 #convert to µmol kg-1
carb.ouptput$CO3 <- carb.ouptput$CO3*1000000 #convert to µmol kg-1
carb.ouptput$DIC <- carb.ouptput$DIC*1000000 #convert to µmol kg-1
carb.ouptput <- carb.ouptput[,-c(1,4,5,8,10:13,19)] #subset variables of interest

#### EXPORT DATA ####
#write.csv(carb.output, "2018-04-30-In-Situ-Carbonate-Chemistry-Parameters.csv", rownames = FALSE)