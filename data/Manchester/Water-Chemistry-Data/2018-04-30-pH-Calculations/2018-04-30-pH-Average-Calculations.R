#I will calculate average and standard error for pH from the discrete samples.

#### IMPORT TOTAL ALKALINITY VALUES ####
pHTanks <- read.csv("2018-04-30-pH-Discrete-Samples-by-Tank.csv", header = TRUE) #Import pH values
head(pHTanks) #Confirm import

#### AVERAGE VALUES ####
samplingDays <- c(2, 2, 5, 5, 7, 7, 12, 12, 14, 14, 19, 19, 21, 21, 26, 26, 28, 28, 33, 33, 36, 36, 42, 42, 44, 44, 48, 48, 52, 52) #Isolate sampling days
averagepH <- data.frame("Day" =  samplingDays,
                        "Treatment" = rep(c("Experiment", "Control"), times = ((length(samplingDays))/2)),
                        "averagepH" = rep(0, times = length(samplingDays)),
                        "standardError" = rep(0, times = length(samplingDays))) #Create an empty dataframe
head(averagepH) #Confirm dataframe creation
npH <- as.numeric(length(averagepH$averagepH)) #Calculate length

for(i in 1:npH){
  averagepH$averagepH[i] <- mean(pHTanks$pH[((3*i)-2):(3*i)])
} #Calculate means and add them to the table
head(averagepH) #Confirm additions

#### CALCULATE STANDARD ERROR ####
for(i in 1:npH){
  averagepH$standardError[i] <- sqrt(var(pHTanks$pH[((3*i)-2):(3*i)]))
} #Calculate standard errors and add them to the table
head(averagepH) #Confirm additions

#### EXPORT DATA ####
write.csv(averagepH, "2018-04-30-Average-pH.csv") #Export dataframe
