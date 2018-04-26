#I will calculate average and standard error values for the total alkaninity measurements Sam gave me on 4/25/2018.

#### IMPORT TOTAL ALKALINITY VALUES ####
totalAlkalinity <- read.csv("2018-04-26-Total-Alkalinity-per-Tank.csv", header = TRUE) #Import TA values
head(totalAlkalinity) #Confirm import

#### AVERAGE VALUES ####
averageAlkalinity <- data.frame("Treatment" = rep(c("Experiment", "Control"), times = 3),
                                "Date" = c("2/20/17", "2/20/17", "3/20/17", "3/20/17", "4/4/17", "4/4/17"),
                                "averageAlkalinity" = rep(0, times = 6),
                                "standardError" = rep(0, times = 6)) #Create an empty dataframe
head(averageAlkalinity) #Confirm dataframe creation
nAlkalinity <- as.numeric(length(averageAlkalinity$averageAlkalinity)) #Calculate length

for(i in 1:nAlkalinity){
  averageAlkalinity$averageAlkalinity[i] <- mean(totalAlkalinity$totalAlkalinity[((3*i)-2):(3*i)])
} #Calculate means and add them to the table
head(averageAlkalinity) #Confirm additions

#### CALCULATE STANDARD ERROR ####
for(i in 1:nAlkalinity){
  averageAlkalinity$standardError[i] <- sqrt(var(totalAlkalinity$totalAlkalinity[((3*i)-2):(3*i)]))
} #Calculate standard errors and add them to the table
head(averageAlkalinity) #Confirm additions

#### EXPORT DATA ####
write.csv(averageAlkalinity, "2018-04-26-Average-Total-Alkalinity.csv") #Export dataframe
