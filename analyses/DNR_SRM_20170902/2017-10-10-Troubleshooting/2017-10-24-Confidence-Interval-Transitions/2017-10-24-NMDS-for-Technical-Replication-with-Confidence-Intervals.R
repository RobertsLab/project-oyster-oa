#In this script, I'll use the slope of my regression line to identify transitions that should be eliminated from analyses.

#### SET WORKING DIRECTORY ####
setwd(dir = "../..") #Set the working directory to the folder with all SRM analysis files (project-oyster-oa/analyses/DNR_SRM_20170902)
getwd()

#### IMPORT DATA ####

SRMDataTargetsOnlyPivotedCorrected <-  read.csv("2017-10-10-Troubleshooting/2017-10-10-Transition-Replicate-Correlations/2017-10-24-Targets-Replicates-Adjacent-Columns.csv", header = TRUE) #Import dataframe
rownames(SRMDataTargetsOnlyPivotedCorrected) <- SRMDataTargetsOnlyPivotedCorrected$X #Set rownames as first column
SRMDataTargetsOnlyPivotedCorrected <- SRMDataTargetsOnlyPivotedCorrected[,-1] #Remove column X
head(SRMDataTargetsOnlyPivotedCorrected) #Confirm changes

#### REFORMAT DATA ####
#Now I'll take my dataframe and split it into two: one for each batch of technical replicates.

SRMDataTargetsReplicateOne <- SRMDataTargetsOnlyPivotedCorrected[, c(seq(from = 1, to = (length(SRMDataTargetsOnlyPivotedCorrected) - 1), by = 2))] #Subset all odd columns (first replicate)
colnames(SRMDataTargetsReplicateOne) #Confirm subset
SRMDataTargetsReplicateTwo <- SRMDataTargetsOnlyPivotedCorrected[, c(seq(from = 2, to = length(SRMDataTargetsOnlyPivotedCorrected), by = 2))] #Subset all even columns (first replicate)
colnames(SRMDataTargetsReplicateTwo) #Confirm subset

#Finally, I'll transpose each dataframe. The resulting dataframes will have transitions in the columns and samples in the rows.

SRMDataTransposedReplicateOne <- t(SRMDataTargetsReplicateOne) #Transpose Replicate 1 dataframe
head(SRMDataTransposedReplicateOne) #Confirm transposition
SRMDataTransposedReplicateTwo <- t(SRMDataTargetsReplicateTwo) #Transpose Replicate 2 dataframe
head(SRMDataTransposedReplicateTwo) #Confirm transposition

#### CREATE A NULL CONFIDENCE INTERVAL MODEL ####
#I FOUND THIS DOESN'T WORK BECAUSE A CONFIDENCE INTERVAL AROUND A PERFECT MODEL DOESN'T EXIST.
#In order to make a 95% confidence interval around a line, I need a model to work with.

#confidenceIntervalData <- data.frame(x = c(1:7000),y = c(1:7000)) #Create a dataframe with data for a confidence interval around a x = y line
#confidenceIntervalModel <- lm(confidenceIntervalData$y ~ confidenceIntervalData$x) #Use x values to predict y values in model

#### CALCULATE CONFIDENCE INTERVAL ####
#I need to know the MSE and SSX to calculate my confidence interval. I can calculate these values outside of my loop function since these values are not going to change between transitions.

#mse <- (summary(confidenceIntervalModel)$sigma)^2 #Calculate MSE using residual standard error from summary statistics
#ssx <- sum((confidenceIntervalData$x - mean(confidenceIntervalData$x))^2) #Calculate the sum of squares for x values
#x <- seq(from = min(confidenceIntervalData$x), to = max(confidenceIntervalData$x), by = 1) #Create a sequence of all possible x values
#y.fitted <- summary(confidenceIntervalModel)$coeff[1] + summary(confidenceIntervalModel)$coeff[2]*x #Formula for prediction line
#y.upper <- y.fitted + sqrt(2*qf(0.95, 2, 44))*sqrt(mse*(1/12 + (x - mean(confidenceIntervalData$x))^2/ssx)) #Formula for upper confidence bound. Generate critical value from F distribution, where alpha = 0.95, numerator df = 2 (estimating intercept and slope), residual standard error df = 44 (46 samples - 2).
#y.lower <- y.fitted - sqrt(2*qf(0.95, 2, 44))*sqrt(mse*(1/12 + (x - mean(confidenceIntervalData$x))^2/ssx)) #Lower confidence bound

#### CREATE PLOTS FOR EACH TRANSITION ####
#I will regress my second batch of technical replicates against my first, plot the regression line, and the x = y line with its confidence interval.

setwd(dir = "2017-10-10-Troubleshooting/2017-10-24-Confidence-Interval-Transitions/") #Change working directory so files are saved in the same directory as the R script
getwd() #Confirm changes

correlationFilenames <- data.frame(filenames = colnames(SRMDataTransposedReplicateOne),
                                   modifier = rep("confint.jpeg", 111)) #Make a dataframe of filenames
correlationFilenames$full <- paste(correlationFilenames$filenames, correlationFilenames$modifier) #Merge the two columns together in a third column. This column has the full filename that will be used
head(correlationFilenames) #Confirm changes

nTransitions <- nrow(SRMDataTargetsReplicateOne) #Number of transitions used
for(i in 1:nTransitions) { #For all transitions
  transitionModel <- lm(SRMDataTransposedReplicateTwo[,i] ~ SRMDataTransposedReplicateOne[,i]) #Predict Replicate 2 from Replicate 1
  fileName <- correlationFilenames$full[i] #Set filename choice as the ith entry
  jpeg(filename = fileName, width = 1000, height = 1000) #Save .jpeg using set filename
  plot(x= SRMDataTransposedReplicateOne[,i], y = SRMDataTransposedReplicateTwo[,i], xlab = "Replicate 1 Area", ylab = "Replicate 2 Area", main = correlationFilenames$filename[i], type = "n") #Create plot, but do not plot points
  text(x = SRMDataTransposedReplicateOne[,i], y = SRMDataTransposedReplicateTwo[,i], labels = rownames(SRMDataTransposedReplicateOne), cex = 0.7) #Plot sample ID instead of points
  abline(transitionModel, col = "red") #Plot regression
  legend("topleft", bty = "n", legend = paste("R2 =", format(summary(transitionModel)$adj.r.squared, digits=4))) #Plot R-squared value
  abline(a = 0, b  = 1, col = "black", lty = 3)
  legend("bottomright", bty = "n", legend = paste("regression line =", format(summary(transitionModel)$coeff[1], digits = 4), "+", format(summary(transitionModel)$coeff[2], digits = 4), "*x"))
  x.slope <- seq(min(SRMDataTransposedReplicateOne[,i]), max(SRMDataTransposedReplicateOne[,i]), by = 1) #Create a sequence of x values
  y.slope <- summary(transitionModel)$coeff[1] + 1*x.slope #Formula for prediction line, taking the intercept from the regression line and setting the slope to 1
  lines(lowess(x.slope, y.slope), col = "blue")
  dev.off() #Turn off plotting mechanism
}