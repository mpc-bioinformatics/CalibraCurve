#################################
### CalibraCurve Version 2.0
#################################

#################################
### Settings
#################################

# Notes: 
# 1. All settings that must be adapted by the user can be given in the chunk 'User settings'.
# 2. The chunk 'Progam settings' comprises settings that are related to data processing and software execution.
# 3. The chunk 'Plot settings' can be used to adapt the look of the CalibraCurve plots. 
# 4. The settings in the chunk named 'Further settings' concern fine tuning and customization of CalibraCurve 
# (e.g. adaption of the diagrams and logging).

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# User settings (MUST BE ADAPTED BY THE USER)
# Notes:
# Absolute path specifications are required. 
# Please use the slash ("/") or the double backslash character(s) ("\\") as 
# delimiting characters for the path components/directories.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

inputDir <- "PATH/TO/YOUR/INPUT/DIRECTORY"
outputDir <- "PATH/TO/YOUR/OUTPUT/DIRECTORY"

# The value for 'fileTypeOfData' specifies the type of the files that contain the data.
# Please note that no other files of the same file type are allowed in the input directory.
fileTypeOfData <- "csv"
# Column delimiter used in the data files
colDelimiterIn = ';'
# Column delimiter used for the result files
colDelimiterOut = '\t'

# The number of the column that contains the concentration values
colNumberConcentration <- 6
# The number of the column that contains the measurement values
colNumberMeasurements <- 7

# Unit that is used to for the concentration quantities
concUnit <- "[fmol/µl]"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Progam settings
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Settings related to quality control 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Minimum number of replications required for each concentration level (value must be greater than 1)
minNumberReplications <- 3

# Settings related to the assessment of precision (used for calculation of the preliminary linear range)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Threshold for the CV [%] that is calculated for each concentration level  
cv_thres <- 20

# Settings related to the assessment of (percentage) bias (used for the calculation of the final linear range)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Thresholds
# Threshold for the maximal allowed percent bias 
perBiasThres <- 20
# Method used for calculation of weights:
# You may apply your own weighting method. To this end change the value of the parameter 'weightingMethod' to a character string 
# of your choice (e.g. 'OwnWeights'). Please note that the values '1/x^2' and '1/x' refer to the weighting methods pre-implemented in CalibraCurve.
# For a detailed description of the complete working steps required for utilization of user defined weights, please refer to appendix 2 of the 
# CalibraCurve manual. 
weightingMethod <- '1/x^2'
# The value for 'centralTendencyMeasure' specifies whether mean or median values should be calculated as a measure for central tendency
centralTendencyMeasure <- 'mean'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Advanced settings
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Advanced Program settings 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Switch for the calculation method for the preliminary linear range:
# TRUE:   All CV values within the preliminary range must pass the CV threshold criteria.
# FALSE:  The preliminary range extends from the lowest to the highest concentration level that pass the CV threshold criteria.
#         Intermediate concentration levelS may fail to pass the threshold criteria.
calcContinuousPrelimRanges <- TRUE
# Switches consideration of CV values calculated from the percent bias for the selection of a 
# concentration level on and off.
considerPerBiasCV<- TRUE
# Threshold for the distance of the average percent bias values calculated for both the currently highest and lowest concentration level. 
# This threshold is used within the algorithm that selects a concentration level for removal.
perBiasDistThres <- 10 

# Settings related to the linear models
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The value for 'finalRangeCalculationMethod' specifies whether the weighted or the unweighted model should be used for the calculation of the final linear range
finalRangeCalculationMethod <- 'weighted_linear_model' 


# Advanced settings calibration graph(s)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The value for 'plotOptReg' specifies whether unweighted and weighted regression lines are drawn in the same plot.
plotOptReg <- 'combined' 
# The value for 'scaling' specifies whether logarithmic or linear scale is used for the plots.
scaling <- 'log'
# Factor that specifies the magnification of the axes annotations relative to the default (represented by the value 1)
magnificAxis = 1
# Factor that specifies the magnification of the axes labels relative to the default (represented by the value 1)
magnificLabels = 1 
# Factor that specifies the magnification of the title relative to the default (represented by the value 1)
magnificTitle = 1
# Factor that specifies the magnification of the equations in the plot relative to the default (represented by the value 1)
magnificEquation = 1
# Resolution of the plots given in dpi
plotResolution <- 300
# Height of the plot given in cm
plotHeight <- 17
# width of the plot given in cm
plotWidth <- 17

# Specific texts for the regression plot
titleOpening <- "Calibration graph for "
titleMiddlePart <-"Linear range extends from "
xLab <- paste("Log concentration ", concUnit, sep = "")
yLab <- "Log area heavy"
legend_unweighted <-"blue: Regression line, equation given for the unweighted fit" 
legend_weighted <- "red: Regression line, equation given for the weighted fit"
# Specifies whether or not a legend is shown in the CalibraCurve plots
showLegend <- TRUE

# Advanced settings response factor graphs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Enables adaption of the y-axis range with user defined settings
yLimAdaptRF <- FALSE
# Lower limit of the y-axis (needs adaption depending on analyzed data!!! if yLimAdaptRF == TRUE)
# minY_lim <- -0.05
minY_lim <- 100000
# Upper limit of the y-axis (needs adaption depending on analyzed data!!! if yLimAdaptRF == TRUE)
# maxY_lim <- 0.15
maxY_lim <- 150000
# Switch that enables restriction of the response factor plots to the data calculated as final linear range
restrictRfToFinalRange <- FALSE
# Scaling option is applied to x-axis of the response factor plots only
scalingRFX <- 'log' 
# Threshold (in percent) used to draw a horizontal line in the response factor plots 
# indicating a downward deviation from the mean response factor.
RfThresL <- 80  
# Threshold (in percent) used to draw a horizontal line in the response factor plots 
# indicating an upward deviation from the mean response factor. 
RfThresU <- 120 
# Factor that specifies the magnification of the axes annotations relative to the default (represented by the value 1)
magnificAxisRF = 1
# Factor that specifies the magnification of the axes labels relative to the default (represented by the value 1)
magnificLabelsRF = 1
# Factor that specifies the magnification of the title relative to the default (represented by the value 1)
magnificTitleRF = 1
# Resolution of the plots given in dpi
plotResolutionRF <- 300
# Height of the plot given in cm
plotHeightRF <- 17
# width of the plot given in cm
plotWidthRF <- 17

# Specific texts for the response factor plot 
titleOpeningRF <- "Response factor plot for "
xLabRF <- "Log concentration [fmol/µl]"
yLabRF <- "Response factor"
legend_RF <- c("blue: concentration levels outside of the linear range","red: concentration levels of the linear range")
# Specifies whether or not a legend is shown in the response factor plots
showLegendRF <- FALSE

# Further settings
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Enables the verbose mode of CalibraCurve. This mode is mainly useful for debugging
verbose <- FALSE 
# Control for the number of used decimal places in the result files
numberDecimals <- 3
# The value for 'fileTypeOfResults' specifies the file extension of the result files
fileTypeOfResults <- "txt"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Several texts used in the context of logging and exception handling
# (can be adapted if CalibraCurve is used different, non-targeted proteomics context)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
logText1 = "Current file name is: "
# Error and warning messages
fatalErrorPathIn <- c("Fatal error. CalibraCurve is shutting down. Something went wrong with the settings for the input directory. Please check.")
fatalErrorPathOut<- c("Fatal error. CalibraCurve is shutting down. Something went wrong with the settings for the output directory. Please check.")
# IO issues
fatalError <- c("Fatal error. Please find a detailed explanation in the log file.")
WritingErrorResGen <- c("Error while creating the result files.")
WritingErrorResF <- c("Error while writing in the result file starting with res_.") 
WritingErrorResRegF <- c("Error while writing in the result file starting with res_regr_.")
WritingErrorLogF <- c("Error while writing into the log file.")
WritingErrorPlot <- c("Error while generating or writing the CalibraCurve plots.")

# Messages for the R console
# Errors are issues that prevent further execution of CalibraCurve
genericWarningOnShell <- gsub("[\r\n]", "", c("An error occurred. Please find more detailed explanation in the log file."))
# Inconsistencies are data related issues, even though CalibraCurve analysis is possible
genericWarningOnShell2 <- gsub("[\r\n]", "", c("Some inconsistencies exist. Please find more detailed explanation (indicated by the prefix Warning) 
in the result file (i.e. the file starting with res_)."))

# Error-Messages given for issues that prevent further execution of CalibraCurve for the current input data (written into the log file)
errorInvalidDataStructure <- gsub("[\r\n]", "", c("Calculation for this sample has been stopped. 
The generated data structure is invalid for processing. Please check if the correct column delimiter is specified 
in the settings, because wrong settings are common causes for this error."))

errorMessCV <- gsub("[\r\n]", "", c("Calculation for this sample has been stopped. 
The preliminary linear range comprises less than two CONSECUTIVE concentration levels 
due to CV values that do not pass the threshold. 
Computing a regression, calculation of percent bias values or plotting is not meaningful."))

errorMessCV2 <- gsub("[\r\n]", "", c("Calculation for this sample has been stopped. The total number of concentration levels 
with CV values below the threshold is less than two. This allows no meaningful 
CalibraCurve analysis."))

errorMessPerBias <- gsub("[\r\n]", "", c("Calculation for this sample has been stopped. 
Because of the application of percent bias thresholds, less than two concentration levels 
remain within the final linear range. 
Computing a regression or plotting is not meaningful."))

# Warning messages for existing inconsistencies (written into the result file) 
warnMessCVRanges <- gsub("[\r\n]", "", c("Warning: More than one range with the SAME number of continuous concentration levels that 
pass the user defined CV threshold exist. 
The COMPLETE range between the lowest and the highest concentration level that have CV values < threshold is selected for further analysis. 
There are maybe intermediate levels that fail the CV threshold criteria."))

warnMessCVRanges2 <- gsub("[\r\n]", "", c("Warning: More than one continuous level range exist with CVs < threshold. 
The range with the highest number of levels is selected for further data analysis. 
It is recommended to review the data in order to identify possible outliers. 
Another possibility is to apply the calcContinuousPrelimRanges setting to 'calcContinuousPrelimRanges = FALSE'. 
This extends the preliminary range from the lowest to the highest concentration level that pass the CV threshold criteria."))

warnMessAcc <- gsub("[\r\n]", "", c("Warning: Some intermediate concentration levels (between LLOQ and ULOQ) show average values 
that fail the percent bias threshold. Please find them highlighted in the table below."))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Texts used for writing the result files (can be adapted if CalibraCurve is used different contexts (e.g. not related to targeted proteomics)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
resFHeader <- "Calculated concentration, CV values and percent bias for "
resRegFHeader <- "Results of the regression analyses for "
resFCVHeader <- "CV values"
resFCVPrefix <- "CV for concentration level "
perBiasOverviewWarning = "Warning: average percent bias is outside of the threshold"
headerPerBias <- "Calculated percent bias and response factor values:"
headerOverviewPerBiasW <- "Descriptive statistics for the percent bias values (based on the weighted linear model):"
headerOverviewPerBiasUW <- "Descriptive statistics for the percent bias values (based on the unweighted linear model):"
headerMeanRf <- "Mean response factors for the concentration levels of the final linear range"
#################################
### CalibraCurve initialization
#################################

tryCatch(setwd(inputDir), error = function(e) stop(fatalErrorPathIn), warning = function(w) stop(fatalErrorPathIn))
# Determination of the data file names.
filenames <- dir(pattern=paste("*.", fileTypeOfData, sep=""))
# Preparation of the log file
tryCatch(setwd(outputDir), error = function(e) stop(fatalErrorPathOut), warning = function(w) stop(fatalErrorPathOut))
logFile <- paste("log.txt", sep="")
tryCatch( 
  write(paste("Start of CalibraCurve at :", date(), sep = ""), 
        file = logFile, append = TRUE),
  error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
setwd(inputDir)

#################################
### Functions
#################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Note: in the first part of the function description the working step that 
# applies the function is specified.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Data preprocessing: Function, which selects all rows with identical concentrations
selConcentrationLevels <- function(x, rawData, colNumbConc){
  result <- rawData[rawData[,colNumbConc]==x,]
  return(result)
}

# Data preprocessing: Function, which checks for sufficent number of replicates
checkNumberReplicates <- function(x, data, minNumber){
  if(dim(data[[x]])[1]<minNumber) {
    result <- FALSE } 
  else {
    result <- TRUE }
  return(result)
}

# Preliminary linear range: Function, which calculates a CV value for a each concentration level
# Return value is a data frame, which is named with the concentration levels
calcCV <- function(x, dataSet, colNumMeas, colNumConc){
  actConcLevData <- dataSet[[x]]
  actSD <-  sd(actConcLevData[, colNumMeas])
  actMean <- mean(actConcLevData[, colNumMeas])
  CV <- actSD/actMean*100
  names(CV) <- unique(actConcLevData[, colNumConc])
  result <- CV
  return(result)
}

# Preliminary linear range: The function retuns a concentration level
currConcentrationLevels <- function(x, dataSet, colNumbConc){
  actConcLevData <- dataSet[[x]]
  result <- unique(actConcLevData[, colNumbConc])
  return(result)
}

# Preliminary linear range: Calculate key features for the existing continuous preliminary ranges
calcContPrelimRanges  <- function(x){
  # Calculating the end points for all ranges
  y <- x
  z <- x
  y <- y[-1]
  z <- z[-length(z)]
  differences <- y-z
  endPoints <- which(differences > 1)
  # Calculating the start points for all ranges
  startPoints <- endPoints+1
  # Add special case: 
  # the highest end point
  endPoints <- c(endPoints, length(x))
  # the lowest start point
  startPoints <- c(1, startPoints)
  # Calculating the extent of the ranges
  extent <- endPoints - startPoints+1
  startPoints <- x[startPoints]
  endPoints <- x[endPoints]
  rangeProperties <- as.data.frame(cbind(startPoints, endPoints, extent))
  result <- rangeProperties
}

# Final linear range: Function, which calculates weights for a specific concentration level
calcWeights <- function(x, dataSet, weightMet, colNumConc){
  currConcLevData <- dataSet[[x]]
  if(weightMet=='1/x^2'){
    weightsCurrLev <- 1/currConcLevData[,colNumConc]^2
  }else{
    if(weightMet=='1/x'){
      weightsCurrLev <- 1/currConcLevData[,colNumConc]
    }else{
      # weightsCurrLev <- # Please define your own weighting method here,
                          # it is also possible to provide a numeric vector with weights for each entry of dataFinal.
    }
  }
  result <- weightsCurrLev
  return(result)
}

# Final linear range: Function, which fits a linear model to the current data
calcLinearModel <- function(x, colNumConc, colNumMeas, applyWeights, w){
  dataSetDF <- do.call(rbind, x) 
  if(applyWeights){
    lmfit <- lm(dataSetDF[,colNumMeas] ~ dataSetDF[,colNumConc], weights=w)
  }else{
    lmfit <- lm(dataSetDF[,colNumMeas] ~ dataSetDF[,colNumConc])
  }
  result <- lmfit
  return(result)
}

# Final linear range: Function, which calculates the percent bias value for a single data point 
calcPerBias <- function(x, LMfit, expConc){
  coeff <- LMfit$coefficients
  calculatedConc <- (x-coeff[1])/coeff[2]
  # Calculation of distances
  d <- abs(calculatedConc-expConc)
  perBias <- 100/expConc*d
  names(perBias) <- NULL
  result <- perBias
  return(result)
}

# Final linear range: Function, which returns a list with percent bias values for a data set (given as list) 
calcPerBiasLevels <- function(x, LMres, colNumMeas){
  concentrations <- as.numeric(names(x))
  perBiasValuesList <- NULL
  for(i in seq_along(x)){
    currConcLevData <- x[[i]]
    expectConc <- concentrations[i]
    currConcLevMeas <- currConcLevData[, colNumMeas]
    perBiasValuesList[i] <- list(sapply(currConcLevMeas, calcPerBias, LMfit=LMres, expConc=expectConc))
  }
  return(perBiasValuesList)
}

# Final linear range: The function returns a data frame consisting of columns for average percent bias 
# and for the related standard deviation and coefficient of variation
calcPerBiasAvgSDCV <- function(x, ctm){
  avgPerBias <- NULL
  stdDevPerBias <- NULL
  CV_PerBias <- NULL
  for(i in seq_along(x)){
    # Average percent bias value
    meanPerBias <- mean(x[[i]])
    if (ctm == 'mean'){
      currAvgPerBias <- meanPerBias
    }
    if (ctm == 'median'){
      currAvgPerBias <- median(x[[i]])
    }
    # Standard deviation of the percent bias values
    currStdDev <- sd(x[[i]])
    currCV <- currStdDev/meanPerBias*100
    avgPerBias <- c(avgPerBias, currAvgPerBias)
    stdDevPerBias <- c(stdDevPerBias, currStdDev) 
    CV_PerBias <- c(CV_PerBias, currCV)
  }
  result <- data.frame(avgPerBias, stdDevPerBias, CV_PerBias)
  return(result)
}

# Final linear range: Function, which calculates a response factor for a single data point
# Formula obtained from:  Green, J. M., A practical guide to analytical method validation. 
#                         Analytical Chemistry 1996, 68, 305A-309A.
calcResponseFactors <- function(x, intercept, expConc){
  result <- (x-intercept)/expConc
  return(result)
}

# Final linear range: Function, which returns a list with response factor values for a data set (given as list) 
calcRFLevels <- function(x, interc, colNumMeas){
  concentrations <- as.numeric(names(x))
  rfValuesList <- NULL
  for(i in seq_along(x)){
    currConcLevData <- x[[i]]
    expectConc <- concentrations[i] 
    currConcLevMeas <- currConcLevData[, colNumMeas]
    rfValuesList[i] <- list(sapply(currConcLevMeas, calcResponseFactors, intercept = interc, expConc=expectConc))
  }
  return(rfValuesList)
}

# Final linear range: The function returns mean response factor values 
calcRFMeans <- function(x){
  avgRF <- NULL
  for(i in seq_along(x)){
    avgRFCurrLevel <- mean(x[[i]])
    avgRF <- c(avgRF, avgRFCurrLevel)
  }
  #result <- data.frame(avgPerBias, stdDev, CV)
  return(avgRF)
}

# Final linear range: Function, which checks whether the final linear range has been reached
checkFinalRange <- function(perBiasInfoWeighted, perBiasInfoUnweighted){
  bothLevelsPassed <- FALSE
  lowLevelPassed <- FALSE
  highLevelPassed <- FALSE
  if(finalRangeCalculationMethod == 'unweighted_linear_model'){ 
    if(perBiasInfoUnweighted$avgPerBias[1]<=perBiasThres){
      lowLevelPassed <- TRUE
    }  
    if(perBiasInfoUnweighted$avgPerBias[length(perBiasInfoUnweighted$avgPerBias)]<=perBiasThres){
      highLevelPassed <- TRUE
    }  
  }else{# Usage of average percent bias values calculated from the weighted linear model 
    if(perBiasInfoWeighted$avgPerBias[1]<=perBiasThres){
      lowLevelPassed <- TRUE
    }  
    if(perBiasInfoWeighted$avgPerBias[length(perBiasInfoWeighted$avgPerBias)]<=perBiasThres){
      highLevelPassed <- TRUE
    }      
  }
  if ((lowLevelPassed == TRUE) && (highLevelPassed == TRUE)){
    bothLevelsPassed <- TRUE
  }
  result <- bothLevelsPassed
  return(result)
}

# Final linear range: Auxillary function - compares the average percent bias
compDistances <- function(hDist, lDist){
  selRes <- NULL
  if(lDist > hDist){
    selRes <- TRUE
  }
  if(lDist < hDist){
    selRes <- FALSE
  }
  if(lDist == hDist){
    # if the values are the same, the lowest level will be removed
    selRes <- TRUE 
  }
  result <- selRes
}

# Final linear range: Function, which selects a concentration level for subsequent removal
selctConcLevel <- function(x, consPerBiasCV, perBiasT, perBiasDistT){
  removeLow <- NULL
  featuresLowestLevel <- x[1,]
  featuresHighestLevel <- x[dim(x)[1],]
  lowErrorPercent <- featuresLowestLevel$avgPerBias
  highErrorPercent <- featuresHighestLevel$avgPerBias
  lowFails <- (lowErrorPercent >= perBiasT)
  highFails <- (highErrorPercent >= perBiasT)
  if(lowFails && highFails){# Case: both the high and low level fails the criteria for the average percent bias value
    dist <- abs(lowErrorPercent-highErrorPercent)
    if(!consPerBiasCV){
      removeLow <- compDistances(hDist = highErrorPercent, lDist = lowErrorPercent)
    }else{
      # Here, the variance of the percent bias values gets a special consideration for the selection of 
      # the concentration level
      if(dist < perBiasDistT){ # CV values are only considered for level selection if the distance is low
        if(featuresLowestLevel$CV > featuresHighestLevel$CV){
          removeLow <- TRUE
        }
        if(featuresLowestLevel$CV < featuresHighestLevel$CV){
          removeLow <- FALSE
        }
        if(featuresLowestLevel$CV == featuresHighestLevel$CV){
          removeLow <- TRUE
        }
      }else{
        removeLow <- compDistances(hDist = highErrorPercent, lDist = lowErrorPercent)  
      }
    }
  }else{
    # Case: only one levels fails the threshold criteria
    if(lowFails){
      removeLow <- TRUE     
    }else{
      removeLow <- FALSE 
    }
  }
  result <- removeLow
}

# Writing result files: Preparation of result files
preparationResultFiles <- function(x, resF, resRegF){
  write(paste(resFHeader, x, sep = ""), 
        file = resF, append = FALSE)
  sink(resRegF, append = FALSE)
  write(paste(resRegFHeader, x, sep = ""), 
        file = resRegF, append = FALSE)
  sink()
}

# Writing result files: CV values for the preliminar linear range
writeCV_forPreliminaryLinearRange  <- function(x, resF, cUnit){
  conc <- names(x)
  write("", resF, append=TRUE)
  write(resFCVHeader, resF, append=TRUE)
  for(i in seq_along(x)){
    textCV <- paste(resFCVPrefix,
                    conc[i], cUnit, ": ",round(x[i], digits = numberDecimals),sep="")
    write(textCV, resF, append=TRUE)
  }
}

# Writing result files: Function, which writes percent bias and response factor values for the final linear range 
# or for writing of intermediate percent bias calculations into the logging file (CalibraCurve verbose mode)
writePerBiasRfFinalLinearRange <- function(x, perBiasWeighted, perBiasUnweighted, Rf, resF, colDelim, noDec){
  for(i in seq_along(dataFinal)){
    PerBias_no_weight<- perBiasUnweighted[[i]]
    PerBias_weight <- perBiasWeighted[[i]]
    ResponseFactors <- Rf[[i]]
    x[[i]] <- cbind(x[[i]], PerBias_no_weight, PerBias_weight, ResponseFactors)
  }
  xDF <- do.call(rbind, x)
  header <- paste(colnames(xDF), collapse = colDelim)
  write(header, resF, append=TRUE)
  xDF <- format(xDF, digits = noDec, nsmall = noDec)
  write.table(xDF, resF, append=TRUE, row.names=FALSE, col.names=FALSE, sep=colDelim)
}

# Writing result files: Function, which writes overview information (central tendency and dispersion
# measures) for the percent bias values of each concentration level
writePerBiasOverviewFinalLinearRange  <- function(x, resF, colDelim, message, noDec, perBiasT){
  indices <- unique(sort(which(x$avgPerBias >= perBiasT), decreasing = FALSE))
  warnings <- NULL
  for(i in seq_along(x$avgPerBias)){
    warnings <- c(warnings, "")
  }
  warnings[indices] <- message 
  x <- cbind(x, warnings)
  header <- paste(colnames(x), collapse = colDelim)
  header <- paste("expectedConc",header, sep = colDelim)
  if(length(indices)>0){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell2, sep="")
    print(mess)
    write(mess, logFile, append=TRUE)
    write(warnMessAcc, resF, append=TRUE)  
  }
  write(header, resF, append=TRUE)
  x <- format(x, digits = noDec, nsmall = noDec)
  write.table(x, resF, append=TRUE, row.names=TRUE, col.names=FALSE, sep=colDelim)
}

# Writing result files: Function, which writes mean response factor values for the 
# concentration levels of the final data range 
writeMeanRFvalues <- function(x, resF, colDelim, noDec){
  c <- paste(format (as.numeric(names(x)), digits = noDec, nsmall = noDec), collapse = colDelim)
  c <- paste("Expected concentration", c, sep = colDelim)
  write(c, resF, append=TRUE)
  x <- format(x, digits = noDec, nsmall = noDec)
  rf <- paste(x, collapse = colDelim)
  rf <- paste("Mean response factor", rf, sep = colDelim)
  write(rf, resF, append=TRUE) 
}

# Writing result files: Function, which writes summary information for the fitted linear model into the second result file
writeLMSummary  <- function(x, weightedLM, resRegF, weightMet){
  # unweighted model
  write("Summary for the unweighted linear model:", resRegF, append = TRUE)
  write("------------------------------------------", resRegF, append = TRUE)
  sink(resRegF, append = TRUE) 
  print(summary(x))
  sink() 
  # Weighted model    
  header <- paste("Summary for the weighted linear model (Weighting method:",weightMet,"):", sep = "")
  write(header, resRegF, append = TRUE)
  write("------------------------------------------", resRegF, append = TRUE)
  sink(resRegF, append = TRUE)
  print(summary(weightedLM))
  sink()
}

# Plotting
generatePlots <- function(x, titlePart1, sampleName, titlePart2, LLOQ, ULOQ, unit, 
                          x_lab, y_lab, lmW, lmNoW, plotOpt, scale,
                          colNumConc, colNumMeas){
  tiff(paste(sampleName, ".tif", sep = ""), width = plotWidth, height = plotHeight, units = "cm", res = plotResolution)
  par(cex.axis = magnificAxis, cex.lab = magnificLabels, cex.main = magnificTitle)
  if(showLegend){
    par(xpd=T, mar=par()$mar+c(6,0,0,0))
  }
  title <- paste(titlePart1, sampleName, '\n',
                 titlePart2, format(as.numeric(LLOQ), digits=numberDecimals, nsmall =numberDecimals),"-", 
                 format(as.numeric(ULOQ), digits=numberDecimals, nsmall =numberDecimals), " ",unit, sep = "")
  if(scale == 'log'){
    plot(x[,colNumConc], x[,colNumMeas], main = title, 
         xlab = x_lab, ylab = y_lab, log = "xy")
  }else{
    plot(x[,colNumConc], x[,colNumMeas], main = title, 
         xlab = x_lab, ylab = y_lab)      
  }
  
  r2NOW=summary(lmNoW)$r.squared
  r2W=summary(lmW)$r.squared
  # Unweighted regression equation
  if (lmNoW$coefficients[1]>0){
    eq <- paste("y=", format(lmNoW$coefficients[2], digits=2, scientific=TRUE),"*x+",
                format(lmNoW$coefficients[1], digits=2, scientific=TRUE),",",
                " R^2=",format(r2NOW, digits=3), sep="")
  }else{
    eq <- paste("y=", format(lmNoW$coefficients[2], digits=2, scientific=TRUE),"*x",
                format(lmNoW$coefficients[1], digits=2, scientific=TRUE),",", 
                " R^2=",format(r2NOW, digits=3), sep="")
  }
  # Weighted regression equation
  if (lmW$coefficients[1]>0){
    eq_w <- paste("y=", format(lmW$coefficients[2], digits=2, scientific=TRUE),"*x+",
                  format(lmW$coefficients[1], digits=2, scientific=TRUE),",",
                  " R^2=",format(r2W, digits=3), sep="")
  }else{
    eq_w <- paste("y=", format(lmW$coefficients[2], digits=2, scientific=TRUE),"*x",
                  format(lmW$coefficients[1], digits=2, scientific=TRUE),",",
                  " R^2=",format(r2W, digits=3), sep="")
  }
  if(plotOpt == "both"){ # including both weighted and unweighted regression lines and corresponding equations into the plot
    # Regression lines
    abline(lmNoW, col="blue", untf=TRUE, xpd=FALSE) 
    abline(lmW, col="red", untf=TRUE, xpd=FALSE)
    # Equations
    # Unweighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq, adj = c( 0, 1 ), col = "blue", cex = magnificEquation)
    if (scale == 'log'){
      # Weighted equation: below the unweighted equation
      text(min(x[,colNumConc]), max(x[,colNumMeas])/2, eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation)   
    }else{
      text(max(x[,colNumConc])/2, max(x[,colNumMeas]), eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation)  
    }
    # Legend
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_weighted, legend_unweighted), 
             inset=c(0,-0.45), # location 2: vertical position
             col=c("red","blue"),
             lty = c(1,1))
    }
    
  }
  if(plotOpt == "regW"){ # including only the weighted regression line into the plot
    abline(lmW, col="red", untf=TRUE, xpd=FALSE)
    # Weighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation )
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_weighted), 
             inset=c(0,-0.45), # location 2: vertical position
             col=c("red"),
             lty = c(1))
    }
  }
  if(plotOpt == "regUW"){ # including only unweighted regression line into the plot
    abline(lmNoW, col="blue", untf=TRUE, xpd=FALSE) 
    # Unweighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq, adj = c( 0, 1 ), col = "blue", cex = magnificEquation )
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_unweighted), 
             inset=c(0,-0.45), # location 2: vertical position
             col=c("blue"),
             lty = c(1))
    }
  }
  dev.off()
}

# Plotting: The function generates a response factor plot (Response factors vs. concentration) and presents
# both response factors for each measurement and mean response factors for each concentration level
# The plot also shows lines indicating user defined thresholds for percent deviation from the mean
# response factor that is calculated from the data of the final linear range.
generateRfPlot <- function(x, RfDataValid,  RfTLowF, RfTUpperF, avgRFDataF, avgRFDataV, restrictRF,
                           adaptYlim, minYlim, maxYlim, 
                           scaleRFX, magnAxisRF, magnLabelsRF,
                           magnTitleRF, plotResRF, plotHRF, plotWRF,
                           titleOpenRF, filename, xLRF, yLRF){
  # Calculation of the mean RF from the RF of data final
  allRf <- NULL
  for (i in seq_along(x)){
    currRf <- x[[i]]
    allRf <- c(allRf, currRf)
  }
  meanRfDF <- mean(allRf)
  # Calculating values for the horizontal threshold lines
  hLineLow <- meanRfDF*RfTLowF
  hLineUpper <- meanRfDF*RfTUpperF
  # Prepare data for plotting from the validated data
  plotReadyDataValid <- NULL
  for (i in seq_along(RfDataValid)){
    plotReadyDataValid <- rbind(plotReadyDataValid, calculateRfDf(RfDataValid[i]))
  }
  # Prepare data for plotting from the final data
  plotReadyDataFinal <- NULL
  for (i in seq_along(x)){
    plotReadyDataFinal <- rbind(plotReadyDataFinal, calculateRfDf(x[i]))
  }
  # Calculation of ylim properties
  if(restrictRF == TRUE){
    ylimMin <- min(plotReadyDataFinal$RF)
    ylimMax <- max(plotReadyDataFinal$RF)  
  }else{
    ylimMin <- min(plotReadyDataValid$RF, plotReadyDataFinal$RF)
    ylimMax <- max(plotReadyDataValid$RF, plotReadyDataFinal$RF)    
  }
  if(hLineLow < ylimMin){
    ylimMin = hLineLow
  }
  if(hLineUpper > ylimMax){
    ylimMax = hLineUpper
  }
  if(adaptYlim){# Enables user defined adaption of the y axis
    ylimMin <- minYlim
    ylimMax <- maxYlim
  }
  rfPlotName <- paste("RFplot_", filename, ".tif", sep="")
  tiff(rfPlotName, width = plotWRF, height = plotHRF, units = 'cm', res = plotResRF)
  par(cex.axis = magnAxisRF, cex.lab = magnLabelsRF, cex.main = magnTitleRF)
  if(showLegendRF){
    par(xpd=T, mar=par()$mar+c(6,0,0,0))
  }
  # Case: Plot response factors for the complete input data
  title <- paste(titleOpenRF, filename, sep="")
  if (restrictRF == FALSE){
    if (scaleRFX == "log"){
      plot(plotReadyDataValid, log = "x", col = 'blue', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)    
    }else{
      plot(plotReadyDataValid, col = 'blue', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)   
    }
  }else{
    # Recalculation of ylim
    if (scaleRFX == "log"){
      plot(plotReadyDataFinal, log = "x", col = 'red', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)    
    }else{
      plot(plotReadyDataFinal, col = 'red', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)   
    }      
  }
  if (restrictRF == FALSE){
    points(plotReadyDataFinal$Concentration, plotReadyDataFinal$RF, col="red")  
    lines(names(avgRFDataV), avgRFDataV, col="blue",lty=1, type = "o", pch = 19)
  }
  lines(names(avgRFDataF), avgRFDataF, col="red",lty=1, type = "o", pch = 19)
  abline(h = hLineLow, col = 'green', lty = 2, xpd=FALSE)
  abline(h = hLineUpper, col = 'green', lty = 2, xpd=FALSE)
  if(showLegendRF){
    legend("bottom", # location 1: centered at the bottom
           legend=c(legend_RF), 
           inset=c(0,-0.45), # location 2: vertical position
           col=c("blue", "red"),
           lty = c(1,1))
  } 
  dev.off()
}

# Auxillary function used for response factor plotting:
# The function generates a data frame with a concentration column and a response factor column
calculateRfDf <- function(x){
  Concentration <- NULL
  RF <- NULL
  conc <- as.numeric(names(x))
  for (i in x){
    Concentration <- c(Concentration, conc)
    RF <- c(RF, i)
  }
  mDF <- data.frame(Concentration, RF)
  return(mDF)
}

#################################
### Data processing
#################################
for (i in filenames){
  setwd(inputDir)
# Boolean variables that control error and warning messages
  errorResFile <- FALSE
  errorResRegFile <- FALSE
  errorPlot <- FALSE
  currFileNameNoExt <- unlist(strsplit(i, "[.]"))
# Reading the data
  data <- read.csv(i, header=TRUE, sep=colDelimiterIn, strip.white = TRUE)
# Checking whether the data structure is valid or not
  if (dim(data)[2]<2){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",fatalError, sep="")
    print(mess)
    tryCatch(write(errorInvalidDataStructure, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
    next
  }
# Removal of empty (NA) columns
  data <- Filter(function(x)!all(is.na(x)), data)
  setwd(outputDir)
  # Writing log information
  tryCatch(write(paste(logText1, currFileNameNoExt[1], sep = ""), file = logFile, append = TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
  # Variable for the file that contains the CV and percent bias calculations
  resFile <- paste("res_", currFileNameNoExt[1], ".", fileTypeOfResults, sep="")
  # Variable for the file that contains summary information for the linear regression 
  resRegFile <- paste("res_regr_", currFileNameNoExt[1], ".", fileTypeOfResults, sep="") 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Data preprocessing
# This working step includes several data quality checks for the columns that contain concentration
# information and measurement values:
# a) Occurrence of 0 values
# b) Occurrence missing values
# c) Sufficent number of replicates for each concentration level
# If necessary the data is cleansed accordingly.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Removing rows that contain unwanted 0 values
  dataValidated <- data[data[colNumberConcentration]!=0 | is.na(data[colNumberConcentration]),]
  dataValidated <- dataValidated[dataValidated[colNumberMeasurements]!=0 | is.na(dataValidated[colNumberMeasurements]),]
  # Removing rows with missing values 
  rowIndicesNA <- unique(c(which(is.na(dataValidated[,colNumberConcentration])),
                           which(is.na(dataValidated[,colNumberMeasurements]))))
  if(length(rowIndicesNA)>0){
    dataValidated <- dataValidated[-rowIndicesNA,]
  }
# Determination of existing concentration levels in the validated data
  concLevels <- unique(dataValidated[,colNumberConcentration]) 
  sort(concLevels, decreasing = FALSE)
# Transforming a data set into a list with entries for each concentration level (and the related data)
  dataValidated <- lapply(concLevels, FUN=selConcentrationLevels, rawData=dataValidated, colNumbConc=colNumberConcentration) 
# Deleting concentration levels with insufficient number of replicates
  dataValidated <- dataValidated[sapply(1:length(dataValidated), FUN=checkNumberReplicates, data=dataValidated, minNumber=minNumberReplications)]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Calculation of CV values and the resulting preliminary linear range
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Coefficient of variation values for all concentration levels of the validated data
  concLevelsCV <- sapply(1:length(dataValidated), calcCV, dataSet=dataValidated, colNumMeas=colNumberMeasurements, colNumConc=colNumberConcentration)
  index <- which(concLevelsCV <= cv_thres)
# Preparation of result files at this point, because writing of results (or info messages) starts
  mess <- paste("Input data file ",currFileNameNoExt[1],": ",WritingErrorResGen, sep="")
  tryCatch(preparationResultFiles(currFileNameNoExt[1], resF = resFile, resRegF =  resRegFile), 
           error = function(e) print(mess), warning = function(w) print(mess))
  if (length(index) < 2){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell, sep="")
    print(mess)
    tryCatch(write(errorMessCV2, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
    next # Calculation stops for this sample
  }
  propertiesOfRanges <- calcContPrelimRanges(index)
  indexStart <- min(propertiesOfRanges[,1])
  indexEnd <- max(propertiesOfRanges[,2])

  if(calcContinuousPrelimRanges==FALSE){# FALSE: Selecting all levels between the lowest and highest concentration level, where CV <threshold
    dataPrelim <- dataValidated[indexStart:indexEnd]
  }else{ # Computing the longest continuous preliminary range
    if (dim(propertiesOfRanges)[1] > 1) {
      mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell2, sep="")
      print(mess)
      tryCatch(write(mess, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
      tryCatch(write(warnMessCVRanges2, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
    }
    longestRange <- propertiesOfRanges[which(propertiesOfRanges$extent==max(propertiesOfRanges$extent)),]
    if(dim(longestRange)[1]>1){# Special case: more than one range with the same number of concentration levels that pass the CV threshold exist
      mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell2, sep="")
      print(mess)
      tryCatch(write(mess, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
      tryCatch(write(warnMessCVRanges, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)     
      dataPrelim <- dataValidated[indexStart:indexEnd]
    }else{
      if(longestRange$extent == 0){
        mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell, sep="")
        print(mess)
        tryCatch(write(errorMessCV, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))       
        next # Calculation stops for this sample
      }
      dataPrelim <- dataValidated[longestRange$startPoints:longestRange$endPoints]
    }
  }
  prelimConcLevels <- sapply(1:length(dataPrelim), currConcentrationLevels, dataSet=dataPrelim, colNumbConc=colNumberConcentration) 
  names(dataPrelim) <- prelimConcLevels
  concLevelsCVSel <- concLevelsCV[as.character(prelimConcLevels)]
# Cancel calculations if less than two remaining concentration levels exist, which pass the CV value threshold
  if (length(concLevelsCVSel) < 2){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell, sep="")
    print(mess)
    tryCatch(write(errorMessCV, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))       
    next # Calculation stops for this sample
  }
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Calculation of regressions, percent bias values and the final linear range
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  dataFinal <- dataPrelim
  finalRangeNotReached <- TRUE
  while(finalRangeNotReached){
    if(length(dataFinal)<2){
      break # Calculation stops for this sample
    }
    allWeights <- as.vector(sapply(1:length(dataFinal), FUN=calcWeights, 
                                   dataSet=dataFinal, weightMet=weightingMethod, 
                                   colNumConc=colNumberConcentration)) 
# Ensure that allWeights is a vector (bugfix: some data sets lead to creation of lists)
    if (class(allWeights)=="list"){
      allWeights <- unlist(allWeights)
    }
    lmWeighted <- calcLinearModel(dataFinal, colNumConc=colNumberConcentration, 
                                  colNumMeas=colNumberMeasurements, applyWeights=TRUE, w=allWeights)
    lmUnweighted <- calcLinearModel(dataFinal, colNumConc=colNumberConcentration, 
                                    colNumMeas=colNumberMeasurements, applyWeights=FALSE)
    perBiasWeighted <- calcPerBiasLevels(dataFinal, LMres=lmWeighted, 
                                          colNumMeas=colNumberMeasurements)
    names(perBiasWeighted) <- names(dataFinal)
    perBiasUnweighted <- calcPerBiasLevels(dataFinal, LMres=lmUnweighted, 
                                            colNumMeas=colNumberMeasurements)
    names(perBiasUnweighted) <- names(dataFinal)
    perBiasAvgSDCVWeighted <- calcPerBiasAvgSDCV(perBiasWeighted, ctm=centralTendencyMeasure)
    rownames(perBiasAvgSDCVWeighted) <- names(dataFinal)
    perBiasAvgSDCVUnweighted <- calcPerBiasAvgSDCV(perBiasUnweighted, ctm=centralTendencyMeasure)
    rownames(perBiasAvgSDCVUnweighted) <- names(dataFinal)
    if(checkFinalRange(perBiasInfoWeighted = perBiasAvgSDCVWeighted, perBiasInfoUnweighted = perBiasAvgSDCVUnweighted)){
      finalRangeNotReached<- FALSE
    }else{ # Removal of concentration levels from the data set
      if (finalRangeCalculationMethod == 'weighted_linear_model'){
        removeLow <- selctConcLevel(perBiasAvgSDCVWeighted, perBiasT = perBiasThres, 
                                    consPerBiasCV = considerPerBiasCV, perBiasDistT = perBiasDistThres)
      }
      if (finalRangeCalculationMethod == 'unweighted_linear_model'){
        removeLow <- selctConcLevel(perBiasAvgSDCVUnweighted, perBiasT = perBiasThres, 
                                    consPerBiasCV = considerPerBiasCV, perBiasDistT = perBiasDistThres)
      } 
      if(removeLow){
        dataFinal[[1]] <- NULL
      }else{
        dataFinal[[length(dataFinal)]] <- NULL
      }
    }
    if(verbose){
      remDir <- getwd()
      setwd(outputDir)
      writePerBiasRfFinalLinearRange(dataFinal, perBiasWeighted = perBiasWeighted,
                                   perBiasUnweighted = perBiasUnweighted, 
                                   Rf = resFacDataF,
                                   resF = logFile,
                                   colDelim = colDelimiterOut,
                                   noDec = numberDecimals)
      setwd(remDir)
    }
  }
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Calculating response factors 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(length(dataFinal)<2){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell, sep="")
    print(mess)    
    tryCatch(write(errorMessPerBias, logFile, append=TRUE), error = function(e) print(WritingErrorLogF), warning = function(w) print(WritingErrorLogF))
    next # Calculation stops for this sample
  }
  dataValConcLevels <- sapply(1:length(dataValidated), currConcentrationLevels, dataSet=dataValidated, colNumbConc=colNumberConcentration) 
  names(dataValidated) <- dataValConcLevels
  intercept <- lmWeighted$coefficients[1]
  names(intercept) <- NULL
  resFacDataF <- calcRFLevels(dataFinal, interc = intercept, colNumMeas = colNumberMeasurements)
  resFacDataV <- calcRFLevels(dataValidated, interc = intercept, colNumMeas =colNumberMeasurements)
  names(resFacDataF) <- names(dataFinal)
  names(resFacDataV) <- names(dataValidated)
# Calculation of mean response factor values
  avgResFacDataF <- calcRFMeans(resFacDataF)
  names(avgResFacDataF) <- names(dataFinal)
  avgResFacDataV <- calcRFMeans(resFacDataV)
  names(avgResFacDataV) <- names(dataValidated)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Writing result files and plotting
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Writing final linear range statement
  finalCL <- names(dataFinal)
  finalRange <- paste("Linear range is from ", format(as.numeric(finalCL[1],concUnit), digits=numberDecimals, nsmall =numberDecimals), 
                      " to ", format(as.numeric(finalCL[length(finalCL)]), digits=numberDecimals, nsmall =numberDecimals), concUnit,".",sep="") 
  tryCatch(write("", resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write(finalRange, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
# Writing data for the preliminary linear range (CV related results)
  tryCatch(writeCV_forPreliminaryLinearRange(concLevelsCVSel, resF = resFile, cUnit = concUnit), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
# Writing calcualted percent bias values (along with the input information)
  tryCatch(write("", resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write(headerPerBias, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(writePerBiasRfFinalLinearRange(dataFinal, perBiasWeighted = perBiasWeighted,
                                            perBiasUnweighted = perBiasUnweighted, 
                                            Rf = resFacDataF,
                                            resF = resFile,
                                            colDelim = colDelimiterOut,
                                            noDec = numberDecimals), 
           error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
# Writing overview information (for percent bias)
  tryCatch(write("", resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write(headerOverviewPerBiasW, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(writePerBiasOverviewFinalLinearRange(perBiasAvgSDCVWeighted,  
                                                resF = resFile,
                                                colDelim = colDelimiterOut,
                                                message = perBiasOverviewWarning,
                                                noDec = numberDecimals,
                                                perBiasT = perBiasThres),
           error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write("", resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write(headerOverviewPerBiasUW, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(writePerBiasOverviewFinalLinearRange(perBiasAvgSDCVUnweighted,  
                                                  resF = resFile,
                                                  colDelim = colDelimiterOut,
                                                  message = perBiasOverviewWarning,
                                                  noDec = numberDecimals,
                                                  perBiasT = perBiasThres), 
           error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write("", resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(write(headerMeanRf, resFile, append=TRUE), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  tryCatch(writeMeanRFvalues(avgResFacDataF, resF = resFile, colDelim = colDelimiterOut, noDec = numberDecimals), error = function(e) errorResFile <<- TRUE, warning = function(w) errorResFile <<- TRUE)
  if (errorResFile==TRUE){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",WritingErrorResF, sep="")
    print(mess)
    write(mess, logFile, append=TRUE)     
  }
  
# Writing summary information for the fitted linear models
  tryCatch(write("", resRegFile, append=TRUE), error = function(e) errorResRegFile <<- TRUE, warning = function(w) errorResRegFile <<- TRUE)
  tryCatch(writeLMSummary(lmUnweighted, lmWeighted, resRegF = resRegFile, weightMet=weightingMethod), error = function(e) errorResRegFile <<- TRUE, warning = function(w) errorResRegFile <<- TRUE)
  lmWeighted$coefficients
  if (errorResRegFile==TRUE){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",WritingErrorResRegF, sep="")
    print(mess)
    write(mess, logFile, append=TRUE)     
  }
# Creation of the plots
  fn <- currFileNameNoExt[1]
  lowerLimit <- finalCL[1]
  upperLimit <- finalCL[length(finalCL)]
  if(plotOptReg == 'combined'){
    # Warnings are ignored for plotting (there are often thrown because of omitting 0 from log plot)
    suppressWarnings(tryCatch(generatePlots(data, titlePart1 = titleOpening, sampleName = fn , titlePart2 = titleMiddlePart,
                               LLOQ = lowerLimit, ULOQ = upperLimit, unit = concUnit,
                               x_lab = xLab,  y_lab = yLab, lmW = lmWeighted, lmNoW = lmUnweighted,
                               plotOpt = 'both', scale = scaling, colNumConc = colNumberConcentration, 
                               colNumMeas = colNumberMeasurements), 
             error = function(e) errorPlot <<- TRUE))
  }
  if(plotOptReg == 'separate'){
    fn_regW <- paste(fn, "_regW", sep="")
    suppressWarnings(tryCatch(generatePlots(data, titlePart1 = titleOpening, sampleName = fn_regW , titlePart2 = titleMiddlePart,
                           LLOQ = lowerLimit, ULOQ = upperLimit, unit = concUnit,
                           x_lab = xLab,  y_lab = yLab, lmW = lmWeighted, lmNoW = lmUnweighted,
                           plotOpt = 'regW', scale = scaling, colNumConc = colNumberConcentration, 
                           colNumMeas = colNumberMeasurements), 
             error = function(e) errorPlot <<- TRUE))
    fn_regUnW <- paste(fn, "_regUnW", sep="")
    suppressWarnings(tryCatch(generatePlots(data, titlePart1 = titleOpening, sampleName = fn_regUnW , titlePart2 = titleMiddlePart,
                           LLOQ = lowerLimit, ULOQ = upperLimit, unit = concUnit,
                           x_lab = xLab,  y_lab = yLab, lmW = lmWeighted, lmNoW = lmUnweighted,
                           plotOpt = 'regUW', scale = scaling, colNumConc = colNumberConcentration, 
                           colNumMeas = colNumberMeasurements), 
             error = function(e) errorPlot <<- TRUE))
  }
# Generating response factor plots
  RfThresUFactor <- RfThresU/100
  RfThresLFactor <- RfThresL/100
  tryCatch(generateRfPlot(resFacDataF, RfDataValid = resFacDataV, RfTLow = RfThresLFactor, RfTUpper = RfThresUFactor, 
                          avgRFDataF = avgResFacDataF, avgRFDataV = avgResFacDataV, restrictRF = restrictRfToFinalRange,
                          adaptYlim =yLimAdaptRF, minYlim = minY_lim, maxYlim = maxY_lim,
                          scaleRFX = scalingRFX, magnAxisRF = magnificAxisRF, magnLabelsRF = magnificLabelsRF,
                          magnTitleRF = magnificTitleRF, plotResRF = plotResolutionRF, plotHRF = plotHeightRF, 
                          plotWRF = plotWidthRF, titleOpenRF = titleOpeningRF, filename = fn, xLRF = xLabRF, yLRF = yLabRF), 
           error = function(e) errorPlot <<- TRUE, warning = function(w) errorPlot <<- TRUE)
  if (errorPlot==TRUE){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",WritingErrorPlot, sep="")
    print(mess)
    write(mess, logFile, append=TRUE)     
  }
}