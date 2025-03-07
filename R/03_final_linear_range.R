

#' FLR: calculate weights for linear model for a specific concentration level
#'
#' @param x data frame
#' @param weightMet method for weighting (currently 1/x and 1/x^2 are supported, default is 1/x^2)
#'
#' @returns
#' @export
#'
#' @examples
calcWeights <- function(x, weightMet = "1/x^2"){
  currConcLevData <- x
  if (weightMet == '1/x^2') {
    weightsCurrLev <- 1/(currConcLevData$Concentration^2)
  }
  if (weightMet == '1/x') {
      weightsCurrLev <- 1/currConcLevData$Concentration
  }
  #else {
      ### TODO
      # weightsCurrLev <- # Please define your own weighting method here,
      # it is also possible to provide a numeric vector with weights for each entry of dataFinal.
  #  }
  result <- weightsCurrLev
  return(result)
}



#' FLR: calculate weighted or unweighted linear model
#'
#' @param x list of data frames (one entry for each concentration level)
#' @param weights vector of weights (default is NULL and will result in an unweighted model)
#'
#' @returns
#' @export
#'
#' @examples
calcLinearModel <- function(x, weights = NULL){

  ## combine list elements to a data set
  dataSetDF <- do.call(rbind, x)

  if (is.null(weights)) {  # unweighted model
    lmfit <- lm(Measurement ~ Concentration, data = dataSetDF)
  } else { # weighted model
    lmfit <- lm(Measurement ~ Concentration, data = dataSetDF, weights = weights)
  }
  return(lmfit)
}


#' FLR: calculate percent bias for a specific concentration level
#'
#' @param x data frame for a specific concentration level
#' @param LMfit linear model fit
#' @param expConc expected (known) concentration level
#'
#' @returns
#' @export
#'
#' @examples
calcPerBias <- function(x, LMfit, expConc){
  coeff <- LMfit$coefficients
  calculatedConc <- (x - coeff[1])/coeff[2]
  # Calculation of distances
  # d <- abs(calculatedConc - expConc)
  perBias <- 100 * abs(calculatedConc - expConc) / expConc
  names(perBias) <- NULL
  return(perBias)
}



#' FLR: calculate list of percent bias values for all concentration levels
#'
#' @param x list of data frames (one list entry for each concentration level)
#' @param LMfit linear model fit
#'
#' @returns list of vectors with percent bias values for each data point per concentration level
#' @export
#'
#' @examples
calcPerBiasLevels <- function(x, LMfit){
  concentrations <- as.numeric(names(x))
  perBiasValuesList <- NULL
  for (i in seq_along(x)) {
    currConcLevData <- x[[i]]
    expectConc <- concentrations[i]
    currConcLevMeas <- currConcLevData$Measurement
    perBiasValuesList[i] <- list(sapply(currConcLevMeas, calcPerBias, LMfit = LMfit, expConc = expectConc))
  }

  names(perBiasValuesList) <- names(x)
  return(perBiasValuesList)
}


#' FLR: calculate average, SD and CV percent bias for each concentration level
#'
#' @param x result of calcPerBiasLevels()
#' @param method method: "mean" or "median", default is "mean"
#'
#' @returns data frame with 3 columns: avgPerBias, stdDevPerBias, CV_PerBias
#'          each row is one concentration level
#' @export
#'
#' @examples
calcPerBiasAvgSDCV <- function(x, method = "mean") {
  avgPerBias <- NULL
  stdDevPerBias <- NULL
  CV_PerBias <- NULL

  for (i in seq_along(x)) {
    # Average percent bias value
    meanPerBias <- mean(x[[i]])
    if (method == 'mean') {
      currAvgPerBias <- meanPerBias
    }
    if (method == 'median') {
      currAvgPerBias <- median(x[[i]])
    }
    # Standard deviation of the percent bias values
    currStdDev <- sd(x[[i]])
    currCV <- currStdDev/meanPerBias*100
    avgPerBias <- c(avgPerBias, currAvgPerBias)
    stdDevPerBias <- c(stdDevPerBias, currStdDev)
    CV_PerBias <- c(CV_PerBias, currCV)
    ### TODO: does CV make sense if method is "median"?
  }
  result <- data.frame(avgPerBias, stdDevPerBias, CV_PerBias)
  rownames(result) <- names(x)
  return(result)
}


# Final linear range: Function, which checks whether the final linear range has been reached
#' FLR: checks if final linear range has been reached (compare average percent bias with threshold)
#'
#' @param perBiasInfo result of calcPerBiasAvgSDCV()
#' @param perBiasThres threshold for average percent bias
#'
#' @returns logical(1) TRUE if bost lowest and highest concentration level passed the check
#'          (and the final linear range is reached), FALSE otherwise
#' @export
#'
#' @examples
checkFinalRange <- function(perBiasInfo, perBiasThres = 10){

  bothLevelsPassed <- FALSE
  lowLevelPassed <- FALSE
  highLevelPassed <- FALSE

  if (perBiasInfo$avgPerBias[1] <= perBiasThres) {
    lowLevelPassed <- TRUE
  }
  if (perBiasInfo$avgPerBias[length(perBiasInfo$avgPerBias)] <= perBiasThres) {
    highLevelPassed <- TRUE
  }

  if (lowLevelPassed && highLevelPassed) {
    bothLevelsPassed <- TRUE
  }
  result <- bothLevelsPassed
  return(result)
}

################################################################################

# Final linear range: Auxillary function - compares the average percent bias
# compDistances <- function(hDist, lDist){
#   selRes <- NULL
#   if(lDist > hDist){
#     selRes <- TRUE
#   }
#   if(lDist < hDist){
#     selRes <- FALSE
#   }
#   if(lDist == hDist){
#     # if the values are the same, the lowest level will be removed
#     selRes <- TRUE
#   }
#   result <- selRes
# }



#' FLR: selects a concentration level for subsequent removal
#'
#' @param x result of calcPerBiasAvgSDCV()
#' @param consPerBiasCV consider CV?
#' @param perBiasT threshold for average percent bias
#' @param perBiasDistT threshold for the difference in average percent bias (for lower differences, CV will be considered)
#'
#' @returns logical(1): TRUE, if lowest concentration will be removed, FALSE if highest will be removed
#'
#' @export
#'
#' @examples
selctConcLevel <- function(x, consPerBiasCV, perBiasT, perBiasDistT) {
  removeLow <- NULL
  featuresLowestLevel <- x[1,]
  featuresHighestLevel <- x[nrow(x),]

  lowErrorPercent <- featuresLowestLevel$avgPerBias
  highErrorPercent <- featuresHighestLevel$avgPerBias

  lowFails <- (lowErrorPercent >= perBiasT)
  highFails <- (highErrorPercent >= perBiasT)

  if (lowFails && highFails) {
    # Case: both the high and low level fails the criteria for the average percent bias value

    if (!consPerBiasCV) {
      removeLow <- ifelse(lowErrorPercent >= highErrorPercent, TRUE, FALSE)
    } else {
      dist <- abs(lowErrorPercent - highErrorPercent)
      # Here, the variance of the percent bias values gets a special consideration for the selection of
      # the concentration level
      if (dist <= perBiasDistT) { # CV values are only considered for level selection if the distance is low
        removeLow <- ifelse(featuresLowestLevel$CV >= featuresHighestLevel$CV, TRUE, FALSE)
      } else {
        removeLow <- ifelse(lowErrorPercent >= highErrorPercent, TRUE, FALSE)
      }
    }

  } else {
    # Case: only one levels fails the threshold criteria
    removeLow <- ifelse(lowFails, TRUE, FALSE)
  }
  return(removeLow)
}




#' Calculate the final linear range
#'
#' @param dataPrelim data frame with concentration levels and corresponding response values (result from calculate_PLR)
#' @param weightingMethod method for weighting (currently "1/x", "1/x^2" and "None" are supported, default is 1/x^2)
#' @param centralTendencyMeasure "mean" or "median" (for calculating average percent bias), default is "mean"
# #@param finalRangeCalculationMethod method for calculating the final linear range ("weighted_linear_model" or "unweighted_linear_model")
#' @param perBiasThres threshold for average percent bias, default is 20%
#' @param considerPerBiasCV consider CV for the selection of the concentration level, default is TRUE
#' @param perBiasDistThres threshold for the difference in average percent bias (for lower differences, CV will be considered), default is 10%
#'
#' @returns
#' @export
#'
#' @examples
calculate_FLR <- function(dataPrelim,
                          weightingMethod = "1/x^2",
                          centralTendencyMeasure = "mean",
                          #finalRangeCalculationMethod = "weighted_linear_model",
                          perBiasThres = 20,
                          considerPerBiasCV = TRUE,
                          perBiasDistThres = 10) {

  dataFinal <- dataPrelim
  finalRangeReached <- FALSE

  while (!finalRangeReached) {
    if (length(dataFinal) < 2) {
      break # Calculation stops because only one concentration level is left
    }

    ## calculate the weights:
    if (weightingMethod != "None") {
      allWeights <- as.vector(sapply(dataFinal,
                                     FUN = calcWeights,
                                     weightMet = weightingMethod))
      # Ensure that allWeights is a vector (bugfix: some data sets lead to creation of lists)
      if (class(allWeights) == "list") {
        allWeights <- unlist(allWeights)
      }
    } else {
      allWeights <- NULL
    }


    ## calculate weighted and unweighted linear model
    mod <- calcLinearModel(dataFinal, weights = allWeights)
    #lmUnweighted <- calcLinearModel(dataFinal, weights = NULL)

    ## calculate the percent bias for each data point
    perBias <- calcPerBiasLevels(dataFinal, LMfit = mod)
    #perBiasUnweighted <- calcPerBiasLevels(dataFinal, LMfit = lmUnweighted)

    ## calculate the average percent bias, standard deviation and CV for each concentration level
    perBiasAvgSDCV <- calcPerBiasAvgSDCV(perBias, method = centralTendencyMeasure)
    #perBiasAvgSDCVUnweighted <- calcPerBiasAvgSDCV(perBiasUnweighted, method = centralTendencyMeasure)


    #if (finalRangeCalculationMethod == 'weighted_linear_model') {
    checkFR <- checkFinalRange(perBiasInfo = perBiasAvgSDCV, perBiasThres = perBiasThres)
    #}
    #if (finalRangeCalculationMethod == 'unweighted_linear_model') {
    #  checkFR <- checkFinalRange(perBiasInfo = perBiasAvgSDCVUnweighted, perBiasThres = perBiasThres)
    #}

    if (checkFR) {
      finalRangeReached <- TRUE
    } else { # Removal of concentration levels from the data set
      #if (finalRangeCalculationMethod == 'weighted_linear_model') {
      removeLow <- selctConcLevel(perBiasAvgSDCV, perBiasT = perBiasThres,
                                  consPerBiasCV = considerPerBiasCV, perBiasDistT = perBiasDistThres)
      #}
      #if (finalRangeCalculationMethod == 'unweighted_linear_model') {
      #  removeLow <- selctConcLevel(perBiasAvgSDCVUnweighted, perBiasT = perBiasThres,
      #                              consPerBiasCV = considerPerBiasCV, perBiasDistT = perBiasDistThres)
      #}
      if (removeLow) {
        dataFinal[[1]] <- NULL
      } else {
        dataFinal[[length(dataFinal)]] <- NULL
      }
    }
  }
  return(list(dataFinal = dataFinal, mod = mod, perBias = perBias, perBiasAvgSDCV = perBiasAvgSDCV))
}

















