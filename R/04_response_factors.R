

#' Calculate Response factors
#'
#' @description
#' Function, which calculates a response factor for a single data point
#'
#'
#' @details
#' Formula obtained from:  Green, J. M., A practical guide to analytical method validation.
#                         Analytical Chemistry 1996, 68, 305A-309A.
#'
#'
#' @param x
#' @param intercept Intercept of the linear model
#' @param expConc
#'
#' @returns
#' @export
#'
#' @examples
calcResponseFactors <- function(x, intercept, expConc) {
  result <- (x - intercept) / expConc
  return(result)
}


#' Calculate Response factors
#'
#' @description
#' Final linear range: Function, which returns a list with response factor values for a data set (given as list)
#'
#' @param x
#' @param mod linear model object
#'
#' @returns
#' @export
#'
#' @examples
calcRFLevels <- function(x, mod) {
  interc <- unname(mod$coefficients[1])
  concentrations <- as.numeric(names(x))

  rfValuesList <- NULL
  for (i in seq_along(x)) {
    currConcLevData <- x[[i]]
    expectConc <- concentrations[i]
    currConcLevMeas <- currConcLevData[, "Measurement"]
    rfValuesList[i] <- list(sapply(currConcLevMeas, calcResponseFactors, intercept = interc, expConc = expectConc))
  }

  names(rfValuesList) <- concentrations
  return(rfValuesList)
}



#' Calculate Response factors
#'
#' @description The function returns mean response factor values
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
calcRFMeans <- function(x) {
  avgRF <- NULL
  for (i in seq_along(x)) {
    avgRFCurrLevel <- mean(x[[i]])
    avgRF <- c(avgRF, avgRFCurrLevel)
  }
  #result <- data.frame(avgPerBias, stdDev, CV)
  names(avgRF) <- names(x)
  return(avgRF)
}


