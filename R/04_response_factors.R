

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
#' @param x *data.frame** \cr Data.frame containing data for a specific concentration level.
#' @param intercept  **numeric(1)** \cr Intercept of the linear model.
#' @param expConc **numeric(1)** \cr  Expected concentration (known concentration value).
#'
#' @returns vector of response factors for this specific concentration level
calcResponseFactors <- function(x, intercept, expConc) {
  result <- (x - intercept) / expConc
  return(result)
}


#' Calculate Response factors
#'
#' @description
#' Final linear range: Function, which returns a list with response factor values for a data set (given as list)
#'
#' @param x **list of data.frames** \cr List of data.frames containing data for each concentration level (result from \code{\link{cleanData}}).
#' @param mod **lm object** \cr Final linear model fit (object "mod" from results of \code{\link{calculate_FLR}}).
#'
#' @returns List of response factor values for each concentration level.
#' @export
#'
#' @examples
#' data(D_MFAP4)
#' D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
#' RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
#'               cv_thres = 10,
#'               calcContinuousPrelimRanges = TRUE)
#' RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)
#'
#' calcRFLevels(D_MFAP4_cleaned, mod = RES_FLR$mod)
#'
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
#' @param x **list** \cr Result of \code{\link{calcRFLevels}} \cr List of response factor values for each concentration level.
#'
#' @returns vector of mean response factor values for the different concentration levels.
#' @export
#'
#' @examples
#' #' data(D_MFAP4)
#' D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
#' RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
#'               cv_thres = 10,
#'               calcContinuousPrelimRanges = TRUE)
#' RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)
#'
#' RES_RF <- calcRFLevels(D_MFAP4_cleaned, mod = RES_FLR$mod)
#' calcRFMeans(RES_RF)
calcRFMeans <- function(x) {
  avgRF <- NULL
  for (i in seq_along(x)) {
    avgRFCurrLevel <- mean(x[[i]])
    avgRF <- c(avgRF, avgRFCurrLevel)
  }
  names(avgRF) <- names(x)
  return(avgRF)
}


