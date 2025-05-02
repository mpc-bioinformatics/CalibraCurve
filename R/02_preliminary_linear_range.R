

#' PLR: calculate CV for one concentration level
#'
#' @param x **data.frame** \cr Data.frame containing data for a specific concentration level.
#'
#' @returns **numeric(1)** \cr Coefficient of variation (CV) for the given concentration level.
#'
#' @examples
calcCV <- function(x) {
  SD <-  sd(x$Measurement)
  Mean <- mean(x$Measurement)
  CV <- SD / Mean * 100
  names(CV) <- x$Concentration[1]
  return(CV)
}


#' PLR: Calculate key features for the existing continuous preliminary ranges
#'
#' @param index **integer** \cr Integer vector containing indices of concentration levels with CV < threshold.
#'
#' @returns Data.frame with one range in each row. Contains the columns "startPoints", "endPoints", and "extent".
#'
#' @examples
calcContPrelimRanges  <- function(index) {
  # Calculating the end points for all ranges
  y <- index
  z <- index
  y <- y[-1] # remove first element
  z <- z[-length(z)] # remove last element
  differences <- y - z
  endPoints <- which(differences > 1) # difference >1 means that there is a gap

  # Calculating the start points for all ranges
  startPoints <- endPoints + 1

  # Add special case:
  # the highest end point
  endPoints <- c(endPoints, length(index))
  # the lowest start point
  startPoints <- c(1, startPoints)

  # Calculating the extent of the ranges
  extent <- endPoints - startPoints + 1

  startPoints <- index[startPoints]
  endPoints <- index[endPoints]

  rangeProperties <- as.data.frame(cbind(startPoints, endPoints, extent))
  return(rangeProperties)
}





#' Calculate preliminary linear range (PLR)
#'
#' @param dataCleaned **data.frame** \cr Data cleaned by \code{\link{cleadData}}.
#' @param cv_thres **numeric(1)** \cr Threshold for CV per concentration level in percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest continuous range is selected (default is TRUE).
#'                                                      If FALSE, gaps with CVs larger than the threshold may be included.
#'
#' @returns List with the following elements:
#' - \code{dataPrelim}: Data frame with the data for the preliminary linear range.
#'
#' @export
#'
#' @examples
calculate_PLR <- function(dataCleaned,
                          cv_thres = 20,
                          calcContinuousPrelimRanges = TRUE) {

  ### check input arguments
  checkmate::assert_numeric(cv_thres, len = 1)
  checkmate::assert_flag(calcContinuousPrelimRanges)

  ### calculate CV for each concentration level
  concLevelsCV <- sapply(dataValidated, calcCV)
  ## which concentration levels have a CV lower than the threshold?
  index <- which(concLevelsCV <= cv_thres)


  ### calculate candidates for the PLR (parts where CV < threshold)
  rangeProperties <- calcContPrelimRanges(index)
  ### minimal and maximal concentration with CV < threshold.
  indexStart <- min(rangeProperties[,1])
  indexEnd <- max(rangeProperties[,2])



  if (!calcContinuousPrelimRanges) {# FALSE: Selecting all levels between the lowest and highest concentration level, where CV < threshold
    dataPrelim <- dataValidated[indexStart:indexEnd]
  } else { # TRUE: Computing the longest continuous preliminary range
    longestRange <- rangeProperties[which(rangeProperties$extent == max(rangeProperties$extent)),]

    if (nrow(longestRange) > 1) {# Special case: more than one range with the same number of concentration levels that pass the CV threshold exist
      warning("More than one preliminary linear range with the same number of concentration levels passing the CV threshold exist. \n The first one will be selected. Set calcContinuousPrelimRanges to TRUE if you want to allow gaps with CV > threshold.")
      dataPrelim <- dataValidated[longestRange$startPoints[1]:longestRange$endPoints[1]]
      ### TODO: does that make any sense?
    } else {
      if (longestRange$extent == 0) {
        stop("No preliminary linear range with (CV < threshold) could be calculated. Please check your data or increase the CV threshold.")
      }
      dataPrelim <- dataValidated[longestRange$startPoints:longestRange$endPoints]  # if nrow(longestRange) == 1
    }
  }

  ### get concentration levels of the preliminary linear range and use them as names of the new, filtered data list
  prelimConcLevels <- sapply(dataPrelim, FUN = function(x) x$Concentration[1])
  names(dataPrelim) <- prelimConcLevels

  ### get CV values of the preliminary linear range
  ##concLevelsCVSel <- concLevelsCV[as.character(prelimConcLevels)]

  # Cancel calculations if less than two remaining concentration levels exist, which pass the CV value threshold
  if (length(prelimConcLevels) < 2) {
    stop("Less than two concentration levels with CV < threshold exist. \n Please check your data or increase the CV threshold.")
  }

  return(list(dataPrelim = dataPrelim,
              concLevelsCV = concLevelsCV,
              prelimConcLevels = prelimConcLevels))
              ##concLevelsCVSel = concLevelsCVSel))
}



