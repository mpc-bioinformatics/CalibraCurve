


#' PLR: calculate CV for each concentration level
#'
#' @param x data frame
#'
#' @returns
#' @export
#'
#' @examples
calcCV <- function(x){
  actSD <-  sd(x$Measurement_1)
  actMean <- mean(x$Measurement_1)
  CV <- actSD/actMean*100
  names(CV) <- unique(x$Concentration)
  result <- CV
  return(result)
}




#' PLR: Calculate key features for the existing continuous preliminary ranges
#'
#' @param x data frame
#'
#' @returns
#' @export
#'
#' @examples
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
  return(result)
}





#' Calulate preliminary linear range (PLR)
#'
#' @param dataValidated data frame cleaned by cleanData()
#' @param cv_thres threshold for CV per concentration level (default is 20)
#' @param calcContinuousPrelimRanges
#'
#' @returns list
#' @export
#'
#' @examples
calculate_PLR <- function(dataValidated, cv_thres = 20, calcContinuousPrelimRanges = TRUE) {

  ### calculate CV for each concentration level
  concLevelsCV <- sapply(dataValidated, calcCV)
  #   0.00025     5e-04    0.0025     0.005     0.025      0.05      0.25       0.5       2.5         5        25
  # 46.343154 52.296684 33.771666 19.763155  8.664128  7.106430  6.271753  6.573934  8.358188 13.292540  7.145171

  ## which concentration levels have a CV lower than the threshold?
  index <- which(concLevelsCV <= cv_thres)
  # 0.005 0.025  0.05  0.25   0.5   2.5     5    25
  #     4     5     6     7     8     9    10    11

  propertiesOfRanges <- calcContPrelimRanges(index)
  indexStart <- min(propertiesOfRanges[,1])
  indexEnd <- max(propertiesOfRanges[,2])
  #### Dies berechnet den range zwischen der kleinsten Konzentration mit CV < 20% und der höchsten Konzentration mit CV < 20%

  calcContinuousPrelimRanges = TRUE
  #### User kann einstellen ob "Lücken" mit höherem CV ok sind oder nicht.
  #### wenn nicht, wird der längere Abschnitt gewählt.
  #### Falls die Abschnitte gleich lang sein sollten????
  if (!calcContinuousPrelimRanges) {# FALSE: Selecting all levels between the lowest and highest concentration level, where CV <threshold
    dataPrelim <- dataValidated[indexStart:indexEnd]
  } else { # Computing the longest continuous preliminary range
    longestRange <- propertiesOfRanges[which(propertiesOfRanges$extent == max(propertiesOfRanges$extent)),]

    if (nrow(longestRange) > 1) {# Special case: more than one range with the same number of concentration levels that pass the CV threshold exist
      dataPrelim <- dataValidated[indexStart:indexEnd]
      ### TODO: does that make any sense?
    } else {
      if (longestRange$extent == 0) {
        ### TODO: warning that no range with CV <20% was found
        next # Calculation stops for this sample
      }
      dataPrelim <- dataValidated[longestRange$startPoints:longestRange$endPoints]
    }
  }

  prelimConcLevels <- sapply(dataPrelim, FUN = function(x) x$Concentration[1])
  names(dataPrelim) <- prelimConcLevels
  concLevelsCVSel <- concLevelsCV[as.character(prelimConcLevels)]
  # Cancel calculations if less than two remaining concentration levels exist, which pass the CV value threshold
  if (length(concLevelsCVSel) < 2) {
    ### TODO: error message (not enough concentration levels with CV < threshold to calculate further)
    next # Calculation stops for this sample
  }

  return(list(dataPrelim = dataPrelim, concLevelsCV = concLevelsCV, prelimConcLevels = prelimConcLevels, concLevelsCVSel = concLevelsCVSel))
}



