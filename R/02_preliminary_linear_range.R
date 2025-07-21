#' PLR: calculate CV for one concentration level
#'
#' @param x **data.frame** \cr Data.frame containing data for a specific concentration level.
#'
#' @returns **numeric(1)** \cr Coefficient of variation (CV) for the given concentration level.
calcCV <- function(x) {
    SD <- stats::sd(x$Measurement)
    Mean <- mean(x$Measurement)
    CV <- SD / Mean * 100
    names(CV) <- NULL # x$Concentration[1]
    return(CV)
}


#' PLR: Calculate key features for the existing continuous preliminary ranges
#'
#' @param index **integer** \cr Integer vector containing indices of concentration levels with CV < threshold.
#'
#' @returns Data.frame with one range in each row. Contains the columns "startPoints", "endPoints", and "extent".
calcContPrelimRanges <- function(index) {
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
#' @param dataCleaned **data.frame** \cr Data cleaned by \code{\link{cleanData}}.
#' @param cv_thres **numeric(1)** \cr Threshold for CV per concentration level in percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest continuous range is selected (default is TRUE).
#'                                                      If FALSE, gaps with CVs larger than the threshold may be included.
#'
#' @returns List with the following elements:
#' - \code{dataPrelim}: List of data.frames containing data only within the preliminary linear range.
#' - \code{concLevelsCV}: Vector with the calculated CV for each concentration level.
#' - \code{prelimConcLevels}: Vector with the concentration levels within the preliminary linear range.
#'
#' @export
#'
#' @examples
#' data(D_MFAP4)
#' D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
#' calculate_PLR(D_MFAP4_cleaned,
#'     cv_thres = 10,
#'     calcContinuousPrelimRanges = TRUE
#' )
#' ## With calcContinuousPrelimRanges = TRUE, only concentration levels 0.025-2.5
#' ## are part of the preliminary linear range.
#'
#' calculate_PLR(D_MFAP4_cleaned,
#'     cv_thres = 10,
#'     calcContinuousPrelimRanges = FALSE
#' )
#' ## With calcContinuousPrelimRanges = FALSE, gaps are allowed, so the preliminary
#' ## linear range is 0.025-25 (although 5 has a CV larger than the threshold 10).
#'
calculate_PLR <- function(dataCleaned,
    cv_thres = 20,
    calcContinuousPrelimRanges = TRUE) {
    ### check input arguments
    checkmate::assert_numeric(cv_thres, len = 1)
    checkmate::assert_flag(calcContinuousPrelimRanges)

    ### calculate CV for each concentration level
    concLevelsCV <- sapply(dataCleaned, calcCV)
    ## which concentration levels have a CV lower than the threshold?
    index <- which(concLevelsCV <= cv_thres)

    if (length(index) <= 1) {
        stop(paste0("No preliminary linear range with CV <= ", cv_thres, " could be calculated. Please check your data or increase the CV threshold."))
    }

    ### calculate candidates for the PLR (parts where CV < threshold)
    rangeProperties <- calcContPrelimRanges(index)
    ### minimal and maximal concentration with CV < threshold.
    indexStart <- min(rangeProperties[, 1])
    indexEnd <- max(rangeProperties[, 2])

    if (!calcContinuousPrelimRanges) { # FALSE: Selecting all levels between the lowest and highest concentration level, where CV < threshold
        dataPrelim <- dataCleaned[indexStart:indexEnd]
    } else { # TRUE: Computing the longest continuous preliminary range
        longestRange <- rangeProperties[which(rangeProperties$extent == max(rangeProperties$extent)), ]

        if (nrow(longestRange) > 1) { # Special case: more than one range with the same number of concentration levels that pass the CV threshold exist
            warning(paste0("More than one preliminary linear ranges with the same number of concentration levels passing the CV threshold exist. \n The first one will be selected. Set calcContinuousPrelimRanges to TRUE if you want to allow gaps with CV >", cv_thres, "."))
            dataPrelim <- dataCleaned[longestRange$startPoints[1]:longestRange$endPoints[1]]
        } else {
            if (longestRange$extent == 0) {
                stop(paste0("No preliminary linear range with CV <= ", cv_thres, " could be calculated. Please check your data or increase the CV threshold."))
            }
            dataPrelim <- dataCleaned[longestRange$startPoints:longestRange$endPoints] # if nrow(longestRange) == 1
        }
    }

    ### get concentration levels of the preliminary linear range and use them as names of the new, filtered data list
    prelimConcLevels <- sapply(dataPrelim, FUN = function(x) x$Concentration[1])
    names(prelimConcLevels) <- NULL
    names(dataPrelim) <- prelimConcLevels

    ### get CV values of the preliminary linear range
    ## concLevelsCVSel <- concLevelsCV[as.character(prelimConcLevels)]

    # Cancel calculations if less than two remaining concentration levels exist, which pass the CV value threshold
    if (length(prelimConcLevels) < 2) {
        stop(paste0("Less than two concentration levels with CV <= ", cv_thres, " exist. \n Please check your data or increase the CV threshold."))
    }

    return(list(
        dataPrelim = dataPrelim,
        concLevelsCV = concLevelsCV,
        prelimConcLevels = prelimConcLevels
    )) ### TODO: could be extracted from dataPrelim. Is this needed?
    ## concLevelsCVSel = concLevelsCVSel))
}
