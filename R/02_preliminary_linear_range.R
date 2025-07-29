#' Calculate preliminary linear range (PLR)
#'
#' @param dataCleaned **data.frame** \cr Data cleaned by \code{\link{cleanData}}.
#' @param cvThres **numeric(1)** \cr Threshold for CV per concentration level in
#'      percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest
#'      continuous range fulfilling the CV threshold is selected (default is
#'      TRUE). If FALSE, gaps with CVs larger than the threshold may be included.
#'
#' @returns List with the following elements:
#' - \code{dataPrelim}: List of data.frames containing data only within the
#'      preliminary linear range.
#' - \code{concLevelsCV}: Vector with the calculated CV for each concentration
#'      level.
#' - \code{prelimConcLevels}: Vector with the concentration levels within the
#'      preliminary linear range.
#'
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#' package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
#'                      assayNumber = 1)
#' data_cleaned <- cleanData(D_list[[1]])
#'
#' calculate_PLR(data_cleaned, calcContinuousPrelimRanges = FALSE)
calculate_PLR <- function(dataCleaned, cvThres = 20,
                          calcContinuousPrelimRanges = TRUE) {
    checkmate::assert_numeric(cvThres, len = 1, lower = 0)
    checkmate::assert_flag(calcContinuousPrelimRanges)

    ### calculate CV for each concentration level
    concLevelsCV <- vapply(dataCleaned, .calcCV, numeric(1))

    ## which concentration levels have a CV lower than the threshold?
    index <- which(concLevelsCV <= cvThres)

    if (length(index) <= 1) {
        stop("No preliminary linear range with CV <= ", cvThres, " could be
             calculated. Please check your data or increase the CV threshold.")
    }

    ### calculate candidates for the PLR (parts where CV < threshold)
    ranges <- .calcContPrelimRanges(index)
    ### minimal and maximal concentration with CV < threshold.
    indexStart <- min(ranges[, 1])
    indexEnd <- max(ranges[, 2])

    if (!calcContinuousPrelimRanges) { # FALSE: Selecting all levels between the lowest and highest concentration level, where CV < threshold
        dataPrelim <- dataCleaned[indexStart:indexEnd]
    } else { # TRUE: Computing the longest continuous preliminary range
        longestRange <- ranges[which(ranges$extent == max(ranges$extent)), ]

        if (nrow(longestRange) > 1) { # Special case: more than one range with the same number of concentration levels that pass the CV threshold exist
            warning("Multiple preliminary linear ranges with the same length
                    exist. The first one will be selected. Set
                    calcContinuousPrelimRanges to FALSE if you want to allow gaps
                    with CV >", cvThres, ".")
            dataPrelim <- dataCleaned[longestRange$start[1]:longestRange$end[1]]
        } else {
            if (longestRange$extent == 0) {
                stop("No preliminary linear range with CV <= ", cvThres,
                     " could be calculated. Please check your data or increase
                     the CV threshold.")
            }
            dataPrelim <- dataCleaned[longestRange$start:longestRange$end]
        }
    }
    ### get concentration levels of the preliminary linear range and use them as names of the new, filtered data list
    prelimConcLevels <- sapply(dataPrelim, FUN = function(x) x$Concentration[1])
    names(prelimConcLevels) <- NULL
    names(dataPrelim) <- prelimConcLevels

    # Cancel calculations if less than two remaining concentration levels exist, which pass the CV value threshold
    if (length(prelimConcLevels) < 2) {
        stop("Less than two concentration levels with CV <= ", cvThres, " exist.
             Please check your data or increase the CV threshold.")
    }

    return(list(dataPrelim = dataPrelim, concLevelsCV = concLevelsCV,
        prelimConcLevels = prelimConcLevels))
}






#' PLR: calculate CV for one concentration level
#'
#' @param x **data.frame** \cr Data.frame containing data for a specific
#'      concentration level.
#'
#' @returns **numeric(1)** \cr Coefficient of variation (CV) for the given
#'      concentration level.
.calcCV <- function(x) {
    SD <- stats::sd(x$Measurement)
    Mean <- mean(x$Measurement)
    CV <- SD / Mean * 100
    names(CV) <- NULL
    return(CV)
}


#' PLR: Calculate key features for the existing continuous preliminary ranges
#'
#' @param index **integer** \cr Integer vector containing indices of
#'      concentration levels with CV < threshold.
#'
#' @returns Data.frame with one range in each row. Contains the columns
#'      "startPoints", "endPoints", and "extent".
.calcContPrelimRanges <- function(index) {
    # Calculating the end points for all ranges
    y <- index
    z <- index
    y <- y[-1] # remove first element
    z <- z[-length(z)] # remove last element
    differences <- y - z
    end <- which(differences > 1) # difference >1 means that there is a gap

    # Calculating the start points for all ranges
    start <- end + 1

    # Add special case:
    # the highest end point
    end <- c(end, length(index))
    # the lowest start point
    start <- c(1, start)

    # Calculating the extent of the ranges
    extent <- end - start + 1

    start <- index[start]
    end <- index[end]

    ranges <- as.data.frame(cbind(start, end, extent))
    return(ranges)
}


