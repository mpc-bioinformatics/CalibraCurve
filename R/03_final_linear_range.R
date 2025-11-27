
#' Calculate the final linear range
#'
#' @param dataPrelim **list of data.frames** \cr List of data.frames containing
#'      data only within the preliminary linear range (result from
#'      \code{\link{calculate_PLR}}).
#' @param weightingMethod **character(1)** \cr Method for weighting (currently
#'      "1/x", "1/x^2" and "None" are supported, default is 1/x^2).
#' @param centralTendencyMeasure **character(1)** \cr Method for calculating
#'      average percent bias, "mean" (default) or "median".
#' @param perBiasThres **numeric(1)** \cr Threshold for average percent bias in
#'      percent, default is 20.
#' @param considerPerBiasCV **logical(1)** \cr If TRUE, CV is considered for the
#'      elimination of the concentration level (default). CV will only be
#'      considered if the difference in percent bias values is lower than
#'      perBiasDistThres.
#' @param perBiasDistThres **numeric(1)** \cr Threshold for the difference in
#'      average percent bias in percent (for lower differences, CV will be
#'      considered), default is 10. Only relevant if considerPerBiasCV = TRUE.
#'
#' @returns List with the following elements:
#' - \code{dataFinal}: List of data.frames containing data only within the
#'      final linear range.
#' - \code{mod}: lm-object containing the final linear model (weighted or
#'      unweighted, depending on weightingMethod).
#' - \code{perBias}: result list of \code{\link{.calcPerBiasLevels}}.
#' - \code{perBiasAvgSDCV}: result data.frame of
#'      \code{\link{.calcPerBiasAvgSDCV}}.
#'
#' @importFrom checkmate assertChoice assertNumeric assertFlag
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
#'         package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol",
#'         substColName = "Substance", assayNumber = 1)
#' data_cleaned <- cleanData(D_list[[1]])
#' RES_PLR <- calculate_PLR(data_cleaned, calcContinuousPrelimRanges = FALSE)
#'
#' calculate_FLR(RES_PLR$dataPrelim)
calculate_FLR <- function(dataPrelim, weightingMethod = "1/x^2",
                        centralTendencyMeasure = "mean", perBiasThres = 20,
                        considerPerBiasCV = TRUE, perBiasDistThres = 10) {
    ### check input arguments
    checkmate::assertChoice(weightingMethod,
                            choices = c("1/x", "1/x^2", "None"))
    checkmate::assertChoice(centralTendencyMeasure,
                            choices = c("mean", "median"))
    checkmate::assertNumeric(perBiasThres, lower = 0, len = 1)
    checkmate::assertFlag(considerPerBiasCV)
    checkmate::assertNumeric(perBiasDistThres, lower = 0, len = 1)
    dataFinal <- dataPrelim
    finalRangeReached <- FALSE
    while (!finalRangeReached) {
        if (length(dataFinal) < 2) {
            stop("Only one concentration level left in linear range. Increasing
                perBiasThres may help, but makes results less accurate.")
        }
        ## calculate the weights for each concentration:
        if (weightingMethod != "None") {
            allWeights <- unlist(lapply(dataFinal, FUN = .calcWeights,
                                        weightingMethod = weightingMethod))
        } else {
            allWeights <- NULL
        }
        mod <- .calcLinearModel(dataFinal, weights = allWeights)
        perBias <- .calcPerBiasLevels(dataFinal, LMfit = mod)
        perBiasAvgSDCV <- .calcPerBiasAvgSDCV(perBias,
                                            method = centralTendencyMeasure)
        checkFLR <- .checkFinalRange(perBiasInfo = perBiasAvgSDCV,
                                    perBiasThres = perBiasThres)
        if (checkFLR) {
            finalRangeReached <- TRUE
        } else {
            # Final linear range is not reached yet, so it is decided if the
            # lowest or highest concentration level will be removed
            removeLow <- .selctConcLevel(perBiasAvgSDCV,
                                        perBiasT = perBiasThres,
                                        consPerBiasCV = considerPerBiasCV,
                                        perBiasDistT = perBiasDistThres)
            if (removeLow) {
                dataFinal[[1]] <- NULL # removes the concentration level
            } else {
                dataFinal[[length(dataFinal)]] <- NULL
            }
        }
    }
    return(list(dataFinal = dataFinal, mod = mod, perBias = perBias,
                perBiasAvgSDCV = perBiasAvgSDCV))
}




#' FLR: calculate weights for linear model for a specific concentration level
#'
#' @param x **data.frame** \cr Data.frame containing data for a specific
#'      concentration level.
#' @param weightingMethod **character(1)** \cr Method for weighting (currently
#'      "1/x", "1/x^2" and "None" are supported, default is "1/x^2").
#'
#' @returns Vector of a constant weights for each measurement in this
#'      concentration level.
.calcWeights <- function(x, weightingMethod = "1/x^2") {
    currConcLevData <- x
    if (weightingMethod == "1/x^2") {
        weightsCurrLev <- 1 / (currConcLevData$Concentration^2)
    }
    if (weightingMethod == "1/x") {
        weightsCurrLev <- 1 / currConcLevData$Concentration
    }
    result <- weightsCurrLev
    return(result)
}



#' FLR: calculate weighted or unweighted linear model
#'
#' @param x **list of data.frames** \cr List of data frames (one entry for each
#'      concentration level), e.g. output "dataPrelim" from
#'      \code{\link{calculate_PLR}}.
#' @param weights **numeric** \cr Vector of weights as calculated by applying
#'      \code{\link{.calcWeights}} (default is NULL and will result in an
#'      unweighted model).
#'
#' @returns Fit of the linear model as an object of class "lm".
#'
#' @importFrom stats lm
#'
.calcLinearModel <- function(x, weights = NULL) {
    ## combine list elements to a data set
    dataSetDF <- do.call(rbind, x)
    lmfit <- stats::lm(Measurement ~ Concentration, data = dataSetDF,
                        weights = weights)
    return(lmfit)
}


#' FLR: calculate percent bias for a specific concentration level
#'
#' @param x **data.frame** \cr Data.frame containing data for a specific
#'      concentration level.
#' @param LMfit **lm object** \cr Linear model fit as calculated by
#'      \code{\link{.calcLinearModel}}.
#' @param expConc **numeric(1)** \cr Expected (known) concentration level.
#'
#' @returns Vector of percent bias values for each data point in this
#'      concentration level.
.calcPerBias <- function(x, LMfit, expConc) {
    x <- x$Measurement
    coeff <- LMfit$coefficients
    calculatedConc <- (x - coeff[1]) / coeff[2]
    perBias <- 100 * abs(calculatedConc - expConc) / expConc
    names(perBias) <- NULL

    return(perBias)
}



#' FLR: calculate list of percent bias values for all concentration levels
#'
#' @param x **list of data.frames** \cr List of data frames (one entry for each
#'      concentration level), e.g. output "dataPrelim" from
#'      \code{\link{calculate_PLR}}.
#' @param LMfit **lm object** \cr Linear model fit as calculated by
#'      \code{\link{.calcLinearModel}}.
#'
#' @returns Matrix with percent bias values for each data point per
#'      concentration level
.calcPerBiasLevels <- function(x, LMfit) {
    concentrations <- as.numeric(names(x))
    perBiasValues <- mapply(FUN = .calcPerBias, x = x, expConc = concentrations,
                            MoreArgs = list(LMfit = LMfit))
    if(!is.matrix(perBiasValues)) perBiasValues <-  t(as.matrix(perBiasValues))
    return(perBiasValues)
}


#' FLR: calculate average, SD and CV percent bias for each concentration level
#'
#' @param x **matrix** \cr Result of
#'      \code{\link{.calcPerBiasLevels}}.
#' @param method **character(1)** \cr Method for calculating the average
#'      percent bias: "mean" (default) or "median".
#'
#' @returns data frame with 3 columns: avgPerBias, stdDevPerBias, CV_PerBias
#'          \cr each row is one concentration level
#'
#' @importFrom stats sd median
#'
.calcPerBiasAvgSDCV <- function(x, method = "mean") {

    if (method == "mean") {
        avgPerBias <- colMeans(x)
    }
    if (method == "median") {
        avgPerBias <- apply(x, 2, stats::median)
    }
    stdDevPerBias <- apply(x, 2, stats::sd)
    CV_PerBias <- stdDevPerBias / avgPerBias * 100

    result <- data.frame(avgPerBias, stdDevPerBias, CV_PerBias)
    rownames(result) <- colnames(x)
    return(result)

}


#' FLR: checks if final linear range has been reached (compare average percent
#'      bias with threshold)
#'
#' @param perBiasInfo **data.frame** Result of
#'      \code{\link{.calcPerBiasAvgSDCV}}.
#' @param perBiasThres **numeric(1)** **numeric(1)** \cr Threshold for average
#'      percent bias in percent, default is 20.
#'
#' @returns  TRUE if both lowest and highest concentration level passed the
#'      check (and the final linear range is reached), FALSE otherwise
.checkFinalRange <- function(perBiasInfo, perBiasThres = 20) {
    bothLevelsPassed <- FALSE
    lowLevelPassed <- FALSE
    highLevelPassed <- FALSE

    if (perBiasInfo$avgPerBias[1] <= perBiasThres) {
        lowLevelPassed <- TRUE
    }
    if (perBiasInfo$avgPerBias[length(perBiasInfo$avgPerBias)] <= perBiasThres){
        highLevelPassed <- TRUE
    }

    if (lowLevelPassed && highLevelPassed) {
        bothLevelsPassed <- TRUE
    }
    result <- bothLevelsPassed
    return(result)
}


#' FLR: selects the highest or lowest concentration level for removal
#'
#' @param x **data.frame** Result of \code{\link{.calcPerBiasAvgSDCV}}.
#' @param perBiasT **numeric(1)** \cr Threshold for average percent bias in
#'      percent, default is 20.
#' @param consPerBiasCV consider CV? default is TRUE
#' @param perBiasDistT **numeric(1)** \cr Threshold for the difference in
#'      average percent bias in percent (for lower differences, CV will be
#'      considered), default is 10. Only used if consPerBiasCV is TRUE.
#'
#' @returns TRUE, if lowest concentration will be removed, FALSE if highest
#'      will be removed
.selctConcLevel <- function(x, perBiasT = 20, consPerBiasCV = TRUE,
                            perBiasDistT = 10) {
    removeLow <- NULL
    featuresLowestLevel <- x[1, ]
    featuresHighestLevel <- x[nrow(x), ]

    lowPerBias <- featuresLowestLevel$avgPerBias
    highPerBias <- featuresHighestLevel$avgPerBias

    lowFails <- (lowPerBias >= perBiasT)
    highFails <- (highPerBias >= perBiasT)

    if (lowFails && highFails) {
        # Case: both the highest and lowest concentration level fails percent
        # bias threshold

        if (!consPerBiasCV) {
            # remove the level with the higher percent bias
            removeLow <- ifelse(lowPerBias >= highPerBias, TRUE, FALSE)
        } else {
            # calculate distance between the two percent bias values
            dist <- abs(lowPerBias - highPerBias)
            if (dist <= perBiasDistT) {
                # CV values are considered because distance between the two
                # percent bias values is low
                removeLow <-
                    ifelse(featuresLowestLevel$CV >= featuresHighestLevel$CV,
                                    TRUE, FALSE)
            } else {
                # only percent bias values are used because distance is high
                removeLow <- ifelse(lowPerBias >= highPerBias, TRUE, FALSE)
            }
        }
    } else {
        # Case: only one levels fails the threshold criteria
        removeLow <- ifelse(lowFails, TRUE, FALSE)
    }
    return(removeLow)
}


