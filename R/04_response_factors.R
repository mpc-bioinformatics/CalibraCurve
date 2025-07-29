

#' Calculate Response factors
#'
#' @description
#' Final linear range: Function, which returns a list with response factor
#'      values for a data set (given as list)
#'
#' @param x **list of data.frames** \cr List of data.frames containing data for
#'      each concentration level (result from \code{\link{cleanData}}).
#' @param mod **lm object** \cr Final linear model fit (object "mod" from
#'      results of \code{\link{calculate_FLR}}).
#'
#' @returns List with the following elements:
#' - \code{RFs}: List with response factors for each concentration level.
#' - \code{meanRFs}: Vector with mean response factors for each concentration
#'      level.
#'
#'
#' response factor values for each concentration level.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#' package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol",
#'         substColName = "Substance", assayNumber = 1)
#' data_cleaned <- cleanData(D_list[[1]])
#' RES_PLR <- calculate_PLR(data_cleaned, calcContinuousPrelimRanges = FALSE)
#' RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)
#'
#' calcRF(data_cleaned, mod = RES_FLR$mod)
calcRF <- function(x, mod) {
    interc <- unname(mod$coefficients[1])
    concentrations <- as.numeric(names(x))

    RFs <- mapply(FUN = .calcResponseFactors, x = x, expConc = concentrations,
           MoreArgs = list(intercept = interc), SIMPLIFY = FALSE)
    names(RFs) <- concentrations

    meanRFs <- vapply(RFs, mean, numeric(1))
    return(list(RFs = RFs, meanRFs = meanRFs))
}


#' Calculate Response factors
#'
#' @description
#' Function, which calculates a response factor for a single data point
#'
#'
#' @details
#' Formula obtained from:  Green, J. M., A practical guide to analytical method
#'      validation. Analytical Chemistry 1996, 68, 305A-309A.
#'
#'
#' @param x *data.frame** \cr Data.frame containing data for a specific
#'      concentration level.
#' @param intercept  **numeric(1)** \cr Intercept of the linear model.
#' @param expConc **numeric(1)** \cr  Expected concentration (known
#'      concentration value).
#'
#' @returns vector of response factors for this specific concentration level
.calcResponseFactors <- function(x, intercept, expConc) {
    result <- (x$Measurement - intercept) / expConc
    return(result)
}


