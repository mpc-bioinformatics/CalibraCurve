#' Assemble result tables
#'
#' @param X **data.frame** \cr Original data set, e.g. result of
#'      \code{\link{readDataTable}}.
#' @param dataCleaned **list of data.frames** \cr Cleaned data, result of
#'      \code{\link{cleanData}}.
#' @param PLR_res **list** \cr Result object of \code{\link{calculate_PLR}}.
#' @param resFacDataV  **list** \cr  List of response factor values for each
#'      concentration level.
#' @param avgResFacDataV **numeric** \cr Vector of mean response factor values
#'      for the different concentration levels.
#' @param FLR_res **list** Result object of \code{\link{calculate_FLR}}.
#' @param mod **lm object** \cr Final linear model fit (object "mod" from
#'      results of \code{\link{calculate_FLR}}).
#' @param cvThres **numeric(1)** \cr Threshold for CV per concentration level
#'      in percent (default is 20).
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in
#'      percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in
#'      percent (default is 120).
#' @param substance **character(1)** \cr Name of the substance (default is
#'      "substance1").
#'
#' @returns List with the following elements:
#' - \code{result_table_conc_levels}: Result table with one line for each
#'      concentration level.
#' - \code{result_table_obs}: Result table with one line per observation (e.g.
#'      individual response factors for each data point).
#'
#' @importFrom checkmate assertNumeric assertCharacter
#' @importFrom stats aggregate predict
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
#' package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol",
#'         substColName = "Substance", assayNumber = 1)
#' data_cleaned <- cleanData(D_list[[1]])
#'
#' RES_PLR <- calculate_PLR(data_cleaned, calcContinuousPrelimRanges = FALSE)
#' RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)
#' RFs <- calcRF(data_cleaned, mod = RES_FLR$mod)
#'
#' assemble_results(
#'     X = D_list[[1]],
#'     dataCleaned = data_cleaned,
#'     PLR_res = RES_PLR,
#'     resFacDataV = RFs$RFs,
#'     avgResFacDataV = RFs$meanRFs,
#'     FLR_res = RES_FLR,
#'     mod = RES_FLR$mod
#' )
assemble_results <- function(X, dataCleaned,  PLR_res, resFacDataV,
    avgResFacDataV, FLR_res, mod, cvThres = 20, RfThresL = 80, RfThresU = 120,
    substance = "substance1") {
    checkmate::assertNumeric(RfThresL, lower = 0, upper = 100, finite = TRUE)
    checkmate::assertNumeric(RfThresU, lower = 100, finite = TRUE)
    checkmate::assertCharacter(substance, len = 1)

    concentrations <- as.numeric(sort(unique(X$Concentration)))
    concentrations_after_cleaning <- as.numeric(names(dataCleaned))

    # mean measurement for each concentration value
    mean_measurement <- stats::aggregate(X$Measurement,
                            by = list(Conc = X$Concentration),
                            FUN = mean)$x

    # predict measurements for each concentration using the final linear model
    estimated_measurement <- stats::predict(object = mod,
                        newdata = data.frame(Concentration = concentrations))

    ### thresholds for response factor
    RfThresUFactor <- RfThresU / 100
    RfThresLFactor <- RfThresL / 100
    concentrations_FLR <- as.numeric(rownames(FLR_res$perBiasAvgSDCV))

    # mean response factor for each concentration value (only within FLR)
    mean_RF <- mean(unlist(resFacDataV[concentrations %in% concentrations_FLR]))
    hLineLow <- mean_RF * RfThresLFactor
    hLineUpper <- mean_RF * RfThresUFactor

    result_table_conc_levels <- .makeResult_table_conc_levels(X, mod, PLR_res,
                                    FLR_res, avgResFacDataV, cvThres,
                                    hLineUpper, hLineLow, concentrations,
                                    concentrations_FLR,
                                    concentrations_after_cleaning, substance)

    result_table_obs <- .makeResult_table_obs(X, concentrations,
                            concentrations_FLR, concentrations_after_cleaning,
                            hLineUpper, hLineLow, FLR_res, resFacDataV,
                            substance)

    return(list(
        result_table_conc_levels = result_table_conc_levels,
        result_table_obs = result_table_obs
    ))
}





#' Save results of CalibraCurve
#'
#' @param CC_res **list** Result object of \code{\link{CalibraCurve}}.
#' @param output_path **character(1)** \cr Path to the output directory.
#' @param suffix **character(1)** \cr Suffix for the output files, ideally
#'      starting with "_" (default is "").
#'
#' @returns Returns nothing, but the function saves the results to the specified
#'      output path.
#'
#' @importFrom openxlsx write.xlsx
#'
.saveCCResult <- function(CC_res, output_path, suffix = "") {
    # save result tables
    openxlsx::write.xlsx(CC_res$result_table_conc_levels,
        file = paste0(output_path, "/result_table_conc_levels", suffix,
                                                                    ".xlsx"),
        rowNames = FALSE, keepNA = TRUE
    )
    openxlsx::write.xlsx(CC_res$result_table_obs,
        file = paste0(output_path, "/result_table_obs", suffix, ".xlsx"),
        rowNames = FALSE, keepNA = TRUE
    )
    # save whole result object
    saveRDS(CC_res, file = paste0(output_path, "/CC_res", suffix, ".rds"))
    return(invisible(NULL))
}











.makeResult_table_conc_levels <- function(X,mod, PLR_res, FLR_res,
        avgResFacDataV, cvThres, hLineUpper, hLineLow, concentrations,
        concentrations_FLR, concentrations_after_cleaning, substance) {
    mean_measurement <- stats::aggregate(X$Measurement,
                            by = list(Conc = X$Concentration), FUN = mean)$x
    estimated_measurement <- stats::predict(object = mod,
                        newdata = data.frame(Concentration = concentrations))

    # initialize the first result table (one line per concentration level)
    tab <- data.frame(substance = substance, concentration = concentrations,
        mean_measurement = mean_measurement,
        estimated_measurement = estimated_measurement,
        removed_while_cleaning = FALSE, CV = NA, CV_within_thres = NA,
        preliminary_linear_range = concentrations %in% PLR_res$prelimConcLevels,
        mean_percentage_bias = NA, SD_percentage_bias = NA,
        CV_percentage_bias = NA, mean_response_factor = NA,
        RF_within_thres = NA, final_linear_range = NA)

    ## result only for concentrations that were not removed during cleaning:
    for (i in seq_along(concentrations)) {
        ind <- which(concentrations[i] == concentrations_after_cleaning)
        if (length(ind) >= 1) {
            tab$removed_while_cleaning[i] <- FALSE
            tab$CV[i] <- PLR_res$concLevelsCV[ind]
            tab$mean_response_factor[i] <- avgResFacDataV[ind]
        } else {
            tab$removed_while_cleaning[i] <- TRUE
        }
    }
    tab$CV_within_thres <- tab$CV <= cvThres
    tab$RF_within_thres <- tab$mean_response_factor <= hLineUpper &
        tab$mean_response_factor >= hLineLow
    # fill table with percent bias information (only within FLR)
    for (i in seq_along(concentrations)) {
        ind <- which(concentrations[i] == concentrations_FLR)
        if (length(ind) >= 1) {
            tab$mean_percentage_bias[i] <-
                FLR_res$perBiasAvgSDCV$avgPerBias[ind]
            tab$SD_percentage_bias[i] <-
                FLR_res$perBiasAvgSDCV$stdDevPerBias[ind]
            tab$CV_percentage_bias[i] <- FLR_res$perBiasAvgSDCV$CV_PerBias[ind]
            tab$final_linear_range[i] <- TRUE
        } else {
            tab$final_linear_range[i] <- FALSE
        }
    }
    return(tab)
}







.makeResult_table_obs <- function(X, concentrations, concentrations_FLR,
                                conce_after_cleaning, hLineUpper,
                                hLineLow, FLR_res, resFacDataV, substance) {
    tab <- data.frame(
        substance = rep(substance, nrow(X)),
        concentration = X$Concentration,
        measurement = X$Measurement,
        removed_while_cleaning = !(X$Concentration %in% conce_after_cleaning) |
            is.na(X$Measurement) | X$Measurement == 0 | X$Concentration == 0,
        percentage_bias = NA,
        response_factor = NA,
        RF_within_thres = NA,
        final_linear_range = X$Concentration %in% concentrations_FLR
    )
    rownames(tab) <- NULL

    ## result only for concentrations that were not removed during cleaning:

    tab$response_factor[!tab$removed_while_cleaning] <- unlist(resFacDataV)
    tab$RF_within_thres <- tab$response_factor <= hLineUpper &
        tab$response_factor >= hLineLow

    # fill table with percent bias information (only within final linear range)
    perBias <- FLR_res$perBias
    for (i in seq_along(concentrations)) {
        if (concentrations[i] %in% concentrations_FLR) {
            ind1 <- which(concentrations_FLR == concentrations[i])
            ind2 <- which(tab$concentration == concentrations[i])
            tab$percentage_bias[ind2] <- perBias[[ind1]]
        }
    }
    return(tab)
}

