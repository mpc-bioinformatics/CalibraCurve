

#' Assemble result tables
#'
#' @param X **data.frame** \cr Original data set, result of \code{\link{readData}}.
#' @param dataCleaned **list of data.frames** \cr Cleaned data, result of \code{\link{cleanData}}.
#' @param cv_thres **numeric(1)** \cr Threshold for CV per concentration level in percent (default is 20).
#' @param PLR_res **list** \cr Result object of \code{\link{calculate_PLR}}.
#' @param resFacDataV  **list** \cr Result of \code{\link{calcRFLevels}}. List of response factor values for each
#' @param avgResFacDataV **numeric** \cr Result of \code{\link{calcRFMeans}}. Vector of mean response factor values for the different concentration levels.
#' @param FLR_res **list** Result object of \code{\link{calculate_FLR}}.
#' @param mod **lm object** \cr Final linear model fit (object "mod" from results of \code{\link{calculate_FLR}}).
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in percent (default is 120).
#' @param substance **character(1)** \cr Name of the substance (default is "substance1").
#'
#' @returns List with the following elements:
#' - \code{result_table_conc_levels}: Result table with one line for each concentration level.
#' - \code{result_table_obs}: Result table with one line per observation (e.g. individual response factors for each data point).
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
#' resFacDataV <- calcRFLevels(D_MFAP4_cleaned, mod = RES_FLR$mod)
#'
#' avgResFacDataV <- calcRFMeans(resFacDataV)
#'
#' assemble_results(X = D_MFAP4,
#'                  dataCleaned = D_MFAP4_cleaned,
#'                  PLR_res = RES_PLR,
#'                  resFacDataV = resFacDataV,
#'                  avgResFacDataV = avgResFacDataV,
#'                  FLR_res = RES_FLR,
#'                  mod = RES_FLR$mod)
#'
assemble_results <- function(X,
                             dataCleaned,
                             cv_thres = 20,
                             PLR_res,
                             resFacDataV,
                             avgResFacDataV,
                             FLR_res,
                             mod,
                             RfThresL = 80,
                             RfThresU = 120,
                             substance = "substance1"
                             ) {

  checkmate::assertNumeric(RfThresL, lower = 0, upper = 100, finite = TRUE)
  checkmate::assertNumeric(RfThresU, lower = 100)
  checkmate::assertCharacter(substance, len = 1)


  concentrations <- as.numeric(sort(unique(X$Concentration)))
  concentrations_after_cleaning <- as.numeric(names(dataCleaned))

  # mean measurement for each concentration value
  mean_measurement <- sapply(split(X, X$Concentration), function(x) mean(x$Measurement))

  #predict measurements for each concentration using the final linear model
  estimated_measurement <- stats::predict(object = mod, newdata = data.frame(Concentration = concentrations))

  ### thresholds for response factor
  RfThresUFactor <- RfThresU/100
  RfThresLFactor <- RfThresL/100
  concentrations_FLR <- as.numeric(rownames(FLR_res$perBiasAvgSDCV))

  # mean response factor for each concentration value (only within final linear range)
  mean_RF <- mean(unlist(resFacDataV[concentrations %in% concentrations_FLR]))
  hLineLow <- mean_RF * RfThresLFactor
  hLineUpper <- mean_RF * RfThresUFactor


  # initialize the first result table (one line per concentration level)
  result_table_conc_levels <- data.frame(
    substance = substance,
    concentration = concentrations,
    mean_measurement = mean_measurement,
    estimated_measurement = estimated_measurement,
    removed_while_cleaning = FALSE,
    CV = NA,
    CV_within_thres = NA,
    preliminary_linear_range = concentrations %in% PLR_res$prelimConcLevels,
    mean_percentage_bias = NA,
    SD_percentage_bias = NA,
    CV_percentage_bias = NA,
    mean_response_factor = NA,
    RF_within_thres = NA,
    final_linear_range = NA
  )

  ## result only for concentrations that were not removed during cleaning:
  for (i in seq_along(concentrations)) {
    ind <- which(concentrations[i] == concentrations_after_cleaning)
    if (length(ind) >= 1) {
      result_table_conc_levels$removed_while_cleaning[i] <- FALSE
      result_table_conc_levels$CV[i] <- PLR_res$concLevelsCV[ind]
      result_table_conc_levels$mean_response_factor[i] <- avgResFacDataV[ind]
    } else {
      result_table_conc_levels$removed_while_cleaning[i] <- TRUE
    }
  }

  result_table_conc_levels$CV_within_thres <- result_table_conc_levels$CV <= cv_thres
  result_table_conc_levels$RF_within_thres <- result_table_conc_levels$mean_response_factor <= hLineUpper & result_table_conc_levels$mean_response_factor >= hLineLow


  # fill table with percent bias information (only within final linear range) and
  # if the concentration level is within the final linear range
  for (i in seq_along(concentrations)) {
    ind <- which(concentrations[i] == concentrations_FLR)
    if (length(ind) >= 1) {
      result_table_conc_levels$mean_percentage_bias[i] <- FLR_res$perBiasAvgSDCV$avgPerBias[ind]
      result_table_conc_levels$SD_percentage_bias[i] <- FLR_res$perBiasAvgSDCV$stdDevPerBias[ind]
      result_table_conc_levels$CV_percentage_bias[i] <- FLR_res$perBiasAvgSDCV$CV_PerBias[ind]
      result_table_conc_levels$final_linear_range[i] <- TRUE
    } else {
      result_table_conc_levels$final_linear_range[i] <- FALSE
    }
  }




  # initialize second result table (one line per observation)
  result_table_obs <- data.frame(
    substance = rep(substance, nrow(X)),
    concentration = X$Concentration,
    measurement = X$Measurement,
    removed_while_cleaning = !(X$Concentration %in% concentrations_after_cleaning) |
      is.na(X$Measurement) | X$Measurement == 0 | X$Concentration == 0,
    percentage_bias = NA,
    response_factor = NA,
    RF_within_thres = NA,
    final_linear_range = X$Concentration %in% concentrations_FLR
  )
  rownames(result_table_obs) <- NULL

  ## result only for concentrations that were not removed during cleaning:

  resFacDataV2 <<- resFacDataV
  result_table_obs2 <<- result_table_obs

  result_table_obs$response_factor[!result_table_obs$removed_while_cleaning] <- unlist(resFacDataV)
  result_table_obs$RF_within_thres <- result_table_obs$response_factor <= hLineUpper & result_table_obs$response_factor >= hLineLow


  # fill table with percent bias information (only within final linear range)
  perBias <- FLR_res$perBias
  for (i in 1:length(concentrations)) {
    if (concentrations[i] %in% concentrations_FLR) {
      ind1 <- which(concentrations_FLR == concentrations[i])
      ind2 <- which(result_table_obs$concentration == concentrations[i])
      result_table_obs$percentage_bias[ind2] <- perBias[[ind1]]
    }
  }



  return(list(result_table_conc_levels = result_table_conc_levels,
         result_table_obs = result_table_obs))

}





#' Save results of CalibraCurve
#'
#' @param CC_res **list** Result object of \code{\link{CalibraCurve}}.
#' @param output_path **character(1)** \cr Path to the output directory.
#' @param suffix **character(1)** \cr Suffix for the output files, ideally starting with "_" (default is "").
#'
#' @returns Returns nothing, but the function saves the results to the specified output path.
#' @export
saveCCResult <- function(CC_res, output_path, suffix = "") {
  # save result tables
  openxlsx::write.xlsx(CC_res$result_table_conc_levels,
            file = paste0(output_path, "/result_table_conc_levels", suffix, ".xlsx"),
            rowNames = FALSE, keepNA = TRUE)
  openxlsx::write.xlsx(CC_res$result_table_obs,
            file = paste0(output_path, "/result_table_obs", suffix, ".xlsx"),
            rowNames = FALSE, keepNA = TRUE)

  # save whole result object
  saveRDS(CC_res, file = paste0(output_path, "/CC_res", suffix, ".rds"))
  return(invisible(NULL))
}









