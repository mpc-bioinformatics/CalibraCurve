

#' Assemle result tables
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

  concentrations <- as.numeric(names(dataCleaned))

  # mean measurement for each concentration value
  mean_measurement <- sapply(dataCleaned, function(x) mean(x$Measurement))

  #precdict measurements for each concentration using the final linear model
  estimated_measurement <- predict(object = mod, newdata = data.frame(Concentration = concentrations))

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
    CV = PLR_res$concLevelsCV,
    CV_within_thres = PLR_res$concLevelsCV < cv_thres,
    preliminary_linear_range = concentrations %in% PLR_res$prelimConcLevels,
    mean_percentage_bias = NA,
    SD_percentage_bias = NA,
    CV_percentage_bias = NA,
    mean_response_factor = avgResFacDataV,
    RF_within_thres = avgResFacDataV <= hLineUpper & avgResFacDataV >= hLineLow,
    final_linear_range = NA
  )

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
    substance = substance,
    concentration = X$Concentration,
    measurement = X$Measurement,
    percentage_bias = NA,
    response_factor = unlist(resFacDataV),
    RF_within_thres = unlist(resFacDataV) <= hLineUpper & unlist(resFacDataV) >= hLineLow,
    final_linear_range = X$Concentration %in% concentrations_FLR
  )
  rownames(result_table_obs) <- NULL

  # fill table with percent bias information (only within final linear range)
  perBias <- FLR_res$perBias
  for (i in 1:length(concentrations)) {
    if (concentrations[i] %in% concentrations_FLR) {
      print(i)
      ind1 <- which(concentrations_FLR == concentrations[i])
      ind2 <- which(result_table_obs$concentration == concentrations[i])
      result_table_obs$percentage_bias[ind2] <- perBias[[ind1]]

      print(ind1)
      print(ind2)
    }
  }



  return(list(result_table_conc_levels = result_table_conc_levels,
         result_table_obs = result_table_obs))

}







