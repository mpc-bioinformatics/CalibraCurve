

#' Title
#'
#' @param X
#' @param dataValidated
#' @param cv_thres
#' @param PLR_res
#' @param avgResFacDataV
#' @param FLR_res
#' @param mod
#'
#' @returns
#' @export
#'
#' @examples
assemble_results <- function(X,
                             dataValidated,
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


  concentrations <- as.numeric(names(dataValidated))
  mean_measurement <- sapply(dataValidated, function(x) mean(x$Measurement))
  estimated_measurement <- predict(object = mod, newdata = data.frame(Concentration = concentrations))

  ### thresholds for response factor
  RfThresUFactor <- RfThresU/100
  RfThresLFactor <- RfThresL/100
  concentrations_FLR <- as.numeric(rownames(FLR_res$perBiasAvgSDCV))

  mean_RF <- mean(unlist(resFacDataV[concentrations %in% concentrations_FLR])) # mean RF (only for final linear range)
  hLineLow <- mean_RF * RfThresLFactor
  hLineUpper <- mean_RF * RfThresUFactor

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


  ##############################################

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







