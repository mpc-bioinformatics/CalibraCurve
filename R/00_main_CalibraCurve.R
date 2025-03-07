

#' CalibraCurve
#'
#' @param path
#' @param conc_col
#' @param meas_col
#' @param filetype
#' @param min_replicates
#' @param weightingMethod
#' @param centralTendencyMeasure
#' @param perBiasThres
#' @param considerPerBiasCV
#' @param perBiasDistThres
#'
#' @returns
#' @export
#'
#' @examples
CalibraCurve <- function(path,
                         conc_col,
                         meas_col,

                         filetype = "xlsx",
                         sep = ",",
                         dec = ";",
                         header = TRUE,
                         na.strings = c("NA", "NaN", "Filtered", "#NV"),
                         sheet = 1,

                         min_replicates = 3,

                         cv_thres = 20,
                         calcContinuousPrelimRanges = TRUE,


                         weightingMethod = "1/x^2",
                         centralTendencyMeasure = "mean",
                         #finalRangeCalculationMethod = "weighted_linear_model",
                         perBiasThres = 20,
                         considerPerBiasCV = TRUE,
                         perBiasDistThres = 10,

                         RfThresL = 80,
                         RfThresU = 120,

                         substance = "substance1"
                         ) {

  ## read in data
  X <- CalibraCurve::readData(path = path,
                              conc_col = conc_col,
                              meas_col = meas_col,
                              filetype = filetype,
                              sep = sep,
                              dec = dec,
                              header = header,
                              na.strings = na.strings,
                              sheet = sheet)

  ## clean data
  dataValidated <- CalibraCurve::cleanData(X,
                                           min_replicates = min_replicates)

  ## calculate preliminary linear range
  PLR_res <- CalibraCurve::calculate_PLR(dataValidated = dataValidated,
                                         cv_thres = cv_thres,
                                         calcContinuousPrelimRanges = calcContinuousPrelimRanges)


  ## calculate final linear range
  FLR_res <- calculate_FLR(PLR_res$dataPrelim,
                       weightingMethod = weightingMethod,
                       centralTendencyMeasure = centralTendencyMeasure,
                       #finalRangeCalculationMethod = "weighted_linear_model",
                       perBiasThres = perBiasThres,
                       considerPerBiasCV = considerPerBiasCV,
                       perBiasDistThres = perBiasDistThres)

  dataFinal <- FLR_res$dataFinal
  mod <- FLR_res$mod


  ### calculate response factors
  #resFacDataF <- calcRFLevels(dataFinal, mod = lmWeighted) ### TODO: wird das überhaupt gebraucht?
  resFacDataV <- calcRFLevels(dataValidated, mod = mod)

  # Calculation of mean response factor values
  #avgResFacDataF <- calcRFMeans(resFacDataF) ### TODO: wird das überhaupt gebraucht?
  avgResFacDataV <- calcRFMeans(resFacDataV)


  #### generate result tables
  tables <- assemble_results(X = X,
                             dataValidated = dataValidated,
                             cv_thres = cv_thres,
                             PLR_res = PLR_res,
                             resFacDataV = resFacDataV,
                             avgResFacDataV = avgResFacDataV,
                             FLR_res = FLR_res,
                             mod = mod,
                             RfThresL = RfThresL,
                             RfThresU = RfThresU
  )



  return(list(mod = mod,
              final_linear_range = as.numeric(names(dataFinal)),
              dataValidated = dataValidated,
              weightingMethod = weightingMethod,
              result_table_conc_levels = tables$result_table_conc_levels,
              result_table_obs = tables$result_table_obs
  ))

}


