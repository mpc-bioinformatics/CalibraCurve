

#' CalibraCurve
#'
#' @param data_path **character(1)** \cr Path to the data file (.csv, .txt or .xlsx file).
#' @param output_path **character(1)** \cr Folder to save results (table and plots).
#' @param conc_col **integer(1)** \cr Column number of the concentration values.
#' @param meas_col **integer(1)** \cr Column number of the concentration values.
#' @param filetype **character(1)** \cr Type of input file: "csv" or "txt" or "xlsx".
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\t" for tab.
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot.
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param na.strings **character** \cr Character vector of strings which are to be interpreted as NA.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
#' @param min_replicates **integer(1)** \cr Minimal number of replicates/data points per concentration level.
#'                                          Concentration levels with too few data points will be removed.
#' @param cv_thres
#' @param calcContinuousPrelimRanges
#' @param weightingMethod
#' @param centralTendencyMeasure
#' @param perBiasThres
#' @param considerPerBiasCV
#' @param perBiasDistThres



#' @param RfThresL
#' @param RfThresU
#' @param substance
#'
#' @returns
#' @export
#'
#' @examples
CalibraCurve <- function(data_path,
                         output_path = NULL,
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
  X <- CalibraCurve::readData(data_path = data_path,
                              conc_col = conc_col,
                              meas_col = meas_col,
                              filetype = filetype,
                              sep = sep,
                              dec = dec,
                              header = header,
                              na.strings = na.strings,
                              sheet = sheet)

  ## clean data
  dataCleaned <- CalibraCurve::cleanData(X,
                                           min_replicates = min_replicates)

  ## calculate preliminary linear range
  PLR_res <- CalibraCurve::calculate_PLR(dataCleaned = dataCleaned,
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
                             dataCleaned = dataCleaned,
                             cv_thres = cv_thres,
                             PLR_res = PLR_res,
                             resFacDataV = resFacDataV,
                             avgResFacDataV = avgResFacDataV,
                             FLR_res = FLR_res,
                             mod = mod,
                             RfThresL = RfThresL,
                             RfThresU = RfThresU,
                             substance = substance
  )



  return(list(mod = mod,
              final_linear_range = as.numeric(names(dataFinal)),
              dataCleaned = dataCleaned,
              weightingMethod = weightingMethod,
              result_table_conc_levels = tables$result_table_conc_levels,
              result_table_obs = tables$result_table_obs
  ))

}


