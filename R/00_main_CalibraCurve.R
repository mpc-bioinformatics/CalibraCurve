
#' CalibraCurve
#'
#' @param data_path **character(1)** \cr Path to the data file (.csv, .txt or .xlsx file).
#' @param output_path **character(1)** \cr Folder to save results (table and plots). If NULL (default), results are not saved.
#' @param conc_col **integer(1)** \cr Column number of the concentration values.
#' @param meas_col **integer(1)** \cr Column number of the concentration values.
#' @param substance **character(1)** \cr Name of the substance (default is "substance1"). Will be added to the result files and may be used when plotting multiple calibration curves in one plot.
#' @param suffix **character(1)** \cr Suffix for the output files, ideally starting with "_" (default is substance with a "_" in front).
#' @param filetype **character(1)** \cr Type of input file: "csv" or "txt" or "xlsx".
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\\t" for tab. Default is ",".
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot. Default is ".".
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param na.strings **character** \cr Character vector of strings which are to be interpreted as NA.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
#' @param min_replicates **integer(1)** \cr Minimal number of replicates/data points per concentration level.
#'                                          Concentration levels with too few data points will be removed.
#' @param cv_thres **numeric(1)** \cr Threshold for CV per concentration level in percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest continuous range is selected (default is TRUE).
#'                                                      If FALSE, gaps with CVs larger than the threshold may be included.
#' @param weightingMethod **character(1)** \cr Method for weighting (currently "1/x", "1/x^2" and "None" are supported, default is 1/x^2).
#' @param centralTendencyMeasure **character(1)** \cr Method for calculating average percent bias, "mean" (default) or "median".
#' @param perBiasThres **numeric(1)** \cr Threshold for average percent bias in percent, default is 20.
#' @param considerPerBiasCV **logical(1)** \cr If TRUE, CV is considered for the elimination of the concentration level (default). CV will only be considered if the difference in
#'                                             percent bias values is lower than perBiasDistThres.
#' @param perBiasDistThres **numeric(1)** \cr Threshold for the difference in average percent bias in percent (for lower differences, CV will be considered), default is 10.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in percent (default is 120).
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range of the calibration curve as a rectangle in the plot.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on the plot.
#' @param point_colour **character(1)** \cr Colour of the data points, default is "black".
#' @param curve_colour **character(1)** \cr Colour of the calibration curve, default is "red".
#' @param linear_range_colour **character(1)** \cr Colour of the linear range background, default is "black" (colour is weakened by alpha = 0.1).
#' @param RF_colour_threshold **character(1)** \cr Response factor plot: Colour for horizontal threshold lines, default is "orange".
#' @param RF_colour_within **character(1)** \cr Response factor plot: Colour for points and lines within the final linear range, default is "#00BFC4" (default ggplot colour).
#' @param RF_colour_outside **character(1)** \cr Response factor plot: Colour for horizontal outside of the final linear range, default is "#F8766D" (default ggplot colour).
#' @param device **character(1)** \cr Device for saving the plot (default is "png"). Other options include "pdf", "jpeg", "tiff", "svg" etc. For details see \code{\link[ggplot2]{ggsave}}.
#' @param plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param plot_dpi **numeric(1)** \cr Plot resolution in dpi (default is 300).
#' @returns List with the following elements:
#' - \code{mod}: lm-object containing the final linear model.
#' - \code{final_linear_range}: vector of concentration levels that are part of the final linear range.
#' - \code{dataCleaned}: list of data.frames, result of \code{\link{cleanData}}.
#' - \code{weightingMethod}: weighting method that was used for the linear model.
#' - \code{result_table_conc_levels}: Result table with one line for each concentration level (generated by \code{\link{assemble_results}}).
#' - \code{result_table_obs}: Result table with one line per observation (e.g. individual response factors for each data point) (generated by \code{\link{assemble_results}}).
#' - \code{plot_CC}: ggplot object of the calibration curve.
#' - \code{plot_RF}: ggplot object of the response factor plot.
#' @export
#'
#' @examples
#'
#' data_path <- system.file("extdata", "ALB_LVNEVTEFAK_y8.xlsx", package = "CalibraCurve")
#' CalibraCurve(data_path = data_path, conc_col = 6, meas_col = 7)
#'
CalibraCurve <- function(data_path,
                         output_path = NULL,
                         conc_col,
                         meas_col,
                         substance = "substance1",
                         suffix = paste0("_", substance),

                         filetype = "xlsx",
                         sep = ",",
                         dec = ".",
                         header = TRUE,
                         na.strings = c("NA", "NaN", "Filtered", "#NV"),
                         sheet = 1,

                         min_replicates = 3,

                         cv_thres = 20,
                         calcContinuousPrelimRanges = TRUE,

                         weightingMethod = "1/x^2",
                         centralTendencyMeasure = "mean",
                         perBiasThres = 20,
                         considerPerBiasCV = TRUE,
                         perBiasDistThres = 10,

                         RfThresL = 80,
                         RfThresU = 120,

                         ### parameters for plotting
                         ylab = "Intensity",
                         xlab = "Concentration",
                         show_regression_info = FALSE,
                         show_linear_range = TRUE,
                         show_data_points = TRUE,
                         point_colour = "black",
                         curve_colour = "red",
                         linear_range_colour = "black",
                         RF_colour_threshold = "orange",
                         RF_colour_within = "#00BFC4",
                         RF_colour_outside = "#F8766D",
                         device = "png",
                         plot_width = 12,
                         plot_height = 10,
                         plot_dpi = 300
                         ) {

  ### check arguments (only those that are not already checked inside one of the low-level functions)
  checkmate::assert_choice(device, c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf"))
  checkmate::assert_numeric(plot_width, lower = 0, len = 1)
  checkmate::assert_numeric(plot_height, lower = 0, len = 1)
  checkmate::assert_numeric(plot_dpi, lower = 0, len = 1)

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
  FLR_res <- CalibraCurve::calculate_FLR(PLR_res$dataPrelim,
                       weightingMethod = weightingMethod,
                       centralTendencyMeasure = centralTendencyMeasure,
                       perBiasThres = perBiasThres,
                       considerPerBiasCV = considerPerBiasCV,
                       perBiasDistThres = perBiasDistThres)

  dataFinal <- FLR_res$dataFinal
  mod <- FLR_res$mod


  ### calculate response factors
  resFacDataV <- CalibraCurve::calcRFLevels(dataCleaned, mod = mod)

  # Calculation of mean response factor values
  avgResFacDataV <- CalibraCurve::calcRFMeans(resFacDataV)


  #### generate result tables
  tables <- CalibraCurve::assemble_results(X = X,
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

  RES <- list(mod = mod,
              final_linear_range = as.numeric(names(dataFinal)),
              dataCleaned = dataCleaned,
              weightingMethod = weightingMethod,
              result_table_conc_levels = tables$result_table_conc_levels,
              result_table_obs = tables$result_table_obs)

  if (!is.null(output_path)) {
  CalibraCurve::saveCCResult(CC_res = RES,
                             output_path = output_path,
                             suffix = suffix)
  }


  ## calibration curve plot
  pl_CC <- CalibraCurve::plotCalibraCurve(RES,
                                          ylab = ylab,
                                          xlab = xlab,
                                          show_regression_info = show_regression_info,
                                          show_linear_range = show_linear_range,
                                          show_data_points = show_data_points,
                                          point_colour = point_colour,
                                          curve_colour = curve_colour,
                                          linear_range_colour = linear_range_colour
                                          )
  if (!is.null(output_path)) {
    ## save the plot
    ggplot2::ggsave(filename = paste0(output_path, "/CalibraCurve", suffix, ".", device),
                    plot = pl_CC,
                    device = device,
                    width = plot_width,
                    height = plot_height,
                    units = "cm",
                    dpi = plot_dpi)
  }
  RES$plot_CC <- pl_CC


  ## response factor plot
  pl_RF <- CalibraCurve::plotResponseFactors(RES,
                                             RfThresL = RfThresL,
                                             RfThresU = RfThresU,
                                             colour_threshold = RF_colour_threshold,
                                             colour_within = RF_colour_within,
                                             colour_outside = RF_colour_outside)

  if (!is.null(output_path)) {
    ## save the plot
    ggplot2::ggsave(filename = paste0(output_path, "/ResponseFactors", suffix, ".", device),
                    plot = pl_RF,
                    device = device,
                    width = plot_width,
                    height = plot_height,
                    units = "cm",
                    dpi = plot_dpi)
  }

  RES$plot_RF <- pl_RF


  return(RES)

}


