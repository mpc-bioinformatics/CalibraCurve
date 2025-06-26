
#' CalibraCurve
#'
#' @param data_path **character(1)** \cr Path to the data file (.csv, .txt or .xlsx file).
#' @param data_folder **character(1)** \cr Path to the folder containing multiple files (.csv, .txt or .xlsx file).
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
#' @param plot_type **character(1)** \cr Type of plot for calibration curves: "single_plots" (default, generate a separate plot for each substance),
#'                                        "multiplot" (generate a graphic with subplots for each substance)
#'                                        or "all_in_one" (generate a single plot with all substances).
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range of the calibration curve as a rectangle in the plot.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on the plot.
#' @param point_colour **character(1)** \cr Colour of the data points, default is "black".
#' @param curve_colour **character(1)** \cr Colour of the calibration curve, default is "red".
#' @param linear_range_colour **character(1)** \cr Colour of the linear range background, default is "black" (colour is weakened by alpha = 0.1).
#' @param multiplot_nrow **integer(1)** \cr Number of rows for the multiplot layout (default is NULL, which means that there is only one row).
#' @param multiplot_ncol **integer(1)** \cr Number of columns for the multiplot layout (default is NULL, which means that the number of columns is determined automatically based on the number of curves).
#' @param multiplot_scales **character(1)** \cr Scales for the multiplot layout, default is "fixed" (which means that each plot has the same scale). Other options are "free", "free_x", "free_y".
#' @param RF_colour_threshold **character(1)** \cr Response factor plot: Colour for horizontal threshold lines, default is "orange".
#' @param RF_colour_within **character(1)** \cr Response factor plot: Colour for points and lines within the final linear range, default is "#00BFC4" (default ggplot colour).
#' @param RF_colour_outside **character(1)** \cr Response factor plot: Colour for horizontal outside of the final linear range, default is "#F8766D" (default ggplot colour).
#' @param device **character(1)** \cr Device for saving the plot (default is "png"). Other options include "pdf", "jpeg", "tiff", "svg" etc. For details see \code{\link[ggplot2]{ggsave}}.
#' @param CC_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param CC_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param RF_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param RF_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param plot_dpi **numeric(1)** \cr Plot resolution in dpi (default is 300).
#' @returns List with the following elements:
#' - \code{RES}: List of CalibraCurve results (one item per substance, output from \code{\link{calc_single_curve}}).
#' - \code{summary_tab}: Data.frame with summary information about the calibration curves, one row per substance.
#' - \code{plot_CC_list}: List of ggplot2 objects of the calibration curves (one per substance for "single_plots", otherwise only one item).
#' - \code{plot_RF_list}: List of ggplot2 objects of the response factor plots (one item).
#' @export
#'
#' @examples
#'
#' ### NOTE: as output_path is not specified here, no files will be saved. Set output_path to a folder
#' ### of your choice to save the results.
#'
#' ### single file:
#' data_path <- system.file("extdata", "xlsx/ALB_LVNEVTEFAK_y8.xlsx", package = "CalibraCurve")
#' RES <- CalibraCurve(data_path = data_path, conc_col = 6, meas_col = 7)
#'#'
#' ### multiple files (in a folder) as multiplot:
#' data_folder <- system.file("extdata", "xlsx", package = "CalibraCurve")
#' RES <- CalibraCurve(data_folder = data_folder, conc_col = 6, meas_col = 7, plot_type = "multiplot")
#'
#' ### multiple files (in a folder) as all in one plot:
#' data_folder <- system.file("extdata", "xlsx", package = "CalibraCurve")
#' RES <- CalibraCurve(data_folder = data_folder, conc_col = 6, meas_col = 7, plot_type = "all_in_one")
CalibraCurve <- function(data_path = NULL,
                         data_folder = NULL,
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

                         min_replicates = 1,

                         cv_thres = 20,
                         calcContinuousPrelimRanges = FALSE,

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
                         plot_type = "single_plots",
                         show_regression_info = FALSE,
                         show_linear_range = TRUE,
                         show_data_points = TRUE,
                         point_colour = "black",
                         curve_colour = "red",
                         linear_range_colour = "black",
                         multiplot_nrow = NULL,
                         multiplot_ncol = NULL,
                         multiplot_scales = "free",
                         RF_colour_threshold = "orange",
                         RF_colour_within = "#00BFC4",
                         RF_colour_outside = "#F8766D",
                         device = "png",
                         CC_plot_width = 15,
                         CC_plot_height = 10,
                         RF_plot_width = 15,
                         RF_plot_height = 10,
                         plot_dpi = 300
) {

  ### check arguments (only those that are not already checked inside one of the low-level functions)
  checkmate::assert_choice(device, c("eps", "ps", "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", "wmf"))
  checkmate::assert_numeric(CC_plot_width, lower = 0, len = 1)
  checkmate::assert_numeric(CC_plot_height, lower = 0, len = 1)
  checkmate::assert_numeric(RF_plot_width, lower = 0, len = 1)
  checkmate::assert_numeric(RF_plot_height, lower = 0, len = 1)
  checkmate::assert_numeric(plot_dpi, lower = 0, len = 1)

  ### TODO: check if exactly one of data_path or data_folder is given

  if (is.null(data_path)) { # if folder (with potentially multiple files) is given
    #### TODO: hier müsste nach relevanten files (xlsx, csv, usw) gefiltert werden
    all_files <- setdiff(list.files(path = data_folder), list.dirs(path = data_folder, recursive = FALSE, full.names = FALSE))

    filetable <- data.frame(file = all_files) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(substance_name = strsplit(basename(file), "\\.")[[1]][1],
                    full_path = paste0(data_folder, "/", file),
                    file_ext = tools::file_ext(file))
  } else { # if data_path is given (single file)
    all_files <- data_path
    data_folder <- dirname(data_path) # extract folder from file path
    filetable <- data.frame(file = all_files) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(substance_name = strsplit(basename(file), "\\.")[[1]][1],
                    full_path = file,
                    file_ext = tools::file_ext(file))
  }

  ### check if all files are valid (correct file type)

  all_files <- all_files[filetable$file_ext == filetype,]



  RES <- list()
  ## calculate calibration curve for each file
  for (i in 1:nrow(filetable)) {
    data_path <- filetable[i, ]$full_path
    substance <- filetable[i, ]$substance_name

    print(paste0("Calculating calibration curve for ", substance, " ..."))

    RES_tmp <- calc_single_curve(data_path = data_path,
                                 output_path = output_path,
                                 conc_col = conc_col,
                                 meas_col = meas_col,
                                 substance = substance,
                                 filetype = filetype,
                                 sep = sep,
                                 dec = dec,
                                 header = header,
                                 na.strings = na.strings,
                                 sheet = sheet,
                                 min_replicates = min_replicates,
                                 cv_thres = cv_thres,
                                 calcContinuousPrelimRanges = calcContinuousPrelimRanges,
                                 weightingMethod = weightingMethod,
                                 centralTendencyMeasure = centralTendencyMeasure,
                                 perBiasThres = perBiasThres,
                                 considerPerBiasCV = considerPerBiasCV,
                                 perBiasDistThres = perBiasDistThres,
                                 RfThresL = RfThresL,
                                 RfThresU = RfThresU)
    RES[[i]] <- RES_tmp
    names(RES)[i] <- substance

    if (!is.null(output_path)) {
      CalibraCurve::saveCCResult(CC_res = RES[[i]],
                                 output_path = output_path,
                                 suffix = paste0("_", substance))
    }
  }




  pl_CC_list <- list()
  if (plot_type == "single_plots") {

    annotation_dat <- NULL

    for (i in 1:length(RES)) {

      RES_tmp <- list(RES[[i]])
      names(RES_tmp) <- names(RES)[i]

      ## generate and ave the calibration curve plot
      pl_CC <- plotCalibraCurve(
        RES = RES_tmp,
        ylab = ylab,
        xlab = xlab,
        plot_type = "multiplot",
        show_regression_info = show_regression_info,
        show_linear_range = show_linear_range,
        show_data_points = show_data_points,
        point_colour = point_colour,
        curve_colour = curve_colour,
        linear_range_colour = linear_range_colour
      )
      pl_CC_list[[i]] <- pl_CC$CC_plot

      if (!is.null(output_path)) {
        ggplot2::ggsave(
          filename = paste0(output_path, "/CalibraCurve_", names(RES)[i], ".", device),
          plot = pl_CC$CC_plot,
          device = device,
          width = CC_plot_width,
          height = CC_plot_height,
          units = "cm",
          dpi = plot_dpi
        )
      }
      annotation_dat <- rbind(annotation_dat, pl_CC$annotation_dat)
    }

  } else { # plot_type == "allinone" or "multiplot"
    pl_CC <- plotCalibraCurve(
      RES = RES,
      ylab = ylab,
      xlab = xlab,
      plot_type = plot_type,
      show_regression_info = show_regression_info,
      show_linear_range = show_linear_range,
      show_data_points = show_data_points,
      point_colour = point_colour,
      curve_colour = curve_colour,
      linear_range_colour = linear_range_colour,
      multiplot_nrow = multiplot_nrow,
      multiplot_ncol = multiplot_ncol,
      multiplot_scales = multiplot_scales
    )
    pl_CC_list <- pl_CC$CC_plot
    if (!is.null(output_path)) {
      suppressWarnings({
        ggplot2::ggsave(
          filename = paste0(output_path, "/CalibraCurve", ".", device),
          plot = pl_CC$CC_plot,
          device = device,
          width = CC_plot_width,
          height = CC_plot_height,
          units = "cm",
          dpi = plot_dpi
        )
      })
    }
    annotation_dat <- pl_CC$annotation_dat
  }



  ## generate and save response factor plots (each as a single plot)
  pl_RF_list <- list()
  for (i in 1:length(RES)) {
    pl_RF <- plotResponseFactors(
      RES = RES[[i]],
      RfThresL = RfThresL,
      RfThresU = RfThresU,
      colour_threshold = RF_colour_threshold,
      colour_within = RF_colour_within,
      colour_outside = RF_colour_outside
    )
    pl_RF_list[[i]] <- pl_RF
    if (!is.null(output_path)) {
      ggplot2::ggsave(
        filename = paste0(output_path, "/ResponseFactors_", names(RES)[i], ".", device),
        plot = pl_RF,
        device = device,
        width = RF_plot_width,
        height = RF_plot_height,
        units = "cm",
        dpi = plot_dpi
      )
    }
  }

  summary_tab <- annotation_dat[, -ncol(annotation_dat)]
  if (!is.null(output_path)) {
    openxlsx::write.xlsx(summary_tab, file = paste0(output_path, "/summarytable_calibration_models.xlsx"))
  }

  return(list(RES = RES, summary_tab = summary_tab, plot_CC_list = pl_CC_list, plot_RF_list = pl_RF_list))
}




#' Calculate all necessary steps for a single calibration curve
#'
#' @param data_path **character(1)** \cr Path to the data file (.csv, .txt or .xlsx file).
#' @param output_path **character(1)** \cr Folder to save results (table and plots). If NULL (default), results are not saved.
#' @param conc_col **integer(1)** \cr Column number of the concentration values.
#' @param meas_col **integer(1)** \cr Column number of the concentration values.
#' @param substance **character(1)** \cr Name of the substance (default is "substance1"). Will be added to the result files and may be used when plotting multiple calibration curves in one plot.
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
#'
#' @returns List with the following elements:
#' - \code{mod}: lm-object containing the final linear model.
#' - \code{final_linear_range}: vector of concentration levels that are part of the final linear range.
#' - \code{weightingMethod}: weighting method that was used for the linear model.
#' - \code{result_table_conc_levels}: Result table with one line for each concentration level (generated by \code{\link{assemble_results}}).
#' - \code{result_table_obs}: Result table with one line per observation (e.g. individual response factors for each data point) (generated by \code{\link{assemble_results}}).
#' @export
#' @examples
#'
#' data_path <- system.file("extdata", "xlsx/ALB_LVNEVTEFAK_y8.xlsx", package = "CalibraCurve")
#' calc_single_curve(data_path = data_path, conc_col = 6, meas_col = 7)
calc_single_curve <- function(data_path,
                              output_path = NULL,
                              conc_col,
                              meas_col,
                              substance = "substance1",

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
                              perBiasThres = 30,
                              considerPerBiasCV = TRUE,
                              perBiasDistThres = 10,

                              RfThresL = 80,
                              RfThresU = 120) {
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
              weightingMethod = weightingMethod,
              result_table_conc_levels = tables$result_table_conc_levels,
              result_table_obs = tables$result_table_obs)

  return(RES)
}




