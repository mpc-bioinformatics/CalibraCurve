#' CalibraCurve
#'
#' @param D_list **data.frame/list** \cr Data.frame or list of data.frames
#'      (one for one specific substance, each must contain two columns:
#'      Concentration and Measurement)
#' @param output_path **character(1)** \cr Folder to save results (table and
#'      plots). If NULL (default), results are not saved.
#' @param substance **character(1)** \cr Name of the substance (default is
#'      "substance1"). Will be added to the result files and may be used when
#'      plotting multiple calibration curves in one plot.
#' @param minReplicates **integer(1)** \cr Minimal number of replicates/data
#'      points per concentration level. Concentration levels with too few data
#'      points will be removed.
#' @param cvThres **numeric(1)** \cr Threshold for CV per concentration level in
#'      percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest
#'      continuous range is selected (default is TRUE). If FALSE, gaps with CVs
#'      larger than the threshold may be included.
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
#'      considered), default is 10.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in
#'      percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in
#'      percent (default is 120).
#' @param plot_type **character(1)** \cr Type of plot for calibration curves:
#'      "single_plots" (default, generate a separate plot for each substance),
#'      "multiplot" (generate a graphic with subplots for each substance) or
#'      "all_in_one" (generate a single plot with all substances).
#' @param RF_colour_threshold **character(1)** \cr Response factor plot: Colour
#'      for horizontal threshold lines, default is "orange".
#' @param RF_colour_within **character(1)** \cr Response factor plot: Colour
#'      for points and lines within the final linear range, default is "#00BFC4"
#'      (default ggplot colour).
#' @param RF_colour_outside **character(1)** \cr Response factor plot: Colour
#'      for horizontal outside of the final linear range, default is "#F8766D"
#'      (default ggplot colour).
#' @param device **character(1)** \cr Device for saving the plot (default is
#'      "png"). Other options include "pdf", "jpeg", "tiff", "svg" etc. For
#'      details see \code{\link[ggplot2]{ggsave}}.
#' @param CC_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param CC_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param RF_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param RF_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param base_size **numeric(1)** \cr Base size for the plot theme,
#'    default is 11.
#' @param RF_legend **logical(1)** \cr If TRUE, a legend will be added to the
#'      response factor plots.
#' @param plot_dpi **numeric(1)** \cr Plot resolution in dpi (default is 300).
#' @param verbose **logical(1)** \cr If FALSE, no messages will be printed.
#' @param ... additional parameters for \code{\link{plotCalibraCurve}}
#' @returns List with the following elements:
#' - \code{RES}: List of CalibraCurve results (one item per substance, output
#'      from \code{\link{calc_single_curve}}).
#' - \code{summary_tab}: Data.frame with summary information about the
#'      calibration curves, one row per substance.
#' - \code{plot_CC_list}: List of ggplot2 objects of the calibration curves
#'      (one per substance for "single_plots", otherwise only one item).
#' - \code{plot_RF_list}: List of ggplot2 objects of the response factor plots
#'      (one item).
#'
#' @importFrom checkmate assert_choice assert_numeric
#'
#' @export
#'
#' @examples
#' ### NOTE: as output_path is not specified here, no files will be saved.
#' ### Set output_path to a folder of your choice to save the results.
#'
#' ### single xlsx file:
#' data_path <- system.file("extdata", "MSQC1_xlsx", "GGPFSDSYR_QTRAP_y5.xlsx",
#'     package = "CalibraCurve")
#' D <- readDataTable(dataPath = data_path, concCol = 16, measCol = 12,
#'     fileType = "xlsx")
#' RES <- CalibraCurve(D_list = D)
#' RES$plot_CC_list
#' RES$plot_RF_list
#'
#' ### multiple xlsx files (in a folder) as multiplot:
#' data_folder <- system.file("extdata", "MSQC1_xlsx", package = "CalibraCurve")
#' D_list <- readMultipleTables(dataFolder = data_folder, fileType = "xlsx",
#'     concCol = 16, measCol = 12)
#' RES <- CalibraCurve(D_list = D_list, plot_type = "multiplot")
#' RES$plot_CC_list
#'
#' ### single rds file (SummarizedExperiment) as an all-in-one plot:
#' file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
#'     package = "CalibraCurve")
#' D_list <- readDataSE(dataPath = file, concColName = "amount_fmol",
#'         substColName = "Substance", assayNumber = 1)
#' RES <- CalibraCurve(D_list, plot_type = "all_in_one")
#' RES$plot_CC_list
#'
CalibraCurve <- function(D_list, output_path = NULL, substance = "substance",
    minReplicates = 3, cvThres = 20, calcContinuousPrelimRanges = FALSE,
    weightingMethod = "1/x^2", centralTendencyMeasure = "mean",
    perBiasThres = 20, considerPerBiasCV = TRUE, perBiasDistThres = 10,
    RfThresL = 80, RfThresU = 120, plot_type = "single_plots",
    RF_colour_threshold = "orange", RF_colour_within = "#00BFC4",
    RF_colour_outside = "#F8766D", device = "png", CC_plot_width = 15,
    CC_plot_height = 10, RF_plot_width = 15, RF_plot_height = 10, base_size = 11,
    RF_legend = FALSE, plot_dpi = 300, verbose = TRUE, ...) {
    checkmate::assert_choice(device, c("eps", "ps", "tex", "pdf", "jpeg",
        "tiff", "png", "bmp", "svg", "wmf"))
    checkmate::assert_numeric(CC_plot_width, lower = 0, len = 1)
    checkmate::assert_numeric(CC_plot_height, lower = 0, len = 1)
    checkmate::assert_numeric(RF_plot_width, lower = 0, len = 1)
    checkmate::assert_numeric(RF_plot_height, lower = 0, len = 1)
    checkmate::assert_numeric(plot_dpi, lower = 0, len = 1)

    if (inherits(D_list, "data.frame")) {
        D_list <- list(D_list)
        names(D_list) <- substance
    }

    RES <- list()
    for (i in seq_along(D_list)) {
        D_tmp <- D_list[[i]]
        substance <- names(D_list)[i]
        D_tmp <- D_tmp[order(D_tmp$Concentration),]

        if (verbose) {
            message("Calculating calibration curve for ", substance, " ...")
        }

        RES_tmp <- calc_single_curve(D = D_tmp, substance = substance, # try({
            minReplicates = minReplicates, cvThres = cvThres,
            calcContinuousPrelimRanges = calcContinuousPrelimRanges,
            weightingMethod = weightingMethod,
            centralTendencyMeasure = centralTendencyMeasure,
            perBiasThres = perBiasThres, considerPerBiasCV = considerPerBiasCV,
            perBiasDistThres = perBiasDistThres,
            RfThresL = RfThresL, RfThresU = RfThresU
        )#}, silent = TRUE)

        if (inherits(RES_tmp, "try-error")) {
            message("Problem while calculating calibration curve for ",
                    substance, ": ", RES_tmp)
            RES[[i]] <- NULL
            next
        } else {
            RES[[i]] <- RES_tmp
            names(RES)[i] <- substance
        }
    }

    pl_CC_list <- list()
    if (plot_type == "single_plots") {
        summary_tab <- NULL
        for (i in seq_along(RES)) {
            if (is.null(RES[[i]])) next
            RES_tmp <- list(RES[[i]])
            names(RES_tmp) <- names(RES)[i]
            ## generate a the calibration curve plot
            pl_CC <- plotCalibraCurve(CC_RES = RES_tmp, plot_type = "multiplot",
                    base_size = base_size, ...)
            pl_CC_list[[i]] <- pl_CC$CC_plot
            summary_tab <- rbind(summary_tab, pl_CC$annotation_dat)
        }
    } else { # plot_type == "allinone" or "multiplot"
        pl_CC <- plotCalibraCurve(CC_RES = RES, plot_type = plot_type,
                                  base_size = base_size, ...)
        pl_CC_list <- pl_CC$CC_plot
        summary_tab <- pl_CC$annotation_dat
    }

    ## generate response factor plots (each as a single plot)
    pl_RF_list <- list()
    for (i in seq_along(RES)) {
        if (is.null(RES[[i]])) next
        pl_RF_list[[i]] <- plotResponseFactors(RES = RES[[i]],
            RfThresL = RfThresL, RfThresU = RfThresU,
            colour_threshold = RF_colour_threshold,
            colour_within = RF_colour_within, colour_outside = RF_colour_outside,
            legend = RF_legend, base_size = base_size, substance = names(RES)[i]
        )
    }
    if (!is.null(output_path)) {
        .saveTablesAndPlots(RES, output_path, pl_CC_list, pl_RF_list,
                            summary_tab, CC_plot_width, CC_plot_height,
                            RF_plot_width, RF_plot_height,
                            plot_type, plot_dpi, device)
    }
        return(list(RES = RES, summary_tab = summary_tab,
                    plot_CC_list = pl_CC_list, plot_RF_list = pl_RF_list))
}




#' Calculate all necessary steps for a single calibration curve
#'
#' @param D **data.frame** data set, e.g. result of \code{\link{readDataTable}}
#'      or \code{\link{readDataSE}}. Has to include exactly two columns:
#'      "Concentration" and "Measurement".
#' @param substance **character(1)** \cr Name of the substance (default is
#'      "substance"). Will be added to the result files and may be used when
#'      plotting multiple calibration curves in one plot.
#' @param minReplicates **integer(1)** \cr Minimal number of replicates per
#'      concentration level. Concentration levels with too few data points will
#'      be removed.
#' @param cvThres **numeric(1)** \cr Threshold for CV per concentration level
#'      in percent (default is 20).
#' @param calcContinuousPrelimRanges **logical(1)** \cr If TRUE, the longest
#'      continuous range is selected (default is TRUE). If FALSE, gaps with CVs
#'      larger than the threshold may be included.
#' @param weightingMethod **character(1)** \cr Method for weighting (currently
#'      "1/x", "1/x^2" and "None" are supported, default is 1/x^2).
#' @param centralTendencyMeasure **character(1)** \cr Method for calculating
#'      average percent bias, "mean" (default) or "median".
#' @param perBiasThres **numeric(1)** \cr Threshold for average percent bias
#'      in percent, default is 20.
#' @param considerPerBiasCV **logical(1)** \cr If TRUE, CV is considered for the
#'       elimination of the concentration level (default). CV will only be
#'       considered if the difference in percent bias values is lower than
#'       perBiasDistThres.
#' @param perBiasDistThres **numeric(1)** \cr Threshold for the difference in
#'      average percent bias in percent (for lower differences, CV will be
#'      considered), default is 10.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in
#'      percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in
#'      percent (default is 120).
#'
#' @returns List with the following elements:
#' - \code{mod}: lm-object containing the final linear model.
#' - \code{final_linear_range}: vector of concentration levels that are part of
#'      the final linear range.
#' - \code{weightingMethod}: weighting method used for the linear model.
#' - \code{result_table_conc_levels}: Result table with one line for each
#'      concentration level (generated by \code{\link{assemble_results}}).
#' - \code{result_table_obs}: Result table with one line per observation (e.g.
#'      individual response factors for each data point) (generated by
#'      \code{\link{assemble_results}}).
#'
#' @export
#' @examples
#'
#' data_path <- system.file("extdata", "MSQC1_xlsx", "GGPFSDSYR_QTRAP_y5.xlsx",
#'         package = "CalibraCurve")
#' D <- readDataTable(dataPath = data_path, concCol = 16, measCol = 12,
#'         fileType = "xlsx")
#' calc_single_curve(D = D)
calc_single_curve <- function(D, substance = "substance", minReplicates = 3,
                            cvThres = 20, calcContinuousPrelimRanges = FALSE,
                            weightingMethod = "1/x^2",
                            centralTendencyMeasure = "mean",
                            perBiasThres = 20, considerPerBiasCV = TRUE,
                            perBiasDistThres = 10, RfThresL = 80,
                            RfThresU = 120) {
    ## clean data
    dataCleaned <- cleanData(D, minReplicates = minReplicates)

    ## calculate preliminary linear range
    PLR_res <- calculate_PLR(dataCleaned = dataCleaned,
        cvThres = cvThres,
        calcContinuousPrelimRanges = calcContinuousPrelimRanges)

    ## calculate final linear range
    FLR_res <- calculate_FLR(PLR_res$dataPrelim,
        weightingMethod = weightingMethod,
        centralTendencyMeasure = centralTendencyMeasure,
        perBiasThres = perBiasThres, considerPerBiasCV = considerPerBiasCV,
        perBiasDistThres = perBiasDistThres)

    ### calculate response factors
    RFs <- calcRF(dataCleaned, mod = FLR_res$mod)

    #### generate result tables
    tables <- assemble_results(X = D, dataCleaned = dataCleaned,
        cvThres = cvThres, PLR_res = PLR_res, resFacDataV = RFs$RFs,
        avgResFacDataV = RFs$meanRFs, FLR_res = FLR_res, mod = FLR_res$mod,
        RfThresL = RfThresL, RfThresU = RfThresU, substance = substance)
    RES <- list(mod = FLR_res$mod,
        final_linear_range = as.numeric(names(FLR_res$dataFinal)),
        weightingMethod = weightingMethod,
        result_table_conc_levels = tables$result_table_conc_levels,
        result_table_obs = tables$result_table_obs)
    return(RES)
}





#' Save result tables and plots
#'
#' @param RES **list** \cr Results of \code{\link{CalibraCurve}}. Each list
#'      element is the result of \code{\link{calc_single_curve}}
#' @param output_path **character(1)** \cr Folder to save results (table and
#'      plots). If NULL (default), results are not saved.
#' @param pl_CC_list **list** \cr List of calibration curves (ggplot objects,
#'      results of \code{\link{plotCalibraCurve}} (CC_plot)).
#' @param pl_RF_list **list** \cr List of response factor plots (ggplot objects,
#'      results of \code{\link{plotResponseFactors}}).
#' @param summary_tab **data.frame** \cr Table with summary information, result
#'      of \code{\link{plotCalibraCurve}} (annotation_dat))
#' @param CC_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param CC_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param RF_plot_width **numeric(1)** \cr Plot width in cm (default is 10).
#' @param RF_plot_height **numeric(1)** \cr Plot height in cm (default is 10).
#' @param plot_type **character(1)** \cr Type of plot for calibration curves:
#'      "single_plots" (default, generate a separate plot for each substance),
#'      "multiplot" (generate a graphic with subplots for each substance) or
#'      "all_in_one" (generate a single plot with all substances).
#' @param plot_dpi **numeric(1)** \cr Plot resolution in dpi (default is 300).
#' @param device **character(1)** \cr Device for saving the plot (default is
#'      "png"). Other options include "pdf", "jpeg", "tiff", "svg" etc. For
#'      details see \code{\link[ggplot2]{ggsave}}.
#'
#' @returns
#' Invisible NULL. Plots and the summary table are saved to the hard drive.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom openxlsx write.xlsx
.saveTablesAndPlots <- function(RES, output_path, pl_CC_list, pl_RF_list,
                                summary_tab, CC_plot_width, CC_plot_height,
                                RF_plot_width, RF_plot_height,
                                plot_type, plot_dpi, device) {

    for (i in seq_along(RES)) {
        .saveCCResult(CC_res = RES[[i]], output_path = output_path,
                    suffix = paste0("_", names(RES)[i]))
        if (plot_type == "single_plots") {
            ggplot2::ggsave(filename = paste0(output_path, "/CalibraCurve_",
                                            names(RES)[i], ".", device),
                plot = pl_CC_list[[i]], device = device, width = CC_plot_width,
                height = CC_plot_height, units = "cm", dpi = plot_dpi)
        }
        ggplot2::ggsave(
            filename = paste0(output_path, "/ResponseFactors_",
                            names(RES)[i], ".", device),
            plot = pl_RF_list[[i]], device = device, width = RF_plot_width,
            height = RF_plot_height, units = "cm", dpi = plot_dpi)
    }

    if (plot_type == "all_in_one" | plot_type == "multiplot") {
        ggplot2::ggsave(
            filename = paste0(output_path, "/CalibraCurve", ".", device),
            plot = pl_CC_list, device = device, width = CC_plot_width,
            height = CC_plot_height, units = "cm", dpi = plot_dpi)
    }
    openxlsx::write.xlsx(summary_tab,
                        file = paste0(output_path,
                                    "/summarytable_calibration_models.xlsx"))
    return(invisible(NULL))
}







