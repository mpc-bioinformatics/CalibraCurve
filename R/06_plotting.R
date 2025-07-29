#' Plot the calibration curve
#'
#' @param CC_RES **list** \cr Results of \code{\link{CalibraCurve}}.
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression
#'    information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range
#'    of the calibration curve as a rectangle in the plot.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on
#'    the plot.
#' @param plot_type **character(1)** \cr Type of plot for calibration curves:
#'    "single_plots" (default, generate a separate plot for each substance),
#'    "multiplot" (generate a graphic with subplots for each substance)
#'    or "all_in_one" (generate a single plot with all substances).
#' @param point_colour **character(1)** \cr Colour of the data points, default
#'    is "black".
#' @param curve_colour **character(1)** \cr Colour of the calibration curve,
#'    default is "red".
#' @param linear_range_colour **character(1)** \cr Colour of the linear range
#'    background, default is "black" (colour is weakened by alpha = 0.1).
#' @param multiplot_nrow **integer(1)** \cr Number of rows for the multiplot
#'    layout (default is NULL, which means that there is only one row).
#' @param multiplot_ncol **integer(1)** \cr Number of columns for the multiplot
#'    layout (default is NULL, which means that the number of columns is
#'    determined automatically based on the number of curves).
#' @param multiplot_scales **character(1)** \cr Scales for the multiplot layout,
#'    default is "free" (which means that each plot has its own scales). Other
#'    options are "fixed", "free_x", "free_y".
#'
#' @importFrom magrittr %>%
#'
#' @returns
#' A ggplot2 object containing the calibration curve plot.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#' package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
#'                      assayNumber = 1)
#' CC_RES <- calc_single_curve(D_list[[1]], calcContinuousPrelimRanges = FALSE)
#'
#' plotCalibraCurve(list(substance = CC_RES))
plotCalibraCurve <- function(CC_RES,
    ylab = "Intensity",
    xlab = "Concentration",
    show_regression_info = FALSE,
    show_linear_range = TRUE,
    show_data_points = TRUE,
    plot_type = "multiplot",
    point_colour = "black",
    curve_colour = "red",
    linear_range_colour = "black",
    multiplot_nrow = NULL,
    multiplot_ncol = NULL,
    multiplot_scales = "free") {
    checkmate::assertCharacter(ylab, len = 1)
    checkmate::assertCharacter(xlab, len = 1)
    checkmate::assertFlag(show_regression_info)
    checkmate::assertFlag(show_linear_range)
    checkmate::assertFlag(show_data_points)
    checkmate::assertCharacter(point_colour, len = 1)
    checkmate::assertCharacter(curve_colour, len = 1)
    checkmate::assertCharacter(linear_range_colour, len = 1)
    # silence notes when checking the package
    concentration <- measurement <- final_linear_range <- predicted <- NULL
    LLOQ <- ULOQ <- eq <- substance <- NULL

    D_calib <- .prepareCalibData(CC_RES)
    annotation_dat <- .prepareAnnotationData(D_calib)
    curve_dat <- .prepareCurveData(CC_RES, D_calib)

    ### initialize plot
    pl <- ggplot2::ggplot(D_calib, ggplot2::aes(x = concentration, y = measurement, alpha = final_linear_range))

    ### log10-transformation of axis
    pl <- pl + ggplot2::scale_x_continuous(trans = "log10", labels = scales::label_comma(drop0trailing = TRUE)) +
        ggplot2::scale_y_continuous(trans = "log10")


    ##################################################################################################
    ###### all in one plot
    if (plot_type == "all_in_one") {
        ### add data points
        if (show_data_points) {
            pl <- pl +
                ggplot2::geom_point(size = 1.7, ggplot2::aes(group = substance, colour = substance), show.legend = c(alpha = TRUE, colour = FALSE)) +
                ggplot2::scale_alpha_manual(values = c("Yes" = 1, "No" = 0.1), name = "Linear range", drop = FALSE)
        }

        ### add calibration curve
        pl <- pl +
            ggplot2::geom_line(
                data = curve_dat,
                ggplot2::aes(x = concentration, y = predicted, color = substance, group = substance),
                inherit.aes = FALSE,
                show.legend = c(colour = TRUE, alpha = FALSE)
            )
    }


    ##################################################################################################
    ###### multiplot (with facets)
    if (plot_type == "multiplot") {
        ### add data points
        if (show_data_points) {
            pl <- pl +
                ggplot2::geom_point(size = 1.7, color = point_colour, show.legend = c(alpha = TRUE, colour = FALSE)) +
                ggplot2::scale_alpha_manual(values = c("Yes" = 1, "No" = 0.1), name = "Linear range", drop = FALSE) +
                ggplot2::facet_wrap(substance ~ ., scales = multiplot_scales, nrow = multiplot_nrow, ncol = multiplot_ncol)
        }

        ### add calibration curve
        pl <- pl +
            ggplot2::geom_line(color = curve_colour, data = curve_dat,
                ggplot2::aes(x = concentration, y = predicted),
                inherit.aes = FALSE, show.legend = FALSE)

        if (show_linear_range) {
                pl <- pl + ggplot2::geom_rect(data = annotation_dat,
                    ggplot2::aes(xmin = LLOQ, xmax = ULOQ, ymin = 0, ymax = Inf),
                    alpha = 0.1, inherit.aes = FALSE, fill = linear_range_colour)
        }


        if (show_regression_info) {
                pl <- pl + ggplot2::geom_text(data = annotation_dat,
                    mapping = ggplot2::aes(y = Inf, x = 0, label = eq),
                    vjust = 1.4, hjust = -0.07, alpha = 0.6, inherit.aes = FALSE,
                    color = curve_colour, size = 3)
        }
    }


    ## theme, guides and labels
    pl <- pl +
        ggplot2::guides(
            alpha = ggplot2::guide_legend(title = "Linear range", reverse = TRUE, order = 2),
            colour = ggplot2::guide_legend(title = "Substance", order = 1)
        ) +
        ggplot2::theme_bw() + ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm"))

    return(list(CC_plot = pl, annotation_dat = annotation_dat))
}







.prepareCalibData <- function(CC_RES) {
  ## remove empty results from list (if calculation stopped)
  is_RES_null <- vapply(CC_RES, is.null, logical(1))
  CC_RES <- CC_RES[!is_RES_null]
  # X: result of 1 CalibraCurveRun
  f <- function(X) {
    tmp <- data.frame(X$result_table_obs,
                      intercept = unname(X$mod$coefficients[1]),
                      coeff = unname(X$mod$coefficients[2]),
                      r2 = summary(X$mod)$r.squared,
                      weight_m = X$weightingMethod)
    return(tmp)
  }
  D_calib <- do.call(rbind, lapply(CC_RES, f))
  D_calib$measurement[D_calib$measurement == 0] <- NA # set 0s to NA to avoid log10(0) for the plot
  D_calib$final_linear_range <- factor(D_calib$final_linear_range, levels = c(FALSE, TRUE), labels = c("No", "Yes")) # convert to logical for ggplot2 aesthetics
  return(D_calib)
}


.prepareAnnotationData <- function(D_calib) {
  D_calib_split <- split(D_calib, D_calib$substance)

  f <- function(X) {
    intercept <- X$intercept[1]
    coeff <- X$coeff[1]
    r2 <- X$r2[1]
    eq <- paste0("y = ", format(intercept, scientific = TRUE, digits = 2), " + ",
                 format(coeff, scientific = TRUE, digits = 2), " * x",
                 " (R2 = ", round(r2, 3), ")")


    tmp <- data.frame(substance = X$substance[1],
                      intercept = intercept,
                      coeff = coeff,
                      r2 = r2,
                      LLOQ = min(X$concentration[X$final_linear_range == "Yes"]),
                      ULOQ = max(X$concentration[X$final_linear_range == "Yes"]),
                      eq = eq)
    return(tmp)
  }

  annotation_dat <- do.call(rbind, lapply(D_calib_split, f))
  return(annotation_dat)
}


.prepareCurveData <- function(CC_RES, D_calib) {
  ## remove empty results from list (if calculation stopped)
  is_RES_null <- vapply(CC_RES, is.null, logical(1))
  CC_RES <- CC_RES[!is_RES_null]
  D_calib_split <- split(D_calib, D_calib$substance)
  substances <- names(D_calib_split)


  f <- function(X, Y, Z) { #X = CC_RES, Y = D_calib, Z = substances
    grid <- seq(log10(min(Y$concentration)), log10(max(Y$concentration)), length.out = 1000)
    pred <- stats::predict(X$mod, newdata = data.frame(Concentration = 10^grid))
    curve_dat_tmp <- data.frame(substance = Z, concentration = 10^grid, predicted = pred)
    curve_dat_tmp <- curve_dat_tmp[curve_dat_tmp$predicted > 0, ] # remove negative values in prediction (causes problems in log10-transformation later)
    #curve_dat <- rbind(curve_dat, curve_dat_tmp)
  }



  curve_dat <- do.call(rbind, mapply(f, CC_RES, D_calib_split, substances, SIMPLIFY = FALSE))
  return(curve_dat)
}




#' Plot response factors
#'
#' @param RES **list** \cr Results of \code{\link{calc_single_curve}}.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in percent (default is 120).
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param colour_threshold **character(1)** \cr Colour for horizontal threshold lines, default is "orange".
#' @param colour_within **character(1)** \cr Colour for points and lines within the threshold, default is "#00BFC4" (default ggplot colour).
#' @param colour_outside **character(1)** \cr Colour for horizontal outside of the threshold, default is "#F8766D" (default ggplot colour).
#'
#' @returns
#' A ggplot2 object containing the response factor plot.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#' package = "CalibraCurve")
#' D_list <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
#'                      assayNumber = 1)
#' CC_RES <- calc_single_curve(D_list[[1]], calcContinuousPrelimRanges = FALSE)
#'
#' plotResponseFactors(RES = CC_RES)
plotResponseFactors <- function(RES, RfThresL = 80, RfThresU = 120,
    ylab = "Response Factor", xlab = "Concentration",
    colour_threshold = "orange", colour_within = "#00BFC4",
    colour_outside = "#F8766D") {
    checkmate::assertNumeric(RfThresL, lower = 0, upper = 100, finite = TRUE)
    checkmate::assertNumeric(RfThresU, lower = 100)
    checkmate::assertCharacter(ylab, len = 1)
    checkmate::assertCharacter(xlab, len = 1)
    checkmate::assertCharacter(colour_threshold, len = 1)
    checkmate::assertCharacter(colour_within, len = 1)
    checkmate::assertCharacter(colour_outside, len = 1)

    # silence notes when checking the package
    concentration <- response_factor <- final_linear_range <- NULL
    mean_response_factor <- RF_within_thres <-  NULL

    range_dat <- RES$result_table_obs
    range_dat <- range_dat[!is.na(range_dat$response_factor), ]
    sum_dat <- RES$result_table_conc_levels
    sum_dat <- sum_dat[!is.na(sum_dat$mean_response_factor), ]
    all_rf_mean <- mean(range_dat$response_factor[range_dat$final_linear_range])

    ### initialize plot
    pl <- ggplot2::ggplot(range_dat, mapping = ggplot2::aes(x = concentration,
            y = response_factor, color = RF_within_thres, fill = RF_within_thres,
            alpha = RF_within_thres, group = 1))

    ### log10 transformation of x-axis
    pl <- pl + ggplot2::scale_x_continuous(trans = "log10", labels = scales::comma)

    ### add data points + mean response factors per concentration level
    pl <- pl + ggplot2::geom_point(size = 1.7, shape = 21) +
        ggplot2::geom_point(data = sum_dat,
          mapping = ggplot2::aes(x = concentration, y = mean_response_factor),
          color = "black", shape = 21, size = 2.5)

    ### adjust colour of lines (only within_colour for consecutive mean response
    ### factors within the threshold)
    sum_dat2 <- sum_dat
    for(i in 1:(nrow(sum_dat)-1)) {
      ### if the mean RF is within the threshold and the next one is too:
      if (sum_dat$RF_within_thres[i] == TRUE & sum_dat$RF_within_thres[i+1] == TRUE) {
        sum_dat2$RF_within_thres[i] <- TRUE
      } else {
        sum_dat2$RF_within_thres[i] <- FALSE
      }
    }

    ### add line between mean response factors
    pl <- pl + ggplot2::geom_line(data = sum_dat2,
            mapping = ggplot2::aes(x = concentration, y = mean_response_factor))

    ## scaling of alpha and colours
    pl <- pl +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
        ggplot2::scale_colour_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside)) +
        ggplot2::scale_fill_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside))

    ### add horizontal lines for response factor thresholds
    pl <- pl +
        ggplot2::geom_hline(yintercept = all_rf_mean * (RfThresU / 100),
            linetype = "dashed", color = colour_threshold) +
        ggplot2::geom_hline(yintercept = all_rf_mean * (RfThresL / 100),
            linetype = "dashed", color = colour_threshold)

    ### theme and axis limits
    pl <- pl +ggplot2::theme_bw() + ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
        ggplot2::theme(legend.position = "none",
                       plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm"))
    return(pl)
}
