
#' Plot the calibration curve
#'
#' @param RES **list** \cr Results of \code{\link{CalibraCurve}}.
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range of the calibration curve as a rectangle in the plot.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on the plot.
#' @param point_colour **character(1)** \cr Colour of the data points, default is "black".
#' @param curve_colour **character(1)** \cr Colour of the calibration curve, default is "red".
#' @param linear_range_colour **character(1)** \cr Colour of the linear range background, default is "black" (colour is weakened by alpha = 0.1).
#'
#' @importFrom magrittr %>%
#'
#' @returns
#' A ggplot2 object containing the calibration curve plot.
#' @export
#'
#' @examples
#'
#' data(RES_ALB)
#' plotCalibraCurve(RES_ALB)
plotCalibraCurve <- function(RES,
                             ylab = "Intensity",
                             xlab = "Concentration",
                             show_regression_info = FALSE,
                             show_linear_range = TRUE,
                             show_data_points = TRUE,
                             point_colour = "black",
                             curve_colour = "red",
                             linear_range_colour = "black") {

  ### TODO: plot multiple curves in one plot
  ### TODO: plot multiple curves with facet_wrap

  checkmate::assertCharacter(ylab, len = 1)
  checkmate::assertCharacter(xlab, len = 1)
  checkmate::assertFlag(show_regression_info)
  checkmate::assertFlag(show_linear_range)
  checkmate::assertFlag(show_data_points)
  checkmate::assertCharacter(point_colour, len = 1)
  checkmate::assertCharacter(curve_colour, len = 1)
  checkmate::assertCharacter(linear_range_colour, len = 1)


  concentration <- measurement <- final_linear_range <- predicted <- NULL # silence notes when checking the package

  mod <- RES$mod
  D <- RES$result_table_obs
  D$measurement[D$measurement == 0] <- NA # set 0s to NA to avoid log10(0) for the plot

  ## generate data for the calibration curve
  grid <- seq(log10(min(D$concentration)), log10(max(D$concentration)), length.out = 1000)
  pred <- stats::predict(mod, newdata = data.frame(Concentration = 10^grid))
  CC_dat <- data.frame(concentration = 10^grid, predicted = pred)
  CC_dat <- CC_dat[CC_dat$predicted > 0, ] # remove negative values in prediction (causes problems in log10-transformation later)


  ### initialize plot
  pl <- ggplot2::ggplot(D, ggplot2::aes(x = concentration, y = measurement))

  ### add data points
  if (show_data_points) {
    pl <- pl +
      ggplot2::geom_point(size = 1.7, ggplot2::aes(alpha = final_linear_range), colour = point_colour) +
      ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1))
  }

  # log10 transformation
  pl <- pl +
    ggplot2::scale_x_continuous(trans = "log10", labels = scales::comma) +
    ggplot2::scale_y_continuous(trans = "log10")


  ### add calibration curve
  pl <- pl +
    ggplot2::geom_line(
      color = curve_colour,
      data = CC_dat,
      ggplot2::aes(x = concentration, y = predicted),
      inherit.aes = FALSE
    )

  pl

  ### add shading for linear range
  if (show_linear_range) {
    pl <- pl +
      ggplot2::annotate(geom = "rect",
               xmin = min(RES$final_linear_range),
               ymin = min(CC_dat$predicted)*0.01, # 0 causes warning because of log10-transformation
               xmax = max(RES$final_linear_range),
               ymax = Inf,
               fill = linear_range_colour,
               alpha = 0.1)
  }

  ## theme and labels
  pl <- pl +
    ggplot2::theme_bw() +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme(legend.position = "bottom",
                   plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm"))



  if (show_regression_info) {
    intercept <- stats::coef(mod)[1]
    slope <- stats::coef(mod)[2]
    r2 <- summary(mod)$r.squared
    eq = paste0(
      "y = ",
      format(intercept, scientific = TRUE, digits = 2),
      " + ",
      format(slope, scientific = TRUE, digits = 2),
      " * x",
      " (R2 = ",
      round(r2, 3),
      ")"
    )

    pl <- pl +
      ggplot2::annotate("text", x = min(D$concentration), y = max(D$measurement), label = eq, hjust = 0)
  }

  return(pl)
}






#' Plot response factors
#'
#' @param RES **list** \cr Results of \code{\link{CalibraCurve}}.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in percent (default is 120).
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param colour_threshold **character(1)** \cr Colour for horizontal threshold lines, default is "orange".
#' @param colour_within **character(1)** \cr Colour for points and lines within the final linear range, default is "#00BFC4" (default ggplot colour).
#' @param colour_outside **character(1)** \cr Colour for horizontal outside of the final linear range, default is "#F8766D" (default ggplot colour).
#'
#' @returns
#' A ggplot2 object containing the response factor plot.
#' @export
#'
#' @examples
#'
#' data(RES_ALB)
#' plotResponseFactors(RES_ALB)
plotResponseFactors <- function(RES,
                                RfThresL = 80,
                                RfThresU = 120,
                                ylab = "Response Factor",
                                xlab = "Concentration",
                                colour_threshold = "orange",
                                colour_within = "#00BFC4",
                                colour_outside = "#F8766D") {


  checkmate::assertNumeric(RfThresL, lower = 0, upper = 100, finite = TRUE)
  checkmate::assertNumeric(RfThresU, lower = 100)
  checkmate::assertCharacter(ylab, len = 1)
  checkmate::assertCharacter(xlab, len = 1)
  checkmate::assertCharacter(colour_threshold, len = 1)
  checkmate::assertCharacter(colour_within, len = 1)
  checkmate::assertCharacter(colour_outside, len = 1)


  concentration <- response_factor <- final_linear_range <- mean_response_factor <- NULL # silence notes when checking the package

  range_dat <- RES$result_table_obs
  range_dat <- range_dat[!is.na(range_dat$response_factor), ]
  sum_dat <- RES$result_table_conc_levels
  sum_dat <- sum_dat[!is.na(sum_dat$mean_response_factor), ]
  all_rf_mean <- mean(range_dat$response_factor[range_dat$final_linear_range])


  ### initialize plot
  pl <- ggplot2::ggplot(
    range_dat,
    ggplot2::aes(
      x = concentration,
      y = response_factor,
      color = final_linear_range,
      fill = final_linear_range,
      alpha = final_linear_range,
      group = 1
    )
  )

  ### log10 transformation of x-axis
  pl <- pl + ggplot2::scale_x_continuous(trans = "log10", labels = scales::comma)

  ### add data points + mean response factors per concentration level
  pl <- pl +
    ggplot2::geom_point(size = 1.7, shape = 21) +
    ggplot2::geom_point(
      data = sum_dat,
      ggplot2::aes(x = concentration, y = mean_response_factor),
      color = "black",
      shape = 21,
      size = 2.5
    )


  sum_dat2 <- sum_dat
  if (any(sum_dat2$final_linear_range)) {
    ind <- which(sum_dat2$final_linear_range)
    ind_last <- ind[length(ind)]
    if (ind_last != length(sum_dat2$final_linear_range))
      sum_dat2$final_linear_range[ind_last] <- FALSE
  }

  ### add line between mean response factors
  pl <- pl +
    ggplot2::geom_line(data = sum_dat2, ggplot2::aes(x = concentration, y = mean_response_factor))

  ## scaling of alpha and colours
  pl <- pl +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    ggplot2::scale_colour_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside))
  #ggplot2::facet_wrap(substance ~ ., scales = "free") +

  ### add horizontal lines for response factor thresholds
  pl <- pl +
    ggplot2::geom_hline(
      yintercept = all_rf_mean * (RfThresU/100),
      linetype = "dashed",
      color = colour_threshold
    ) +
    ggplot2::geom_hline(
      yintercept = all_rf_mean * (RfThresL/100),
      linetype = "dashed",
      color = colour_threshold
    )

  ### theme and axis limits
  pl <- pl +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm")) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab)

  return(pl)
}



