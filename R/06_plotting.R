
#' Plot the calibration curve
#'
#' @param RES **list** \cr Results of \code{\link{CalibraCurve}}.
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range of the calibration curve.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on the plot.
#'
#' @importFrom magrittr %>%
#'
#' @returns
#' @export
#'
#' @examples
plotCalibraCurve <- function(RES,
                             ylab = "Intensity",
                             xlab = "Concentration",
                             show_regression_info = FALSE,
                             show_linear_range = TRUE,
                             show_data_points = TRUE) {


  mod <- RES$mod
  D <- RES$result_table_obs
  D$measurement <- D$measurement + 1 # add 1 to avoid log10(0) in the plot



  ## generate data for the calibration curve
  grid <- seq(log10(min(D$concentration)), log10(max(D$concentration)), length.out = 1000)
  pred <- predict(mod, newdata = data.frame(Concentration = 10^grid))
  CC_dat <- data.frame(Concentration = 10^grid, predicted = pred)
  CC_dat <- CC_dat[CC_dat$predicted > 0, ] # remove negative values in prediction (causes problems in log10-transformation later)


  ### initialize plot
  pl <- ggplot2::ggplot(D, ggplot2::aes(x = concentration, y = measurement))

  ### add data points
  if (show_data_points) {
    pl <- pl +
      ggplot2::geom_point(size = 1.7, ggplot2::aes(alpha = final_linear_range))
  }

  # log10 transformation
  pl <- pl +
    ggplot2::scale_x_continuous(trans = "log10", labels = scales::comma) +
    ggplot2::scale_y_continuous(trans = "log10")


  ### add calibration curve
  pl <- pl +
    ggplot2::geom_line(
      color = "red",
      data = CC_dat,
      ggplot2::aes(x = Concentration, y = predicted),  #log10
      inherit.aes = FALSE
    )

  pl


  ### TODO: plot multiple curves in one plot
  ### TODO: plot multiple curves with facet_wrap


  ### add shading for linear range
  if (show_linear_range) {
    pl <- pl +
      ggplot2::annotate(geom = "rect",
               xmin = min(RES$final_linear_range),
               ymin = min(CC_dat$predicted)*0.01, # 0 causes warning because of log10-transformation
               xmax = max(RES$final_linear_range),
               ymax = Inf,
               fill = "black",
               alpha = 0.1)

  }

  ## theme and labels
  pl <- pl +
    ggplot2::theme_bw() +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme(legend.position = "bottom",
                   plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm") ) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1))


  ### TODO: r^2 + model formula
  # if (show_regression_info) {
  #   pl <- pl + ggplot2::geom_text(
  #     data = eq_dat,
  #     ggplot2::aes(-Inf, Inf, label = eq),
  #     vjust = 3,
  #     hjust = -0.1
  #   )
  # }

  return(pl)
}






plotResponseFactors <- function(RES,
                                ylab = "Response Factor",
                                xlab = "Concentration") {

  range_dat <- calib_list$result_table_obs
  sum_dat <- calib_list$result_table_conc_levels

  all_rf_mean <- mean(range_dat$response_factor[range_dat$final_linear_range])

  ggplot(
    range_dat,
    aes(
      x = log10(concentration),
      y = response_factor,
      color = final_linear_range,
      fill = final_linear_range,
      alpha = final_linear_range,
      group = 1
    )
  ) +
    geom_point(size = 1.7, shape = 21) +
    geom_point(
      data = sum_dat,
      aes(x = log10(concentration), y = mean_response_factor),
      color = "black",
      shape = 21,
      size = 2.5
    ) +
    geom_line(data = sum_dat, aes(x = log10(concentration), y = mean_response_factor)) +
    scale_x_continuous(
      labels = function(x)
        sub(
          "\\.?0+$",
          "",
          format(10^x, scientific = FALSE, zero.print = FALSE)
        )
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    facet_wrap(substance ~ ., scales = "free") +
    geom_hline(
      yintercept = all_rf_mean * 1.2,
      linetype = "dashed",
      color = "orange"
    ) +
    geom_hline(
      yintercept = all_rf_mean * 0.8,
      linetype = "dashed",
      color = "orange"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab(ylab) +
    xlab(xlab)
}



