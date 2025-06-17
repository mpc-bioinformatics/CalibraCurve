
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
                             conc_unit = "[mg/L]",
                             point_colour = "black",
                             in_facet_wrap = TRUE,
                             curve_colour = "red",
                             linear_range_colour = "black") {



  checkmate::assertCharacter(ylab, len = 1)
  checkmate::assertCharacter(xlab, len = 1)
  checkmate::assertFlag(show_regression_info)
  checkmate::assertFlag(show_linear_range)
  checkmate::assertFlag(show_data_points)
  checkmate::assertCharacter(point_colour, len = 1)
  checkmate::assertCharacter(curve_colour, len = 1)
  checkmate::assertCharacter(linear_range_colour, len = 1)


  D_calib <- NULL
  for (i in 1:length(RES)) {
    print(i)

    mod2 <<- RES[[i]]$mod

    tmp <- data.frame(RES[[i]]$result_table_obs,
               intercept = RES[[i]]$mod$coefficients[1],
               coeff = RES[[i]]$mod$coefficients[2],
               r2 = summary(RES[[i]]$mod)$r.squared,
               weight_m = RES[[i]]$weightingMethod)
    D_calib <- rbind(D_calib, tmp)

  }

  D_calib$measurement[D_calib$measurement == 0] <- NA # set 0s to NA to avoid log10(0) for the plot


  annotation_dat <- D_calib %>%
    select(substance, concentration, intercept, coeff, r2,final_linear_range) %>%
    group_by(substance, intercept, coeff,r2) %>%
    summarise(LLOQ = min(concentration[final_linear_range]),
              ULOQ = max(concentration[final_linear_range])) %>%
    mutate(lin_range = paste0("linear range: (", LLOQ, ", ", ULOQ, ") ", conc_unit),
           eq = paste0(
             "y = ",
             format(intercept, scientific = TRUE, digits = 2),
             " + ",
             format(coeff, scientific = TRUE, digits = 2),
             " * x",
             " (R2 = ",
             round(r2, 3),
             ")"
           ),
           r2_eq = paste0("R^2 = ", round(r2, 4)))


  curve_dat <- NULL
  for (i in 1:length(RES)) {
    substance <- names(RES)[i]
    D_calib_tmp <- D_calib[D_calib$substance == substance, ]
    ## generate data for the calibration curve
    grid <- seq(log10(min(D_calib_tmp$concentration)), log10(max(D_calib_tmp$concentration)), length.out = 1000)
    pred <- stats::predict(RES[[i]]$mod, newdata = data.frame(Concentration = 10^grid))
    curve_dat_tmp <- data.frame(substance = substance, concentration = 10^grid, predicted = pred)
    curve_dat_tmp <- curve_dat_tmp[curve_dat_tmp$predicted > 0, ] # remove negative values in prediction (causes problems in log10-transformation later)
    curve_dat <- rbind(curve_dat, curve_dat_tmp)
  }

  curve_dat2 <<- curve_dat

  # curve_dat <- D_calib %>%
  #   select(substance, concentration, measurement, intercept, coeff, final_linear_range) %>%
  #   group_by(substance, intercept, coeff) %>%
  #   reframe(concentration = list(10^seq(log10(min(concentration)), log10(max(concentration)), length.out = 1000))) %>%
  #   unnest(concentration) %>%
  #   mutate(predicted = intercept + coeff * concentration) %>%
  #   filter(predicted >= 0)


  ### initialize plot
  pl <- ggplot2::ggplot(D_calib, ggplot2::aes(x = concentration, y = measurement)) +
    ggplot2::theme_bw() +
    # log10 transformation only for x-values to keep geom_rect()-functionalities
    ggplot2::scale_x_continuous(trans = "log10", labels = scales::label_comma(drop0trailing = TRUE)) +
    ggplot2::scale_y_continuous(trans = "log10")

  #+
    #ggplot2::scale_y_continuous(labels = function(y) paste0(round(10^y, 5)))

  if(!in_facet_wrap) {

    ### add data points
    if (show_data_points) {
      pl <- pl +
        ggplot2::geom_point(size = 1.7, ggplot2::aes(alpha = final_linear_range, color = substance, group = substance)) +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1))
    }

    ### add calibration curve
    pl <- pl +
      ggplot2::geom_line(
        data = curve_dat,
        ggplot2::aes(x = concentration, y = predicted, color = substance, group = substance),
        inherit.aes = FALSE
      )

    if (show_linear_range) {

      range_y <- layer_scales(pl)$y$range$range
      print(range_y)

      pl <- pl +
        ggplot2::annotate(geom = "rect",
                          xmin = annotation_dat$LLOQ,
                          ymin = min(curve_dat$predicted)*0.00001, # 0 causes warning because of log10-transformation
                          xmax = annotation_dat$ULOQ,
                          ymax = Inf,
                          fill = linear_range_colour,
                          alpha = 0.1) +
        coord_cartesian(ylim = 10^range_y)
    }

  }

  if(in_facet_wrap) {

    ### add data points
    if (show_data_points) {
      pl <- pl +
        ggplot2::geom_point(size = 1.7, ggplot2::aes(alpha = final_linear_range), color = point_colour) +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1)) +
        ggplot2::facet_wrap(substance ~ ., scales = "free")
    }

    ### add calibration curve
    pl <- pl +
      ggplot2::geom_line(
        color = curve_colour,
        data = curve_dat,
        ggplot2::aes(x = concentration, y = predicted),
        inherit.aes = FALSE
      )

  }





  if (show_regression_info) {
    pl <- pl + ggplot2::geom_text(
      data = annotation_dat,
      ggplot2::aes(y = -Inf, x = Inf, label = eq),
      vjust = -1.2,
      hjust = 1.1,
      alpha = 0.6,
      inherit.aes = FALSE,
      color = curve_colour,
      size = 3
    )
  }



  ## theme and labels
  pl <- pl +
    ggplot2::guides(alpha = "none") +
    ggplot2::theme_bw() +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme(title = ggplot2::element_text("Substance"),
                   plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm"))

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

  response_dat <- RES %>%
    mutate(result_dat = list(res$result_table_obs)) %>%
    unnest(result_dat) %>%
    select(substance_name, concentration, response_factor, final_linear_range) %>%
    drop_na(concentration, response_factor)

  mean_rf <- response_dat %>%
    group_by(substance_name, concentration, final_linear_range) %>%
    reframe(mean_rf = mean(response_factor))

  all_mean_rf <- mean_rf %>%
    filter(final_linear_range) %>%
    group_by(substance_name)%>%
    reframe(all_mean_rf = mean(mean_rf))


  ### initialize plot
  pl <- ggplot2::ggplot(
    response_dat,
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
  pl <- pl + ggplot2::scale_x_continuous(trans = "log10", labels = scales::label_comma(drop0trailing=TRUE))

  ### add data points + mean response factors per concentration level
  pl <- pl +
    ggplot2::geom_point(size = 1.7, shape = 21) + facet_wrap(substance_name ~ ., scales = "free") +
    ggplot2::geom_point(
      data = mean_rf,
      ggplot2::aes(x = concentration, y = mean_rf),
      color = "black",
      shape = 21,
      size = 2.5
    )

  ### add line between mean response factors
  pl <- pl +
    ggplot2::geom_line(data = mean_rf, ggplot2::aes(x = concentration, y = mean_rf))

  ## scaling of alpha and colours
  pl <- pl +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    ggplot2::scale_colour_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside)) +
    ggplot2::scale_fill_manual(values = c("TRUE" = colour_within, "FALSE" = colour_outside))

  ### add horizontal lines for response factor thresholds
  pl <- pl +
    ggplot2::geom_hline(data = all_mean_rf, aes(
      yintercept = all_mean_rf * (RfThresU/100)),
      linetype = "dashed",
      color = colour_threshold
    ) +
    ggplot2::geom_hline(data = all_mean_rf, aes(
      yintercept = all_mean_rf * (RfThresL/100)),
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



