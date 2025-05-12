
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


  range_dat <- RES$result_table_obs
  final_mod <- RES$mod
  # weight_m <- RES$weightingMethod

  ### extract model coefficients and r2 per substance
  ### TODO: each substance may have a separate model, how to deal with that?
  ### TODO: there may be multiple models per substance (due to different weighting methods), how to deal with that?
  range_dat2 <- range_dat %>%
    dplyr::filter(final_linear_range) %>%
    dplyr::group_by(substance) %>%
    dplyr::reframe(predd = list(final_mod$coefficients),
                   r2 = summary(final_mod)$r.squared)

  range_dat2 <- range_dat2 %>%
    tidyr::unnest(predd) %>%
    dplyr::mutate(partype = rep(c("intercept", "coeff"), dplyr::n() / 2)) %>%
    tidyr::pivot_wider(names_from = partype, values_from = predd)

  ### extract minimum and maximum concentration per substance (independent from linear range)
  subst_range <- range_dat %>%
    dplyr::group_by(substance) %>%
    dplyr::summarise(mincon = min(concentration),
              maxcon = max(concentration))

  ### generate equation text (model formula and R squared)
  eq_dat <- range_dat2 %>%
    dplyr::mutate(
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
      r2_eq = paste0("R^2 = ", round(r2, 4))
    )


  ### generate predicted values for the calibration curve and add them to range_dat
  range_dat3 <- dplyr::left_join(range_dat, range_dat2, by = "substance") %>%
    dplyr::mutate(predicted = intercept + coeff * concentration)

  ### data frame with lower and upper limits of quantification for each substance
  an_dat <- range_dat3 %>% dplyr::group_by(substance) %>%
    dplyr::summarise(LLOQ = min(concentration[final_linear_range]),
              ULOQ = max(concentration[final_linear_range]))


  ### generate grid for predicting values of the calibration curve from the model formula
  ### xo are the values on the x-axis, predicted are on the y-axis
  range_dat4 <- dplyr::left_join(range_dat2, subst_range, by = "substance") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(xo = list(10^(seq(
      log10(mincon), log10(maxcon), length.out = 1000
    )))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(xo) %>%
    dplyr::mutate(predicted = intercept + coeff * xo)



  ### initialize plot
  pl <- ggplot2::ggplot(range_dat3, ggplot2::aes(x = concentration, y = measurement))

  ### add calibration curve
  pl <- pl +
    ggplot2::geom_line(
      color = "red",
      data = range_dat4,
      ggplot2::aes(x = xo, y = predicted),  #log10
      inherit.aes = FALSE
    )

  ### add data points
  pl <- pl +
    ggplot2::geom_point(size = 1.7, ggplot2::aes(alpha = final_linear_range))

  ### add shading for linear range
  if (show_linear_range) {
    pl <- pl +
      ggplot2::geom_rect(
        data = an_dat,
        ggplot2::aes(
          xmin = LLOQ,
          xmax = ULOQ,
          ymin = 0,
          ymax = Inf
        ),
        alpha = 0.1,
        inherit.aes = FALSE
      )
  }

  ## theme and labels
  pl <- pl +
    ggplot2::theme_bw() +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme(legend.position = "bottom",
                   plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm") )


  pl <- pl +
    #theme(legend.position = "none") +
    ggplot2::facet_wrap(substance ~ ., scales = "free") +
    ggplot2::scale_x_continuous(trans = "log10", labels = scales::comma) +
    ggplot2::scale_y_continuous(trans = "log10") +

    # ggplot2::scale_x_continuous(
    #   labels = function(x)
    #     sub(
    #       "\\.?0+$",
    #       "",
    #       format(10^x, scientific = FALSE, zero.print = FALSE)
    #     )
    # ) +
    # ggplot2::scale_y_continuous(
    #   labels = function(y)
    #     paste0(round(10^y, 5))
    # ) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1))




  if (show_regression_info) {
    pl <- pl + ggplot2::geom_text(
      data = eq_dat,
      ggplot2::aes(-Inf, Inf, label = eq),
      vjust = 3,
      hjust = -0.1
    )
  }

  return(pl)
}
#plotCalibCurve()























plotRespFactor <- function(calib_list = calibres,
                           ylab = "Response Factor",
                           xlab = "Concentration[mg/L]") {
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
#plotRespFactor()

