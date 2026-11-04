#' Plot the calibration curve
#' @details
#' Can be used also on pre-calculated curves, if you import the .rds file
#' using [readRDS()]. Does not save the plot, this can be done using
#' [ggplot2::ggsave].
#'
#'
#'
#' @param CC_RES **list** \cr Results of \code{\link{CalibraCurve}}. Each list
#'      element is the result of \code{\link{calc_single_curve}}
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param log10_x_axis **logical(1)** \cr If TRUE (default), use a log10
#'    transformation for the x-axis. If FALSE, keep original concentration
#'    values on a linear x-axis.
#' @param log_tick_label_style **character(1)** \cr Label style for log-scale
#'    axes. One of \code{"auto"} (default), \code{"math"}, \code{"plain"}.
#'    \code{"math"} uses labels like \eqn{10^{-5}}. \code{"plain"} uses
#'    numeric labels (for example, 0.00001). \code{"auto"} uses \code{"math"}
#'    for very small/large ranges and \code{"plain"} otherwise. For
#'    \code{"math"}, tick positions are restricted to powers of 10.
#' @param show_regression_info **logical(1)** \cr If TRUE, show regression
#'    information (R2, slope, intercept) on the plot.
#' @param show_linear_range **logical(1)** \cr If TRUE, show the linear range
#'    of the calibration curve as a rectangle in the plot.
#' @param show_data_points **logical(1)** \cr If TRUE, show the data points on
#'    the plot.
#' @param plot_zeros **logical(1)** \cr If TRUE, show concentration 0 values
#'    in a separate "Blank" side panel for single-plot and all-in-one
#'    calibration curves.
#' @param blank_panel_ratio **numeric(1)** \cr Target width ratio of Blank
#'    panel relative to Calibration panel in linear x-axis mode when
#'    \code{plot_zeros = TRUE}. For example, \code{0.3} means the Blank panel
#'    is about 30% as wide as the Calibration panel when measured on the
#'    rendered plot (default is 0.1). Set to NULL to use automatic panel
#'    widths.
#' @param blank_panel_width_factor **numeric(1)** \cr Relative width of the
#'    Blank panel for linear x-axis plots when \code{plot_zeros = TRUE}.
#'    Interpreted as fraction of the positive calibration concentration range
#'    used to build the Blank panel x-span (default is 0.12).
#' @param plot_type **character(1)** \cr Type of plot for calibration curves:
#'    "single_plot" (generate one plot for one substance),
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
#' @param regression_info_size **numeric(1)** \cr Size of the regression
#'    information text on the plot, default is 3.
#' @param base_size **numeric(1)** \cr Base size for the plot theme,
#'    default is 11.
#'
#'
#' @importFrom checkmate assertCharacter assertFlag
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 geom_point scale_alpha_manual geom_line facet_wrap
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_rect geom_text guides guide_legend
#' @importFrom ggplot2  theme theme_bw xlab ylab unit
#' @importFrom scales label_comma
#'
#' @returns
#' List of a ggplot2 object containing the calibration curve plot (CC_plot) and
#' a data.frame containing summary data for annotations (annotation_dat).
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
#'     package = "CalibraCurve"
#' )
#' D_list <- readDataSE(file,
#'     concColName = "amount_fmol",
#'     substColName = "Substance", assayNumber = 1
#' )
#' CC_RES <- calc_single_curve(D_list[[1]], calcContinuousPrelimRanges = FALSE)
#'
#' plotCalibraCurve(list(substance = CC_RES))
plotCalibraCurve <- function(
        CC_RES, ylab = "Intensity", xlab = "Concentration",
    log10_x_axis = TRUE, log_tick_label_style = "auto",
        show_regression_info = FALSE, show_linear_range = TRUE,
    show_data_points = TRUE, plot_zeros = FALSE, blank_panel_ratio = 0.1,
        blank_panel_width_factor = 0.12, plot_type = "multiplot",
        point_colour = "black", curve_colour = "red",
        linear_range_colour = "black", multiplot_nrow = NULL,
        multiplot_ncol = NULL, multiplot_scales = "free",
        regression_info_size = 3, base_size = 11) {
    checkmate::assertCharacter(ylab, len = 1)
    checkmate::assertCharacter(xlab, len = 1)
    checkmate::assertFlag(log10_x_axis)
    checkmate::assertChoice(log_tick_label_style,
        choices = c("auto", "math", "plain")
    )
    checkmate::assertFlag(show_regression_info)
    checkmate::assertFlag(show_linear_range)
    checkmate::assertFlag(show_data_points)
    checkmate::assertFlag(plot_zeros)
    checkmate::assertNumber(blank_panel_ratio, lower = 0, upper = 1,
        null.ok = TRUE
    )
    checkmate::assertNumber(blank_panel_width_factor, lower = 0, upper = 1)
    checkmate::assertChoice(plot_type,
        choices = c("single_plot", "multiplot", "all_in_one")
    )
    if (plot_zeros && plot_type == "multiplot") {
        stop("plot_zeros = TRUE is not supported for plot_type = 'multiplot'. Use 'single_plot' or 'all_in_one'.")
    }
    checkmate::assertCharacter(point_colour, len = 1)
    checkmate::assertCharacter(curve_colour, len = 1)
    checkmate::assertCharacter(linear_range_colour, len = 1)
    concentration <- measurement <- final_linear_range <- predicted <- NULL
    LLOQ <- ULOQ <- eq <- substance <- NULL # silence notes when checking

    D_calib <- .prepareCalibData(CC_RES)
    annotation_dat <- .prepareAnnotationData(D_calib)
    curve_dat <- .prepareCurveData(CC_RES, D_calib)
    if (nrow(D_calib) == 0) {
        D_calib <- data.frame(
            substance = character(0),
            concentration = numeric(0),
            measurement = numeric(0),
            final_linear_range = factor(character(0), levels = c("No", "Yes"))
        )
    }

    D_plot <- D_calib
    use_blank_panel <- FALSE
    zero_x_display <- NA_real_
    blank_panel_span_dat <- NULL
    blank_panel_cutoff <- NA_real_
    blank_panel_upper_max <- NA_real_
    xintercept <- NULL

    if (plot_zeros && plot_type != "multiplot") {
        pos_idx <- is.finite(D_calib$concentration) & D_calib$concentration > 0
        zero_idx <- is.finite(D_calib$concentration) & D_calib$concentration == 0 &
            !is.na(D_calib$measurement)

        if (nrow(D_calib) > 0 && any(pos_idx) && any(zero_idx)) {
            D_plot <- D_calib[pos_idx | zero_idx, ]
            if (nrow(D_plot) > 0) {
                use_blank_panel <- TRUE
                zero_x_display <- min(D_calib$concentration[pos_idx], na.rm = TRUE) / 10
                pos_min <- min(D_calib$concentration[pos_idx], na.rm = TRUE)
                pos_max <- max(D_calib$concentration[pos_idx], na.rm = TRUE)
                blank_panel_cutoff <- (zero_x_display + pos_min) / 2

                D_plot$panel <- ifelse(D_plot$concentration == 0, "Blank", "Calibration")
                D_plot$concentration[D_plot$panel == "Blank"] <- zero_x_display
                D_plot$panel <- factor(D_plot$panel, levels = c("Blank", "Calibration"))

                if (nrow(curve_dat) > 0) {
                    curve_dat$panel <- factor("Calibration", levels = c("Blank", "Calibration"))
                }
                if (nrow(annotation_dat) > 0) {
                    annotation_dat$panel <- factor("Calibration", levels = c("Blank", "Calibration"))
                }

                y_ref <- min(D_plot$measurement[D_plot$measurement > 0], na.rm = TRUE)
                blank_x_span <- if (log10_x_axis) {
                    c(zero_x_display / 1.1, zero_x_display * 1.1)
                } else {
                    # Keep a visible Blank facet in linear-x mode with free_x spacing.
                    range_pos <- pos_max - pos_min
                    half_width <- max(zero_x_display * 0.1,
                        range_pos * blank_panel_width_factor
                    )
                    c(zero_x_display - half_width, zero_x_display + half_width)
                }

                # Dynamic upper bound for Blank facet detection. This tracks the
                # actual computed span and allows for ggplot panel expansion.
                blank_panel_upper_max <- blank_x_span[2] +
                    (blank_x_span[2] - blank_x_span[1]) * 0.25

                blank_panel_span_dat <- data.frame(
                    concentration = blank_x_span,
                    measurement = y_ref,
                    panel = factor("Blank", levels = c("Blank", "Calibration"))
                )
            }
        }
    }

    is_blank_panel_range <- function(x_vals) {
        if (!is.finite(blank_panel_cutoff) || !is.finite(blank_panel_upper_max) ||
            length(x_vals) != 2) {
            return(FALSE)
        }
        tol <- max(abs(c(x_vals, blank_panel_upper_max)), 1) * 1e-8
        x_vals[2] <= (blank_panel_upper_max + tol)
    }

    format_concentration_labels <- function(x_vals) {
        finite_vals <- x_vals[is.finite(x_vals)]
        non_zero_vals <- finite_vals[finite_vals != 0]
        if (length(non_zero_vals) > 0 &&
            (min(abs(non_zero_vals)) < 1e-3 || max(abs(non_zero_vals)) >= 1e5)) {
            return(scales::label_scientific(digits = 2)(x_vals))
        }
        scales::label_comma(drop0trailing = TRUE)(x_vals)
    }

    format_log10_labels_math <- function(x_vals) {
        lab <- rep(NA_character_, length(x_vals))
        pos_idx <- is.finite(x_vals) & x_vals > 0
        if (any(pos_idx)) {
            exp_vals <- log10(x_vals[pos_idx])
            exp_chr <- ifelse(
                abs(exp_vals - round(exp_vals)) < 1e-8,
                as.character(as.integer(round(exp_vals))),
                format(round(exp_vals, 2), trim = TRUE, scientific = FALSE)
            )
            lab[pos_idx] <- paste0("10^{", exp_chr, "}")
        }
        lab
    }

    use_math_log_labels <- function(x_vals) {
        if (log_tick_label_style == "math") {
            return(TRUE)
        }
        if (log_tick_label_style == "plain") {
            return(FALSE)
        }
        finite_vals <- x_vals[is.finite(x_vals) & x_vals > 0]
        if (length(finite_vals) == 0) {
            return(FALSE)
        }
        min(finite_vals) < 1e-3 || max(finite_vals) >= 1e5
    }

    format_log_labels <- function(x_vals) {
        if (use_math_log_labels(x_vals)) {
            lab <- format_log10_labels_math(x_vals)
            return(parse(text = ifelse(is.na(lab), "NA", lab)))
        }
        format_concentration_labels(x_vals)
    }

    power10_breaks <- function(x_vals) {
        finite_vals <- x_vals[is.finite(x_vals) & x_vals > 0]
        if (length(finite_vals) == 0) {
            return(numeric(0))
        }
        exp_seq <- seq(
            floor(log10(min(finite_vals))),
            ceiling(log10(max(finite_vals))),
            by = 1
        )
        10^exp_seq
    }

    use_math_x_labels <- use_math_log_labels(D_plot$concentration)
    use_math_y_labels <- use_math_log_labels(D_plot$measurement)

    if (log10_x_axis && !use_blank_panel) {
        # Non-positive x values cannot be displayed on a log10 axis.
        D_plot <- D_plot[D_plot$concentration > 0, , drop = FALSE]
    }

    # Non-positive or missing y values cannot be displayed on a log10 y-axis.
    D_plot <- D_plot[is.finite(D_plot$measurement) & D_plot$measurement > 0,
        , drop = FALSE
    ]

    pl <- ggplot2::ggplot(D_plot, ggplot2::aes(
        x = concentration,
        y = measurement, alpha = final_linear_range
    ))

    pl <- pl + ggplot2::scale_y_continuous(
        trans = "log10",
        breaks = if (use_math_y_labels) power10_breaks else ggplot2::waiver(),
        labels = function(y) {
            format_log_labels(y)
        }
    )
    if (use_blank_panel) {
        if (log10_x_axis) {
            pl <- pl + ggplot2::scale_x_continuous(
                trans = "log10",
                breaks = function(x) {
                    if (is_blank_panel_range(x)) {
                        return(zero_x_display)
                    }
                    if (use_math_x_labels) {
                        return(power10_breaks(x))
                    }
                    return(scales::log_breaks()(x))
                },
                labels = function(x) {
                    tol <- max(abs(c(x, zero_x_display)), na.rm = TRUE) * 1e-8
                    if (use_math_log_labels(x)) {
                        lab <- format_log10_labels_math(x)
                        lab[abs(x - zero_x_display) < tol] <- "0"
                        return(parse(text = ifelse(is.na(lab), "NA", lab)))
                    }
                    zero_lab <- ifelse(abs(x - zero_x_display) < tol,
                        "0", NA_character_
                    )
                    base_lab <- format_concentration_labels(x)
                    ifelse(is.na(zero_lab), base_lab, zero_lab)
                },
                minor_breaks = NULL
            )
        } else {
            pl <- pl + ggplot2::scale_x_continuous(
                breaks = function(x) {
                    if (is_blank_panel_range(x)) {
                        return(zero_x_display)
                    }
                    return(scales::breaks_pretty()(x))
                },
                labels = function(x) {
                    tol <- max(abs(c(x, zero_x_display)), na.rm = TRUE) * 1e-8
                    zero_lab <- ifelse(abs(x - zero_x_display) < tol, "0", NA_character_)
                    base_lab <- format_concentration_labels(x)
                    ifelse(is.na(zero_lab), base_lab, zero_lab)
                },
                minor_breaks = NULL
            )
        }
    } else {
        if (log10_x_axis) {
            pl <- pl + ggplot2::scale_x_continuous(
                trans = "log10",
                breaks = if (use_math_x_labels) power10_breaks else ggplot2::waiver(),
                labels = function(x) {
                    format_log_labels(x)
                }
            )
        } else {
            pl <- pl + ggplot2::scale_x_continuous(
                labels = format_concentration_labels
            )
        }
    }

    if (plot_type == "single_plot") {
        if (show_data_points) {
            pl <- pl +
                ggplot2::geom_point(
                    size = 1.7, color = point_colour,
                    show.legend = c(alpha = TRUE, colour = FALSE)
                ) +
                ggplot2::scale_alpha_manual(
                    values = c("Yes" = 1, "No" = 0.1),
                    name = "Linear range", drop = FALSE
                )
        }
        if (nrow(curve_dat) > 0) {
            pl <- pl + ggplot2::geom_line(
                color = curve_colour, data = curve_dat,
                ggplot2::aes(x = concentration, y = predicted),
                inherit.aes = FALSE, show.legend = FALSE
            )
        }

        if (use_blank_panel) {
            pl <- pl + ggplot2::facet_grid(
                . ~ panel,
                scales = "free_x",
                space = "free_x"
            )

            if (!is.null(blank_panel_span_dat)) {
                pl <- pl + ggplot2::geom_blank(
                    data = blank_panel_span_dat,
                    ggplot2::aes(x = concentration, y = measurement),
                    inherit.aes = FALSE
                )
            }

            pl <- pl + ggplot2::geom_vline(
                data = data.frame(
                    panel = factor("Blank", levels = c("Blank", "Calibration")),
                    xintercept = zero_x_display
                ),
                ggplot2::aes(xintercept = xintercept),
                inherit.aes = FALSE,
                color = "grey80",
                linewidth = 0.3
            )
        }

        if (show_linear_range) {
            pl <- pl + ggplot2::geom_rect(
                data = annotation_dat,
                ggplot2::aes(xmin = LLOQ, xmax = ULOQ),
                ymin = -Inf,
                ymax = Inf,
                alpha = 0.1, inherit.aes = FALSE, fill = linear_range_colour
            )
        }
        if (show_regression_info) {
            pl <- pl + ggplot2::geom_text(
                data = annotation_dat,
                mapping = ggplot2::aes(y = Inf, x = 0, label = eq),
                vjust = 1.4, hjust = -0.07, alpha = 0.6, inherit.aes = FALSE,
                color = curve_colour, size = regression_info_size
            )
        }
    }

    if (plot_type == "all_in_one") {
        if (show_data_points) {
            pl <- pl +
                ggplot2::geom_point(
                    size = 1.7,
                    ggplot2::aes(group = substance, colour = substance),
                    show.legend = c(alpha = TRUE, colour = FALSE)
                ) +
                ggplot2::scale_alpha_manual(
                    values = c("Yes" = 1, "No" = 0.1),
                    name = "Linear range", drop = FALSE
                )
        }
        if (nrow(curve_dat) > 0) {
            pl <- pl +
                ggplot2::geom_line(
                    data = curve_dat,
                    ggplot2::aes(
                        x = concentration, y = predicted,
                        color = substance, group = substance
                    ),
                    inherit.aes = FALSE,
                    show.legend = c(colour = TRUE, alpha = FALSE)
                )
        }

        if (use_blank_panel) {
            pl <- pl + ggplot2::facet_grid(
                . ~ panel,
                scales = "free_x",
                space = "free_x"
            )

            if (!is.null(blank_panel_span_dat)) {
                pl <- pl + ggplot2::geom_blank(
                    data = blank_panel_span_dat,
                    ggplot2::aes(x = concentration, y = measurement),
                    inherit.aes = FALSE
                )
            }

            pl <- pl + ggplot2::geom_vline(
                data = data.frame(
                    panel = factor("Blank", levels = c("Blank", "Calibration")),
                    xintercept = zero_x_display
                ),
                ggplot2::aes(xintercept = xintercept),
                inherit.aes = FALSE,
                color = "grey80",
                linewidth = 0.3
            )
        }
    }

    if (plot_type == "multiplot") {
        if (show_data_points) {
            pl <- pl +
                ggplot2::geom_point(
                    size = 1.7, color = point_colour,
                    show.legend = c(alpha = TRUE, colour = FALSE)
                ) +
                ggplot2::scale_alpha_manual(
                    values = c("Yes" = 1, "No" = 0.1),
                    name = "Linear range", drop = FALSE
                )
        }
        if (nrow(curve_dat) > 0) {
            pl <- pl + ggplot2::geom_line(
                color = curve_colour, data = curve_dat,
                ggplot2::aes(x = concentration, y = predicted),
                inherit.aes = FALSE, show.legend = FALSE
            )
        }
        pl <- pl + ggplot2::facet_wrap(substance ~ .,
            scales = multiplot_scales,
            nrow = multiplot_nrow, ncol = multiplot_ncol
        )
        if (show_linear_range) {
            pl <- pl + ggplot2::geom_rect(
                data = annotation_dat,
                ggplot2::aes(xmin = LLOQ, xmax = ULOQ),
                ymin = -Inf,
                ymax = Inf,
                alpha = 0.1, inherit.aes = FALSE, fill = linear_range_colour
            )
        }
        if (show_regression_info) {
            pl <- pl + ggplot2::geom_text(
                data = annotation_dat,
                mapping = ggplot2::aes(y = Inf, x = 0, label = eq),
                vjust = 1.4, hjust = -0.07, alpha = 0.6, inherit.aes = FALSE,
                color = curve_colour, size = regression_info_size
            )
        }
    }

    pl <- pl +
        ggplot2::guides(
            alpha = ggplot2::guide_legend(
                title = "Linear range",
                reverse = TRUE, order = 2
            ),
            colour = ggplot2::guide_legend(title = "Substance", order = 1)
        ) +
        ggplot2::theme_bw(base_size = base_size) + ggplot2::ylab(ylab) +
        ggplot2::xlab(xlab) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm"))

    if (use_blank_panel && !is.null(blank_panel_ratio)) {
        pl <- pl + ggplot2::theme(
            panel.widths = grid::unit(c(blank_panel_ratio, 1), "null")
        )
    }

    return(list(CC_plot = pl, annotation_dat = annotation_dat))
}







#' Prepare calibration data
#'
#' @param CC_RES **list** Result object of \code{\link{CalibraCurve}}.
#'
#' @returns Dataframe with calibration data for plotting.
#' @keywords internal
.prepareCalibData <- function(CC_RES) {
    ## remove empty results from list (if calculation stopped)
    is_RES_null <- vapply(CC_RES, is.null, logical(1))
    CC_RES <- CC_RES[!is_RES_null]
    # X: result of 1 CalibraCurveRun
    f <- function(X) {
        if (is.null(X$result_table_obs)) {
            return(data.frame(
                substance = character(0),
                concentration = numeric(0),
                measurement = numeric(0),
                removed_while_cleaning = logical(0),
                percentage_bias = numeric(0),
                response_factor = numeric(0),
                RF_within_thres = logical(0),
                final_linear_range = logical(0),
                intercept = numeric(0),
                coeff = numeric(0),
                r2 = numeric(0),
                weight_m = character(0)
            ))
        }
        has_model <- !is.null(X$mod) && inherits(X$mod, "lm")
        tmp <- data.frame(X$result_table_obs,
            intercept = if (has_model) unname(X$mod$coefficients[1]) else NA_real_,
            coeff = if (has_model) unname(X$mod$coefficients[2]) else NA_real_,
            r2 = if (has_model) summary(X$mod)$r.squared else NA_real_,
            weight_m = if (!is.null(X$weightingMethod)) X$weightingMethod else NA_character_
        )
        return(tmp)
    }
    if (length(CC_RES) == 0) {
        return(data.frame(
            substance = character(0),
            concentration = numeric(0),
            measurement = numeric(0),
            removed_while_cleaning = logical(0),
            percentage_bias = numeric(0),
            response_factor = numeric(0),
            RF_within_thres = logical(0),
            final_linear_range = logical(0),
            intercept = numeric(0),
            coeff = numeric(0),
            r2 = numeric(0),
            weight_m = character(0)
        ))
    }
    D_calib <- do.call(rbind, lapply(CC_RES, f))
    if (nrow(D_calib) == 0) {
        return(D_calib)
    }
    D_calib$measurement[D_calib$measurement == 0] <- NA # avoid log10(0)
    D_calib$final_linear_range <- factor(D_calib$final_linear_range,
        levels = c(FALSE, TRUE), labels = c("No", "Yes")
    )
    return(D_calib)
}


#' Prepare annotation data
#'
#' @param D_calib **data.frame** \cr Result of \code{\link{.prepareCalibData}}
#'
#' @returns Dataframe with annotation data for the plot.
#' @keywords internal
.prepareAnnotationData <- function(D_calib) {
    if (nrow(D_calib) == 0) {
        return(data.frame(
            substance = character(0),
            intercept = numeric(0),
            coeff = numeric(0),
            r2 = numeric(0),
            LLOQ = numeric(0),
            ULOQ = numeric(0),
            eq = character(0)
        ))
    }
    D_calib_split <- split(D_calib, D_calib$substance)

    f <- function(X) {
        intercept <- X$intercept[1]
        coeff <- X$coeff[1]
        r2 <- X$r2[1]
        eq <- if (is.na(intercept) || is.na(coeff) || is.na(r2)) {
            "Calibration failed"
        } else {
            paste0(
                "y = ", format(intercept, scientific = TRUE, digits = 2),
                " + ", format(coeff, scientific = TRUE, digits = 2), " * x",
                " (R2 = ", round(r2, 3), ")"
            )
        }

        yes_idx <- X$final_linear_range == "Yes"
        tmp <- data.frame(
            substance = X$substance[1],
            intercept = intercept,
            coeff = coeff,
            r2 = r2,
            LLOQ = if (any(yes_idx)) min(X$concentration[yes_idx]) else NA_real_,
            ULOQ = if (any(yes_idx)) max(X$concentration[yes_idx]) else NA_real_,
            eq = eq
        )
        return(tmp)
    }

    annotation_dat <- do.call(rbind, lapply(D_calib_split, f))
    return(annotation_dat)
}


#' Prepare curve data
#'
#' @param CC_RES **list** \cr Results of \code{\link{CalibraCurve}}. Each list
#'      element is the result of \code{\link{calc_single_curve}}
#' @param D_calib **data.frame** \cr Result of \code{\link{.prepareCalibData}}
#'
#' @returns Dataframe with data points on a grid to plot the calibration curve.
#'
#' @importFrom stats predict
#' @keywords internal
.prepareCurveData <- function(CC_RES, D_calib) {
    ## remove empty results from list (if calculation stopped)
    is_RES_null <- vapply(CC_RES, is.null, logical(1))
    CC_RES <- CC_RES[!is_RES_null]
    if (length(CC_RES) == 0) {
        return(data.frame(
            substance = character(0),
            concentration = numeric(0),
            predicted = numeric(0)
        ))
    }
    D_calib_split <- split(D_calib, D_calib$substance)
    substances <- names(D_calib_split)


    f <- function(X, Y, Z) { # X = CC_RES, Y = D_calib, Z = substances
        if (is.null(X$mod) || !inherits(X$mod, "lm")) {
            return(data.frame(
                substance = character(0),
                concentration = numeric(0),
                predicted = numeric(0)
            ))
        }
        conc_pos <- Y$concentration[is.finite(Y$concentration) & Y$concentration > 0]
        if (length(conc_pos) == 0) {
            return(data.frame(
                substance = character(0),
                concentration = numeric(0),
                predicted = numeric(0)
            ))
        }

        grid <- seq(log10(min(conc_pos)), log10(max(conc_pos)),
            length.out = 1000
        )
        pred <- stats::predict(X$mod, newdata =
                                data.frame(Concentration = 10^grid))
        curve_dat_tmp <- data.frame(
            substance = Z, concentration = 10^grid,
            predicted = pred
        )
        # remove negative values in prediction (causes problems in log10-trans)
        curve_dat_tmp <- curve_dat_tmp[curve_dat_tmp$predicted > 0, ]
    }

    curve_dat <- do.call(rbind, mapply(f, CC_RES, D_calib_split, substances,
        SIMPLIFY = FALSE
    ))
    return(curve_dat)
}




#' Plot response factors
#'
#' @param RES **list** \cr Results of \code{\link{calc_single_curve}}.
#' @param RfThresL **numeric(1)** \cr Lower threshold for response factor in
#'    percent (default is 80).
#' @param RfThresU **numeric(1)** \cr Upper threshold for response factor in
#'    percent (default is 120).
#' @param ylab **character(1)** \cr y-axis label.
#' @param xlab **character(1)** \cr x-axis label.
#' @param log10_x_axis **logical(1)** \cr If TRUE (default), use a log10
#'    transformation for the x-axis. If FALSE, keep original concentration
#'    values on a linear x-axis.
#' @param colour_threshold **character(1)** \cr Colour for horizontal threshold
#'    lines, default is "orange".
#' @param colour_within **character(1)** \cr Colour for points and lines within
#'    the threshold, default is "#00BFC4" (default ggplot colour).
#' @param colour_outside **character(1)** \cr Colour for horizontal outside of
#'    the threshold, default is "#F8766D" (default ggplot colour).
#' @param legend **logical(1)** \cr If TRUE, a legend is added to the plot,
#'      default is FALSE.
#' @param base_size **numeric(1)** \cr Base size for the plot theme,
#'      default is 11.
#' @param plot_type **character(1)** \cr Type of response factor plot:
#'      "single_plot" (default, one plot for one substance),
#'      "multiplot" (facetted plot for multiple substances)
#'      or "all_in_one" (all substances in one panel).
#' @param multiplot_nrow **integer(1)** \cr Number of rows for the multiplot
#'      layout (default is NULL, determined automatically).
#' @param multiplot_ncol **integer(1)** \cr Number of columns for the
#'      multiplot layout (default is NULL, determined automatically).
#' @param multiplot_scales **character(1)** \cr Scales for the multiplot
#'      layout, default is "free". Other options are "fixed", "free_x",
#'      "free_y".
#' @param substance **character(1)** \cr Name of the substance (default is
#'      "substance1"). Will be added to the result files and will be used as
#'      labels in the plots.
#'
#' @returns
#' A ggplot2 object containing the response factor plot.
#'
#' @importFrom checkmate assertNumeric assertCharacter
#' @importFrom ggplot2 ggplot scale_x_continuous geom_point aes geom_line
#' @importFrom ggplot2 scale_alpha_manual scale_colour_manual scale_fill_manual
#' @importFrom ggplot2 geom_hline theme theme_bw xlab ylab unit
#' @importFrom scales comma
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
#'     package = "CalibraCurve"
#' )
#' D_list <- readDataSE(file,
#'     concColName = "amount_fmol",
#'     substColName = "Substance", assayNumber = 1
#' )
#' CC_RES <- calc_single_curve(D_list[[1]], calcContinuousPrelimRanges = FALSE)
#'
#' plotResponseFactors(RES = CC_RES)
plotResponseFactors <- function(
        RES, RfThresL = 80, RfThresU = 120,
        ylab = "Response Factor", xlab = "Concentration",
    log10_x_axis = TRUE,
        colour_threshold = "orange", colour_within = "#00BFC4",
        colour_outside = "#F8766D", legend = TRUE, base_size = 11,
    plot_type = "single_plot", multiplot_nrow = NULL,
    multiplot_ncol = NULL, multiplot_scales = "free",
        substance = NULL) {
    checkmate::assertNumeric(RfThresL, lower = 0, upper = 100, finite = TRUE)
    checkmate::assertNumeric(RfThresU, lower = 100)
    checkmate::assertCharacter(ylab, len = 1)
    checkmate::assertCharacter(xlab, len = 1)
    checkmate::assertFlag(log10_x_axis)
    checkmate::assertCharacter(colour_threshold, len = 1)
    checkmate::assertCharacter(colour_within, len = 1)
    checkmate::assertCharacter(colour_outside, len = 1)
    checkmate::assertChoice(plot_type,
        choices = c("single_plot", "multiplot", "all_in_one")
    )
    checkmate::assertChoice(multiplot_scales,
        choices = c("fixed", "free", "free_x", "free_y")
    )
    if (!is.null(multiplot_nrow)) {
        checkmate::assertInt(multiplot_nrow, lower = 1)
    }
    if (!is.null(multiplot_ncol)) {
        checkmate::assertInt(multiplot_ncol, lower = 1)
    }

    # silence notes when checking the package
    concentration <- response_factor <- final_linear_range <- NULL
    mean_response_factor <- RF_within_thres <- NULL
    substance_rf <- rf_upper <- rf_lower <- threshold_y <- threshold_type <- NULL

    if (plot_type != "single_plot") {
        checkmate::assertList(RES, min.len = 1)
        all_names <- names(RES)
        if (is.null(all_names)) {
            all_names <- rep(NA_character_, length(RES))
        }

        req_obs_cols <- c(
            "concentration", "response_factor",
            "final_linear_range", "RF_within_thres"
        )
        req_sum_cols <- c(
            "concentration", "mean_response_factor", "RF_within_thres"
        )

        .norm_rf_cols <- function(dat, req_cols) {
            if (is.null(dat)) {
                dat <- data.frame()
            }
            dat <- as.data.frame(dat)
            miss_cols <- setdiff(req_cols, names(dat))
            if (length(miss_cols) > 0) {
                for (col_nm in miss_cols) {
                    dat[[col_nm]] <- NA
                }
            }
            dat[, req_cols, drop = FALSE]
        }

        .rf_single <- function(X, nm, idx) {
            dat_obs <- .norm_rf_cols(X$result_table_obs, req_obs_cols)
            dat_obs <- dat_obs[!is.na(dat_obs$response_factor), , drop = FALSE]
            dat_obs$substance_rf <- ifelse(
                !is.na(nm) && nzchar(nm), nm, paste0("substance_", idx)
            )

            dat_sum <- .norm_rf_cols(X$result_table_conc_levels, req_sum_cols)
            dat_sum <- dat_sum[!is.na(dat_sum$mean_response_factor), , drop = FALSE]
            dat_sum$substance_rf <- ifelse(
                !is.na(nm) && nzchar(nm), nm, paste0("substance_", idx)
            )

            list(obs = dat_obs, sum = dat_sum)
        }

        rf_parts <- mapply(.rf_single, RES, all_names, seq_along(RES),
            SIMPLIFY = FALSE
        )
        range_dat <- do.call(rbind, lapply(rf_parts, `[[`, "obs"))
        sum_dat <- do.call(rbind, lapply(rf_parts, `[[`, "sum"))
    } else {
        range_dat <- RES$result_table_obs
        range_dat <- range_dat[!is.na(range_dat$response_factor), ]
        sum_dat <- RES$result_table_conc_levels
        sum_dat <- sum_dat[!is.na(sum_dat$mean_response_factor), ]
        range_dat$substance_rf <- ifelse(
            is.null(substance) || !nzchar(substance), "substance1", substance
        )
        sum_dat$substance_rf <- ifelse(
            is.null(substance) || !nzchar(substance), "substance1", substance
        )
    }

    threshold_dat <- split(range_dat, range_dat$substance_rf)
    threshold_dat <- lapply(names(threshold_dat), function(nm) {
        dat <- threshold_dat[[nm]]
        data.frame(
            substance_rf = nm,
            rf_mean = mean(dat$response_factor[dat$final_linear_range])
        )
    })
    threshold_dat <- do.call(rbind, threshold_dat)
    threshold_dat$rf_upper <- threshold_dat$rf_mean * (RfThresU / 100)
    threshold_dat$rf_lower <- threshold_dat$rf_mean * (RfThresL / 100)
    threshold_dat <- threshold_dat[is.finite(threshold_dat$rf_mean), ]
    threshold_lines <- rbind(
        data.frame(
            substance_rf = threshold_dat$substance_rf,
            threshold_y = threshold_dat$rf_upper,
            threshold_type = "Threshold range"
        ),
        data.frame(
            substance_rf = threshold_dat$substance_rf,
            threshold_y = threshold_dat$rf_lower,
            threshold_type = "Threshold range"
        )
    )

    ### initialize plot
    pl <- ggplot2::ggplot(range_dat, mapping = ggplot2::aes(
        x = concentration,
        y = response_factor, color = RF_within_thres,
        fill = RF_within_thres, alpha = RF_within_thres
    ))

    if (plot_type == "single_plot") {
        pl <- pl + ggplot2::aes(group = 1)
    } else if (plot_type == "all_in_one") {
        pl <- pl + ggplot2::aes(group = interaction(substance_rf, RF_within_thres))
    } else {
        pl <- pl + ggplot2::aes(group = 1)
    }

    if (log10_x_axis) {
        pl <- pl + ggplot2::scale_x_continuous(
            trans = "log10",
            labels = scales::comma
        )
    } else {
        pl <- pl + ggplot2::scale_x_continuous(labels = scales::comma)
    }

    ### add data points + mean response factors per concentration level
    if (plot_type != "all_in_one") {
        pl <- pl + ggplot2::geom_point(size = 1.7, shape = 21)
    }
    pl <- pl + ggplot2::geom_point(
        data = sum_dat,
        mapping = ggplot2::aes(x = concentration, y = mean_response_factor),
        color = "black", shape = 21, size = 2.5
    )

    ### adjust colour of lines (only within_colour for consecutive mean response
    ### factors within the threshold)
    sum_dat2 <- sum_dat
    for (i in seq_len(nrow(sum_dat) - 1)) {
        ### if the mean RF is within the threshold and the next one is too:
        if (sum_dat$RF_within_thres[i] == TRUE &
            sum_dat$RF_within_thres[i + 1] == TRUE) {
            sum_dat2$RF_within_thres[i] <- TRUE
        } else {
            sum_dat2$RF_within_thres[i] <- FALSE
        }
    }

    ### add line between mean response factors
    pl <- pl + ggplot2::geom_line(
        data = sum_dat2,
        mapping = ggplot2::aes(
            x = concentration,
            y = mean_response_factor,
            group = if (plot_type == "all_in_one") substance_rf else 1
        )
    )

    ## scaling of alpha and colours
    pl <- pl +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
        ggplot2::scale_colour_manual(values = c(
            "TRUE" = colour_within,
            "FALSE" = colour_outside
        )) +
        ggplot2::scale_fill_manual(values = c(
            "TRUE" = colour_within,
            "FALSE" = colour_outside
        ))

    if (plot_type != "all_in_one") {
        pl <- pl + ggplot2::scale_linetype_manual(values = c(
            "Threshold range" = "dashed"
        ))
    }

    ### add horizontal lines for response factor thresholds
    if (plot_type != "all_in_one") {
        pl <- pl +
            ggplot2::geom_hline(
                data = threshold_lines,
                mapping = ggplot2::aes(
                    yintercept = threshold_y,
                    linetype = threshold_type
                ),
                color = colour_threshold,
                inherit.aes = FALSE
            )
    }

    ### theme and axis limits
    pl <- pl + ggplot2::theme_bw(base_size = base_size) + ggplot2::ylab(ylab) +
        ggplot2::xlab(xlab)

    if (plot_type == "multiplot") {
        pl <- pl + ggplot2::facet_wrap(substance_rf ~ ., scales = multiplot_scales,
            nrow = multiplot_nrow, ncol = multiplot_ncol
        )
    }

    if (!legend) {
        pl <- pl + ggplot2::theme(
            legend.position = "none",
            plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm")
        )
    } else {
        legend_name <- "RF within threshold?"
        pl <- pl + ggplot2::theme(
            plot.margin = ggplot2::unit(c(0.5, 0.7, 0.5, 0.5), "cm")
        ) +
            ggplot2::guides(
                color = ggplot2::guide_legend(title = legend_name, order = 1),
                fill = ggplot2::guide_legend(title = legend_name, order = 1),
                alpha = ggplot2::guide_legend(title = legend_name, order = 1),
                linetype = ggplot2::guide_legend(title = "Threshold", order = 2)
            )
    }

    if (plot_type == "single_plot") {
        pl <- pl + ggplot2::ggtitle(substance)
    }
    return(pl)
}


