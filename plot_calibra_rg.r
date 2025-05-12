library(tidyverse)
library(ggrepel)
library(scales)

## install.packages("/home/grugelro/Documents/kleinkram/calibracurve_dev/CalibraCurve-master",
##                  repos = NULL,
##                  type = "source")
# library(CalibraCurve)
# calibres <- CalibraCurve(path = "/home/grugelro/Documents/kleinkram/calibracurve_dev/Clonidin_0224.xlsx",
#              conc_col = 1,
#              meas_col = 3)


# plotConcVar <- function(calib_list = calibres) {
#   range_dat <- calib_list$result_table_obs
#
#   sum_dat <- range_dat %>%
#     group_by(substance, concentration) %>%
#     reframe(sigma = sd(measurement),
#             sigma2 = sigma^2) %>%
#     pivot_longer(cols = c(sigma, sigma2))
#
#   ggplot(sum_dat, aes(x = concentration, y = value, col = name)) +
#     geom_point() +
#     geom_line() + theme_bw()
#
# }
# plotConcVar()

plotCalibCurve <- function(calib_list = calibres,
                           ylab = "log10(Peak area)",
                           xlab = "log10(Concentration[mg/L])",
                           reg_eq = FALSE) {
  range_dat <- calib_list$result_table_obs
  final_mod <- calib_list$mod
  weight_m <- calib_list$weightingMethod
  
  range_dat2 <- range_dat %>%
    filter(final_linear_range) %>%
    group_by(substance) %>%
    reframe(predd = list(final_mod$coefficients),
            r2 = summary(final_mod)$r.squared)
  
  subst_range <- range_dat %>%
    group_by(substance) %>%
    summarise(mincon = min(concentration),
              maxcon = max(concentration))
  
  range_dat2 <- range_dat2 %>%
    unnest(predd) %>%
    mutate(partype = rep(c("intercept", "coeff"), n() / 2)) %>%
    pivot_wider(names_from = partype, values_from = predd)
  
  eq_dat <- range_dat2 %>%
    mutate(
      eq = paste0(
        "y = ",
        format(intercept, scientific = TRUE, digits = 2),
        " + ",
        format(coeff, scientific = TRUE, digits = 2),
        " * x",
        " (R² = ",
        round(r2, 3),
        ")"
      ),
      r2_eq = paste0("R^2 = ", round(r2, 4))
    )
  
  range_dat3 <- left_join(range_dat, range_dat2) %>%
    mutate(predicted = intercept + coeff * concentration)
  
  an_dat <- range_dat3 %>% group_by(substance) %>%
    summarise(LLOQ = min(concentration[final_linear_range]),
              ULOQ = max(concentration[final_linear_range]))
  
  range_dat4 <- left_join(range_dat2, subst_range) %>%
    rowwise() %>%
    mutate(xo = list(10^(seq(
      log10(mincon), log10(maxcon), length.out = 1000
    )))) %>%
    ungroup() %>%
    unnest(xo) %>%
    mutate(predicted = intercept + coeff * xo)
  
  res_plot <- ggplot(range_dat3, aes(x = log10(concentration), log10(measurement))) +
    geom_line(
      color = "red",
      data = range_dat4,
      aes(x = log10(xo), y = log10(predicted)),
      inherit.aes = FALSE
    ) +
    geom_point(size = 1.7, aes(alpha = final_linear_range)) + theme_bw() +
    geom_rect(
      data = an_dat,
      aes(
        xmin = log10(LLOQ),
        xmax = log10(ULOQ),
        ymin = -Inf,
        ymax = Inf
      ),
      alpha = 0.1,
      inherit.aes = FALSE
    ) +
    theme(legend.position = "none") +
    facet_wrap(substance ~ ., scales = "free") +
    scale_x_continuous(
      labels = function(x)
        sub(
          "\\.?0+$",
          "",
          format(10^x, scientific = FALSE, zero.print = FALSE)
        )
    ) +
    scale_y_continuous(
      labels = function(y)
        paste0(round(10^y, 5))
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1)) +
    ylab(ylab) +
    xlab(xlab)
  
  if (reg_eq) {
    res_plot <- res_plot + geom_text(
      data = eq_dat,
      aes(-Inf, Inf, label = eq),
      vjust = 3,
      hjust = -0.1
    )
  }
  
  res_plot
}
# plotCalibCurve()

plotRespFactor <- function(calib_list = calibres,
                           ylab = "Response Factor",
                           xlab = "log10(Concentration[mg/L])") {
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
# plotRespFactor()
