



#' Predict concentrations from intensity values based on a given calibration curve
#'
#' @param CC_res **list** \cr Results of \code{\link{CalibraCurve}}.
#' @param newdata **numeric** \cr A vector of intensity values for which to predict concentrations.
#'
#' @returns A data frame with the following columns:
#' - `intensity`: The input intensity values.
#' - `predicted_concentrations`: The predicted concentrations based on the calibration curve.
#' - `linear_range`: A logical vector indicating whether the predicted concentrations are within the final linear range.
#' @export
#'
#' @details
#' The function will give a warning if any of the predicted concentrations are outside the final linear range.
#' This is important to ensure that the predictions are reliable and within the linear range of the calibration curve.
#'
#' @examples
predictConcentration <- function(CC_res, newdata) {

  mod <- CC_res$mod

  FLR <- CC_res$final_linear_range
  min_FLR <- min(FLR)
  max_FLR <- max(FLR)

  coeffs <- stats::coefficients(mod)
  intercept <- coeffs[1]
  slope <- coeffs[2]

  # linear model: intensity = intercept + slope * concentration + e
  # prediction of concentration: concentration = (intensity - intercept) / slope
  predictedConcentration = (intercept + slope * newdata) / slope

  # check if predicted concentration is within the final linear range
  linear_range <- predictedConcentration >= min_FLR & predictedConcentration <= max_FLR
  if (any(!linear_range)) {
    warning("At least one predicted concentration is outside the final linear range. Results may be unreliable.")
  }

  res <- data.frame(intensity = newdata, predicted_concentrations = predictedConcentration,
                    linear_range = linear_range)

  return(res)

}

