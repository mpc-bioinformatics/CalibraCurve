% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_final_linear_range.R
\name{calculate_FLR}
\alias{calculate_FLR}
\title{Calculate the final linear range}
\usage{
calculate_FLR(
  dataPrelim,
  weightingMethod = "1/x^2",
  centralTendencyMeasure = "mean",
  perBiasThres = 20,
  considerPerBiasCV = TRUE,
  perBiasDistThres = 10
)
}
\arguments{
\item{dataPrelim}{data frame with concentration levels and corresponding response values (result from calculate_PLR)}

\item{weightingMethod}{method for weighting (currently "1/x", "1/x^2" and "None" are supported, default is 1/x^2)}

\item{centralTendencyMeasure}{"mean" or "median" (for calculating average percent bias), default is "mean"}

\item{perBiasThres}{threshold for average percent bias, default is 20\%}

\item{considerPerBiasCV}{consider CV for the selection of the concentration level, default is TRUE}

\item{perBiasDistThres}{threshold for the difference in average percent bias (for lower differences, CV will be considered), default is 10\%}
}
\description{
Calculate the final linear range
}
