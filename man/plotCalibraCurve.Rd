% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_plotting.R
\name{plotCalibraCurve}
\alias{plotCalibraCurve}
\title{Plot the calibration curve}
\usage{
plotCalibraCurve(
  RES,
  ylab = "Intensity",
  xlab = "Concentration",
  show_regression_info = FALSE,
  show_linear_range = TRUE,
  show_data_points = TRUE,
  point_colour = "black",
  curve_colour = "red",
  linear_range_colour = "black"
)
}
\arguments{
\item{RES}{\strong{list} \cr Results of \code{\link{CalibraCurve}}.}

\item{ylab}{\strong{character(1)} \cr y-axis label.}

\item{xlab}{\strong{character(1)} \cr x-axis label.}

\item{show_regression_info}{\strong{logical(1)} \cr If TRUE, show regression information (R2, slope, intercept) on the plot.}

\item{show_linear_range}{\strong{logical(1)} \cr If TRUE, show the linear range of the calibration curve as a rectangle in the plot.}

\item{show_data_points}{\strong{logical(1)} \cr If TRUE, show the data points on the plot.}

\item{point_colour}{\strong{character(1)} \cr Colour of the data points, default is "black".}

\item{curve_colour}{\strong{character(1)} \cr Colour of the calibration curve, default is "red".}

\item{linear_range_colour}{\strong{character(1)} \cr Colour of the linear range background, default is "black" (colour is weakened by alpha = 0.1).}
}
\value{
A ggplot2 object containing the calibration curve plot.
}
\description{
Plot the calibration curve
}
\examples{

data(RES_ALB)
plotCalibraCurve(RES_ALB)
}
