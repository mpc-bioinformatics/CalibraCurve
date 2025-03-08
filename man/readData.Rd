% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_data_preprocessing.R
\name{readData}
\alias{readData}
\title{Read in data in different input formats}
\usage{
readData(
  path,
  filetype,
  conc_col,
  meas_col,
  sep = ",",
  dec = ";",
  header = TRUE,
  na.strings = c("NA", "NaN", "Filtered", "#NV"),
  sheet = 1
)
}
\arguments{
\item{path}{\strong{character(1)} \cr path to the file}

\item{filetype}{\strong{character(1)} \cr csv or txt or xlsx}

\item{conc_col}{\strong{integer(1)} \cr column number of concentration levels}

\item{sep}{\strong{character(1)} \cr separator}

\item{dec}{\strong{character(1)} \cr decimal separator}

\item{header}{\strong{logical(1)} \cr if TRUE, first line is counted as a heading}

\item{na.strings}{\strong{character} \cr character vector of strings which are to be interpreted as NA}

\item{sheet}{number of the sheet (only needed for xlsx files, default is to use the first sheet)}

\item{meas_cols}{\strong{integer} \cr column number(s) of measurements}
}
\description{
Read in data in different input formats
}
