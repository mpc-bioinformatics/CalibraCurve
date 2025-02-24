



#' Read in data in different input formats
#'
#' @param path \strong{character(1)} \cr path to the file
#' @param filetype \strong{character(1)} \cr csv or txt or xlsx
#' @param conc_col \strong{integer(1)} \cr column number of concentration levels
#' @param meas_cols \strong{integer} \cr column number(s) of measurements
#' @param sep \strong{character(1)} \cr separator
#' @param dec \strong{character(1)} \cr decimal separator
#' @param header \strong{logical(1)} \cr if TRUE, first line is counted as a heading
#' @param na.strings \strong{character} \cr character vector of strings which are to be interpreted as NA
#' @param sheet number of the sheet (only needed for xlsx files, default is to use the first sheet)
#'
#' @returns
#' @export
#'
#' @examples
readData <- function(path,
                     filetype,
                     conc_col,
                     meas_cols,
                     sep = ",",
                     dec = ";",
                     header = TRUE,
                     na.strings = c("NA", "NaN", "Filtered", "#NV"),
                     sheet = 1) {
  ### TODO: error messages if data doesn't fit to the expected format

  if (filetype == "csv" | filetype == "txt") {
    rawData <- read.table(path,
                          sep = sep,
                          header = header,
                          dec = dec)
  }
  if (filetype == "xlsx") {
    rawData <- openxlsx::read.xlsx(path, colNames = header, sheet = sheet)
  }


  ### extract relevant columns:
  rawData <- rawData[, c(conc_col, meas_cols)]
  colnames(rawData) <- c("Concentration", paste("Measurement", 1:length(meas_cols), sep = "_"))

  return(rawData)
}




#' Data preprocessing: Select all rows with identical concentrations
#'
#' @param x concentration level to check
#' @param rawData data set to be filtered
#'
#' @returns
#' @export
#'
#' @examples
selConcentrationLevels <- function(x, rawData) {
  result <- rawData[rawData$Concentration == x, ]
  return(result)
}



#' Data preprocessing: Check for sufficient number of replicates
#'
#' @param x concentration level to check
#' @param data data set to be filtered
#' @param minNumber minimal number of data points per concentration level
#'
#' @returns
#' @export
#'
#' @examples
checkNumberReplicates <- function(x, data, minNumber) {
  if (nrow(data[[x]]) < minNumber) {
    result <- FALSE
  }
  else {
    result <- TRUE
  }
  return(result)
}



#' Clean data
#'
#' @param rawData data set to be cleaned, result of readData.
#' @param min_replicates minimal number of data points per concentration level
#'
#' @returns
#' @export
#'
#' @examples
cleanData <- function(rawData, min_replicates) {
  ### TODO: make possible to keep more than one measurement column
  ### TODO: remove 0 values (problems later with log-transformation?)

  ### remove rows with unknown concentration
  dataValidated <- rawData[!is.na(rawData$Concentration),]

  ### remove rows with unknown measurement
  dataValidated <<- rawData[!is.na(rawData$Measurement_1),]

  # Determination of existing concentration levels in the validated data
  concLevels <- unique(dataValidated$Concentration)
  concLevels <- sort(concLevels, decreasing = FALSE)

  # Transforming a data set into a list with entries for each concentration level (and the related data)
  dataValidated <- lapply(concLevels, FUN=selConcentrationLevels, rawData=dataValidated)

  # Deleting concentration levels with insufficient number of replicates
  dataValidated <- dataValidated[sapply(1:length(dataValidated), FUN = checkNumberReplicates, data = dataValidated, minNumber = min_replicates)]
  ### TODO: keep the data to show them in the plot, maybe as grey/faded data points

  return(dataValidated)
}

