
#' Read in data in different input formats
#'
#' @param data_path **character(1)** \cr Path to the data file (.csv, .txt or .xlsx file).
#' @param filetype **character(1)** \cr Type of input file: "csv" or "txt" or "xlsx".
#' @param conc_col **integer(1)** \cr Column number of the concentration values.
#' @param meas_cols **integer** \cr Column number of the concentration values.
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\t" for tab.
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot.
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param na.strings **character** \cr Character vector of strings which are to be interpreted as NA.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
#'
#' @returns data.frame with two numeric columns: Concentration and Measurement
#' @export
#'
#' @examples
readData <- function(data_path,
                     filetype,
                     conc_col,
                     meas_col,
                     sep = ",",
                     dec = ";",
                     header = TRUE,
                     na.strings = c("NA", "NaN", "Filtered", "#NV"),
                     sheet = 1) {

  ### check input arguments
  checkmate::assert_character(data_path, len = 1)
  checkmate::assert_choice(filetype, c("csv", "txt", "xlsx"))
  checkmate::assert_int(conc_col)
  checkmate::assert_int(meas_col)
  checkmate::assert_character(sep, len = 1)
  checkmate::assert_character(dec, len = 1)
  checkmate::assert_flag(header)
  checkmate::assert_character(na.strings)
  checkmate::assert_int(sheet)
  if (conc_col == meas_col) stop("Concentration and measurement columns cannot be identical.")



  if (filetype == "csv" | filetype == "txt") {
    rawData <- read.table(data_path,
                          sep = sep,
                          header = header,
                          dec = dec)
  }
  if (filetype == "xlsx") {
    rawData <- openxlsx::read.xlsx(data_path, colNames = header, sheet = sheet)
  }

  ### check if column numbers are valid
  if (meas_col > ncol(rawData)) {
    stop("Number of measurement column cannot be larger than number of columns in data set.")
  }
  if (conc_col > ncol(rawData)) {
    stop("Number of concentration column cannot be larger than number of columns in data set.")
  }

  ### extract relevant columns:
  rawData <- rawData[, c(conc_col, meas_col)]
  colnames(rawData) <- c("Concentration", "Measurement")

  ### check if relevant columns are numeric:
  if (!is.numeric(rawData[["Concentration"]])) {
    stop("Concentration column must be numeric. Issue may come from non-fitting decimal separator or na.strings.")
  }
  if (!is.numeric(rawData[["Measurement"]])) {
    stop("Measurement column must be numeric. Issue may come from non-fitting decimal separator or na.strings.")
  }

  return(rawData)
}




#' Data preprocessing: Helper function to select all rows from a specific concentration level
#'
#' @param x **numeric(1)** \cr concentration level to select
#' @param rawData **data.frame** \cr data set to be filtered (e.g. result of \code{\link{readData}})
#'
#' @returns data.frame
#'
#' @examples
filterConcentrationLevel <- function(x,
                                     rawData) {
  result <- rawData[rawData$Concentration == x, ]
  return(result)
}



#' Data preprocessing: Helper function to check for sufficient number of replicates for a specific concentration level
#'
#' @param x **numeric(1)** \cr concentration level to check
#' @param data **list of data.frames** \cr list of data.frames (each dataframe contains data for a specific concentration level)
#' @param minNumber **integer(1)** \cr minimal number of data points per concentration level
#'
#' @returns **locgical(1)** \cr TRUE if there are enough replicates or FALSE if not
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



#' Clean data (remove 0s and NAs, remove concentration levels with insufficient number of replicates)
#'
#' @param rawData **data.frame** \cr data set to be cleaned, result of readData.
#' @param min_replicates **integer(1)** \cr Minimal number of replicates/data points per concentration level.
#'                                          Concentration levels with too few data points will be removed.
#'
#' @returns list of data.frames, each data.frame contains data for a specific concentration level
#' @export
#'
#' @examples
cleanData <- function(rawData, min_replicates) {

  ### check input arguments
  checkmate::assert_int(min_replicates)

  # Removing rows that contain unwanted 0 values (problems with log-transform later) or NA values in either
  # the concentration or measurement column
  dataValidated <- rawData[rawData$Concentration != 0 & !is.na(rawData$Concentration) & rawData$Measurement != 0 & !is.na(rawData$Measurement),]

  # Determination of existing concentration levels in the validated data
  concLevels <- unique(dataValidated$Concentration)
  concLevels <- sort(concLevels, decreasing = FALSE)

  # Transforming a data set into a list with entries for each concentration level (and the related data)
  dataValidated <- lapply(concLevels, FUN = filterConcentrationLevel, rawData = dataValidated)

  # Deleting concentration levels with insufficient number of replicates
  dataValidated <- dataValidated[sapply(1:length(dataValidated), FUN = checkNumberReplicates, data = dataValidated, minNumber = min_replicates)]

  dataValConcLevels <- sapply(dataValidated, FUN = function(x) x$Concentration[1])
  names(dataValidated) <- dataValConcLevels

  return(dataValidated)
}

