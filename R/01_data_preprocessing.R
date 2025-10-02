#' Read data in different table input formats (xlsx, csv or txt).
#' Extracts the two relevant columns (concentration and measurement) and
#' orders the data by increasing concentration.
#'
#' @param dataPath **character(1)** \cr Path to the data file
#'      (.csv, .txt or .xlsx file).
#' @param fileType **character(1)** \cr Type of file: "csv", "txt" or "xlsx".
#' @param concCol **integer(1)** \cr Column number of the concentration values.
#' @param measCol **integer** \cr Column number of the concentration values.
#' @param sep **character(1)** \cr The field separator, default is ",".
#' @param dec **character(1)** \cr Decimal separator, default is ".".
#' @param header **logical(1)** \cr If TRUE, first line is counted as column
#'      names. The default is TRUE.
#' @param naStrings **character** \cr Vector of strings which are to be
#'      interpreted as NA. The default is c("NA", "NaN", "Filtered", "#NV").
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files,
#'      default is to use the first sheet).
#'
#' @returns Data.frame with two numeric columns: Concentration and Measurement
#' @export
#'
#' @examples
#' ### xlsx file:
#' file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx",
#'     package = "CalibraCurve"
#' )
#' D <- readDataTable(file, fileType = "xlsx", concCol = 16, measCol = 12)
#'
readDataTable <- function(
    dataPath, fileType, concCol, measCol, sep = ",", dec = ".", header = TRUE,
    naStrings = c("NA", "NaN", "Filtered", "#NV"), sheet = 1) {
    ### check input parameters
    checkmate::assert_file_exists(dataPath)
    checkmate::assert_character(dataPath, len = 1)
    checkmate::assert_choice(fileType, c("csv", "txt", "xlsx"))
    checkmate::assert_int(concCol, lower = 1)
    checkmate::assert_int(measCol)
    checkmate::assert_character(sep, len = 1)
    checkmate::assert_character(dec, len = 1)
    checkmate::assert_flag(header)
    checkmate::assert_character(naStrings)
    checkmate::assert_int(sheet)
    if (concCol == measCol) stop("Concentration and measurement columns cannot
                            be identical.")
    if (fileType == "csv" | fileType == "txt") {
        rawData <- utils::read.table(dataPath, sep = sep, header = header,
            dec = dec, na.strings = naStrings
        )
    }
    if (fileType == "xlsx") {
        rawData <- openxlsx::read.xlsx(dataPath,
            colNames = header,
            sheet = sheet, na.strings = naStrings
        )
    }
    ### check if column numbers are valid
    if (measCol > ncol(rawData)) {
        stop("Number of measurement column cannot be larger than number of
            columns in data set.")
    }
    if (concCol > ncol(rawData)) {
        stop("Number of concentration column cannot be larger than number of
            columns in data set.")
    }
    ### extract relevant columns:
    rawData <- data.frame(
        "Concentration" = rawData[, concCol],
        "Measurement" = rawData[, measCol]
    )
    ### check if relevant columns are numeric:
    if (!is.numeric(rawData[, 1]) | !is.numeric(rawData[, 2])) {
        stop("Concentration and measurement columns must be numeric.
            Issue may come from non-fitting decimal separator or na.strings.")
    }
    ### sort by concentration level (from lowest to highest)
    rawData <- rawData[order(rawData$Concentration), ]
    return(rawData)
}


#' Read folder of files in different table input formats (xlsx, csv or txt).
#'
#' @param dataFolder **character(1)** \cr Folder containing either xlsx, csv or
#'    txt files
#' @param fileType **character(1)** \cr Type of file: "csv", "txt" or "xlsx".
#' @param concCol **integer(1)** \cr Column number of the concentration values.
#' @param measCol **integer** \cr Column number of the concentration values.
#' @param ... additional parameters to \code{\link{readDataTable}}
#'
#' @returns List of data.frame, each with two numeric columns:
#'    Concentration and Measurement
#' @export
#'
#' @examples
#' data_folder <- system.file("extdata", "MSQC1_xlsx",
#'     package = "CalibraCurve")
#' D_list <- readMultipleTables(
#'     dataFolder = data_folder, fileType = "xlsx",
#'     concCol = 16, measCol = 12
#' )
readMultipleTables <- function(dataFolder, fileType, concCol, measCol, ...) {
    allFiles <- setdiff(
        list.files(path = dataFolder),
        list.dirs(path = dataFolder, recursive = FALSE, full.names = FALSE)
    )

    fileTable <- data.frame(file = allFiles) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            substanceName = strsplit(basename(file), "\\.")[[1]][1],
            fullPath = paste0(dataFolder, "/", file),
            fileExt = tools::file_ext(file)
        )

    ### filter files for correct filetype
    fileTable <- fileTable[fileTable$fileExt == fileType, ]

    rawDataList <- lapply(fileTable$fullPath,
        FUN = readDataTable,
        fileType = fileType, concCol = concCol, measCol = measCol, ...
    )

    names(rawDataList) <- fileTable$substanceName
    return(rawDataList)
}









#' Read data stored as an SummarizedExperiment object (directly or stored in an
#' .rds file).
#' Extracts the two relevant columns (concentration and measurement) and
#' orders the data by increasing concentration level.
#'
#' @details
#' The SummarizedExperiments object may contain quantitative values from
#' targeted proteomics, lipidomics or metabolomics experiments.
#' The colData has to contain a column with the concentration levels
#' (concColName).
#' The rowData has to contain a column with the substance names (e.g. peptide
#' sequence, name of lipid or metabolite etc).
#'
#' @param dataPath **character(1)** \cr Path to the data file (.rds file)
#' @param rawDataSE **SummarizedExperiment** \cr SummarizedExperiment object
#' @param concColName **character(1)** \cr Name of the column in the colData()
#'    containing the concentration levels.
#' @param substColName **character(1)** \cr column name of rowData() containing
#' the substance name (must be a unique value in each row)
#' @param assayNumber **integer(1)** \cr Number of assay to be extracted
#'    from the SummarizedExperiment object
#' @param rowNumbers **integer** \cr Row numbers to extract from the
#'    SummarizedExperiment object. Default is NULL, which means that all rows
#'    in the object will be used.
#'
#' @returns List of data.frame, each with two numeric columns:
#'    Concentration and Measurement
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#'     package = "CalibraCurve"
#' )
#' D_list <- readDataSE(file,
#'     concColName = "amount_fmol",
#'     substColName = "Substance", assayNumber = 1
#' )
#'
#' # Alternative: import SummarizedExperiment object directly
#' rawDataSE <- readRDS(file)
#'
#' D_list2 <- readDataSE(rawDataSE = rawDataSE,
#'     concColName = "amount_fmol",
#'     substColName = "Substance", assayNumber = 1
#' )
readDataSE <- function(dataPath = NULL, rawDataSE = NULL, concColName,
                       substColName, assayNumber = 1, rowNumbers = NULL) {
    if (is.null(rawDataSE)) {
        checkmate::assert_file_exists(dataPath)
        rawDataSE <- readRDS(dataPath)
    }
    checkmate::assert_class(rawDataSE, "SummarizedExperiment")
    if (!is.null(rowNumbers)) rawDataSE <- rawDataSE[rowNumbers, ]

    Data <- SummarizedExperiment::assays(rawDataSE)[[assayNumber]]
    Data$Substance <- SummarizedExperiment::rowData(rawDataSE)[[substColName]]
    concentrations <- SummarizedExperiment::colData(rawDataSE)[, concColName]
    concentrations <- as.numeric(concentrations)
    colNames <- colnames(Data)

    rawData <- tidyr::pivot_longer(Data,
        cols = !tidyr::last_col(),
        names_to = "Concentration",
        values_to = "Measurement"
    )
    rawData$Concentration <- concentrations[match(
        rawData$Concentration,
        colNames
    )]
    rawData <- as.data.frame(rawData)

    rawData <- rawData[order(rawData$Concentration), ]
    rawDataList <- split(rawData, rawData$Substance)
    rawDataList <- lapply(rawDataList, function(x) x[, -1]) # rm substance col

    return(rawDataList)
}




#' Clean data (remove 0s and NAs, remove concentration levels with insufficient
#'  number of replicates)
#'
#' @param rawData **data.frame** \cr data set to be cleaned, result of
#'  \code{\link{readDataTable}} or \code{\link{readDataSE}}.
#' @param minReplicates **integer(1)** \cr Minimal number of replicates
#'  per concentration level. Concentration levels with too few data points will
#'  be removed.
#'
#' @returns list of data.frames, each element contains data for a specific
#'  concentration level
#' @export
#'
#' @examples
#' file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
#'     package = "CalibraCurve"
#' )
#' D_list <- readDataSE(file,
#'     concColName = "amount_fmol",
#'     substColName = "Substance", assayNumber = 1
#' )
#' cleanData(D_list[[1]])
cleanData <- function(rawData, minReplicates = 3) {
    ### check input arguments
    checkmate::assert_int(minReplicates, lower = 1)

    # Removing rows that contain unwanted 0 values (problems with log-transform
    # later) or NA values in either the concentration or measurement column
    dataCleaned <- rawData[rawData$Concentration != 0 &
        !is.na(rawData$Concentration) &
        rawData$Measurement != 0 &
        !is.na(rawData$Measurement), ]

    # Determination of existing concentration levels in the validated data
    concLevels <- unique(dataCleaned$Concentration)

    # Transforming a data set into a list with entries for each concentration
    # level (and the related data)
    dataCleaned <- lapply(concLevels,
        FUN = .filterConcentrationLevel,
        data = dataCleaned
    )

    # Deleting concentration levels with insufficient number of replicates
    ind <- vapply(seq_along(dataCleaned),
        FUN = .checkNumberReplicates,
        logical(1), data = dataCleaned, minReplicates = minReplicates
    )
    dataCleaned <- dataCleaned[ind]

    if (length(dataCleaned) <= 1) {
        stop(
            "One or less concentration level(s) with at least ", minReplicates,
            " replicates found. Please check your data or lower minReplicates."
        )
    }

    concLevelsCleaned <- vapply(dataCleaned,
        FUN = function(x) x$Concentration[1], numeric(1)
    )
    names(dataCleaned) <- concLevelsCleaned

    return(dataCleaned)
}





#' Data preprocessing: Helper function to select all rows from a specific
#' concentration level
#'
#' @param x **numeric(1)** \cr concentration level to select
#' @param data **data.frame** \cr data set to be filtered
#'    (e.g., result of \code{\link{readDataTable}})
#'
#' @returns data.frame
.filterConcentrationLevel <- function(x, data) {
    result <- data[data$Concentration == x, ]
    return(result)
}



#' Data preprocessing: Helper function to check for sufficient number of
#' replicates for a specific concentration level
#'
#' @param x **numeric(1)** \cr concentration level to check
#' @param data **list of data.frames** \cr list of data.frames (each dataframe
#'    contains data for a specific concentration level)
#' @param minReplicates **integer(1)** \cr minimal number of data points per
#'    concentration level
#'
#' @returns **logical(1)** \cr TRUE if there are enough replicates, else FALSE
.checkNumberReplicates <- function(x, data, minReplicates) {
    if (nrow(data[[x]]) < minReplicates) {
        result <- FALSE
    } else {
        result <- TRUE
    }
    return(result)
}
