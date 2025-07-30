test_that("readData xlsx and SummarizedExperiment", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx", package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
        measCol = 12)

    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D_list <- readDataSE(dataPath = file, concColName = "amount_fmol",
                         substColName = "Substance", assayNumber = 1)

    data_folder <- system.file("extdata", "MSQC1_xlsx", package = "CalibraCurve")
    D_list2 <- readMultipleTables(dataFolder = data_folder, fileType = "xlsx",
        concCol = 16, measCol = 12)

    rownames(D) <- NULL
    rownames(D_list[[4]]) <- NULL
    rownames(D_list2[[4]]) <- NULL

    expect_equal(D, D_list[[4]])
    expect_equal(D, D_list2[[4]])

    expect_equal(nrow(D), 18)
    expect_equal(ncol(D), 2)
    expect_equal(colnames(D), c("Concentration", "Measurement"))

    expect_equal(length(D_list), 9)
    expect_equal(unname(vapply(D_list, nrow, numeric(1))), rep(18,9))
    expect_equal(unname(vapply(D_list, ncol, numeric(1))), rep(2,9))

    expect_equal(length(D_list2), 9)
    expect_equal(unname(vapply(D_list2, nrow, numeric(1))), rep(18,9))
    expect_equal(unname(vapply(D_list2, ncol, numeric(1))), rep(2,9))
})


test_that("cleanData", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx", package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)

    D_cleaned <- cleanData(D)
    ### should result in a list of length 7, each containing a data.frame with 3 rows
    expect_equal(length(D_cleaned), 6)
    expect_equal(class(D_cleaned), "list")
    expect_equal(nrow(D_cleaned[[1]]), 3)

    ### error because no concentration level is left
    expect_error(cleanData(D_cleaned, minReplicates = 5))
})


test_that("cleanData with zeros and NAs", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx", package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)

    ### introduce zeros and NAs to the data
    D[1, 1] <- 0
    D[2, 2] <- NA

    D_cleaned <- cleanData(D, minReplicates = 2)

    expect_equal(length(D_cleaned), 5)
    expect_equal(class(D_cleaned), "list")
    expect_equal(nrow(D_cleaned[[1]]), 3)

    D_cleaned <- cleanData(D, minReplicates = 1)
    ### first concentration should be removed (because only 3 replicates are left)
    expect_equal(length(D_cleaned), 6)
    expect_equal(class(D_cleaned), "list")
    expect_equal(nrow(D_cleaned[[1]]), 1)
    expect_equal(nrow(D_cleaned[[2]]), 3)
})
