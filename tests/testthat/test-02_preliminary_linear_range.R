test_that("preliminary linear range", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx", package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)
    D_cleaned <- cleanData(D)
    RES_PLR <- calculate_PLR(D_cleaned, cvThres = 4,
        calcContinuousPrelimRanges = FALSE)

    expect_equal(names(RES_PLR), c("dataPrelim", "concLevelsCV", "prelimConcLevels"))
    expect_equal(length(RES_PLR$dataPrelim), 6)
    expect_equal(class(RES_PLR$dataPrelim), "list")
    expect_equal(length(RES_PLR$concLevelsCV), 6)
    expect_equal(RES_PLR$prelimConcLevels, c(25,50,200,1000,2000,5000))


    RES_PLR2 <- calculate_PLR(D_cleaned, cvThres = 4,
                             calcContinuousPrelimRanges = TRUE)
    expect_equal(RES_PLR2$prelimConcLevels, c(1000, 2000, 5000))

    expect_error(calculate_PLR(D_cleaned, cvThres = 1.4, calcContinuousPrelimRanges = FALSE))
})
