test_that("response factors", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx",
                        package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)
    D_cleaned <- cleanData(D)
    RES_PLR <- calculate_PLR(D_cleaned, cvThres = 10,
                             calcContinuousPrelimRanges = TRUE)
    RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

    RFs <- CalibraCurve::calcRF(D_cleaned, mod = RES_FLR$mod)

    expect_equal(length(RFs$RFs), 6)
    expect_equal(class(RFs$RFs), "list")
    expect_equal(length(RFs$meanRFs), 6)
})
