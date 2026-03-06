test_that("final linear range", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx",
                        package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)
    D_cleaned <- cleanData(D)
    RES_PLR <- calculate_PLR(D_cleaned, cvThres = 10,
                             calcContinuousPrelimRanges = TRUE)
    RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

    expect_equal(names(RES_FLR),
                 c("dataFinal", "mod", "perBias", "perBiasAvgSDCV"))
    expect_equal(length(RES_FLR$dataFinal), 4)
    expect_equal(as.numeric(names(RES_FLR$dataFinal)), c(25,50,200,1000))
    expect_equal(class(RES_FLR$dataFinal), "list")
    expect_equal(class(RES_FLR$mod), "lm")
    expect_equal(class(RES_FLR$perBias), "list")
    expect_equal(length(RES_FLR$perBias), 4)
    expect_equal(unname(lengths(RES_FLR$perBias)), rep(3, 4))
    expect_equal(nrow(RES_FLR$perBiasAvgSDCV), 4)
    expect_equal(ncol(RES_FLR$perBiasAvgSDCV), 3)
    expect_equal(colnames(RES_FLR$perBiasAvgSDCV),
                 c("avgPerBias", "stdDevPerBias", "CV_PerBias"))
})
