test_that("final linear range", {
    data(D_MFAP4)
    D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
    RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
        cv_thres = 10,
        calcContinuousPrelimRanges = TRUE
    )
    RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

    expect_equal(names(RES_FLR), c("dataFinal", "mod", "perBias", "perBiasAvgSDCV"))
    expect_equal(length(RES_FLR$dataFinal), 3)
    expect_equal(class(RES_FLR$dataFinal), "list")
    expect_equal(class(RES_FLR$mod), "lm")
    expect_equal(length(RES_FLR$perBias), 3)
    expect_equal(length(RES_FLR$perBias[[1]]), 5)
    expect_equal(nrow(RES_FLR$perBiasAvgSDCV), 3)
    expect_equal(ncol(RES_FLR$perBiasAvgSDCV), 3)
})
