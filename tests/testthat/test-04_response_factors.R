test_that("response factors", {
    data(D_MFAP4)
    D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
    RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
        cv_thres = 20,
        calcContinuousPrelimRanges = TRUE
    )
    RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

    resFacDataV <- CalibraCurve::calcRFLevels(D_MFAP4_cleaned, mod = RES_FLR$mod)

    expect_equal(length(resFacDataV), 11)
    expect_equal(class(resFacDataV), "list")

    # Calculation of mean response factor values
    avgResFacDataV <- CalibraCurve::calcRFMeans(resFacDataV)

    expect_equal(length(avgResFacDataV), 11)
})
