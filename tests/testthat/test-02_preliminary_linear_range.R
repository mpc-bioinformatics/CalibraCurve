test_that("preliminary linear range", {
    data(D_MFAP4)
    D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
    RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
        cv_thres = 10,
        calcContinuousPrelimRanges = TRUE
    )
    ## With calcContinuousPrelimRanges = TRUE, only concentration levels 0.025-2.5
    ## are part of the preliminary linear range.

    expect_equal(names(RES_PLR), c("dataPrelim", "concLevelsCV", "prelimConcLevels"))
    expect_equal(length(RES_PLR$dataPrelim), 5)
    expect_equal(class(RES_PLR$dataPrelim), "list")
    expect_equal(length(RES_PLR$concLevelsCV), 11)
    expect_equal(RES_PLR$prelimConcLevels, c(0.025, 0.05, 0.25, 0.5, 2.5))



    RES_PLR2 <- calculate_PLR(D_MFAP4_cleaned,
        cv_thres = 10,
        calcContinuousPrelimRanges = FALSE
    )
    ## With calcContinuousPrelimRanges = FALSE, gaps are allowed, so the preliminary
    ## linear range is 0.025-25 (although 5 has a CV larger than the threshold 10).
    expect_equal(RES_PLR2$prelimConcLevels, c(0.025, 0.05, 0.25, 0.5, 2.5, 5, 25))

    expect_error(calculate_PLR(D_MFAP4_cleaned, cv_thres = 5, calcContinuousPrelimRanges = FALSE))
})
