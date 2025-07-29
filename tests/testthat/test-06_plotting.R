test_that("plot calibration curve", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
    package = "CalibraCurve")
    D_list <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                         assayNumber = 1)
    CC_RES <- calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = TRUE, substance = "QTRAP_y5")

    suppressWarnings({
        pl_CC <- CalibraCurve::plotCalibraCurve(list("QTRAP_y5" = CC_RES),
            show_regression_info = TRUE,
            show_linear_range = TRUE,
            show_data_points = TRUE
        )
        vdiffr::expect_doppelganger("CC_testfile_1", pl_CC$CC_plot)
    })
})




test_that("plot response factor plot", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
                        package = "CalibraCurve")
    D_list <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                         assayNumber = 1)
    CC_RES <- calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = FALSE)


    pl_CC <- CalibraCurve::plotResponseFactors(RES = CC_RES)
    vdiffr::expect_doppelganger("RF_testfile_1", pl_CC)

})
