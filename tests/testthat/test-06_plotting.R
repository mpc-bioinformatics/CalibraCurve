test_that("plot calibration curve", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
    package = "CalibraCurve")
    D_list <- CalibraCurve::readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                                       assayNumber = 1)
    CC_RES <- CalibraCurve::calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = TRUE, substance = "QTRAP_y5")

    suppressWarnings({
        pl_CC <- CalibraCurve::plotCalibraCurve(list("QTRAP_y5" = CC_RES),
            show_regression_info = TRUE,
            show_linear_range = TRUE,
            show_data_points = TRUE
        )
        vdiffr::expect_doppelganger("CC_testfile_1", pl_CC$CC_plot)
    })
})


test_that("plot calibration curve with linear x-axis", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
    package = "CalibraCurve")
    D_list <- CalibraCurve::readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                                       assayNumber = 1)
    CC_RES <- CalibraCurve::calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = TRUE, substance = "QTRAP_y5")

    suppressWarnings({
        pl_CC <- CalibraCurve::plotCalibraCurve(
            list("QTRAP_y5" = CC_RES),
            log10_x_axis = FALSE
        )
        vdiffr::expect_doppelganger("CC_testfile_1_linear_x_axis", pl_CC$CC_plot)
    })

    x_scale <- pl_CC$CC_plot$scales$get_scales("x")
    expect_equal(x_scale$trans$name, "identity")
})




test_that("plot response factor plot", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
                        package = "CalibraCurve")
    D_list <- CalibraCurve::readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                                       assayNumber = 1)
    CC_RES <- CalibraCurve::calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = FALSE)


    pl_CC <- CalibraCurve::plotResponseFactors(RES = CC_RES)
    vdiffr::expect_doppelganger("RF_testfile_1", pl_CC)

})


test_that("plot response factor plot with linear x-axis", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
                        package = "CalibraCurve")
    D_list <- CalibraCurve::readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                                       assayNumber = 1)
    CC_RES <- CalibraCurve::calc_single_curve(D_list[[4]], calcContinuousPrelimRanges = FALSE)

    pl_CC <- CalibraCurve::plotResponseFactors(
        RES = CC_RES,
        log10_x_axis = FALSE
    )

    vdiffr::expect_doppelganger("RF_testfile_1_linear_x_axis", pl_CC)

    x_scale <- pl_CC$scales$get_scales("x")
    expect_equal(x_scale$trans$name, "identity")
})


test_that("CalibraCurve propagates linear x-axis to RF plot", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds",
                        package = "CalibraCurve")
    D_list <- CalibraCurve::readDataSE(file, concColName = "amount_fmol", substColName = "Substance",
                                       assayNumber = 1)

    res <- suppressWarnings(CalibraCurve::CalibraCurve(
        D_list = D_list[4],
        plot_type = "single_plots",
        log10_x_axis = FALSE,
        verbose = FALSE
    ))

    x_scale <- res$plot_RF_list[[1]]$scales$get_scales("x")
    expect_equal(x_scale$trans$name, "identity")
})
