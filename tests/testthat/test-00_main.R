test_that("main_multiplot", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance")
    RES <- CalibraCurve(D, verbose = FALSE, plot_type = "multiplot")

    suppressWarnings({
        pl_CC <- RES$plot_CC_list
        vdiffr::expect_doppelganger("CC_testfile_2", pl_CC)
    })
})


test_that("main_allin1", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance")
    RES <- CalibraCurve(D, verbose = FALSE, plot_type = "all_in_one")

    suppressWarnings({
        pl_CC <- RES$plot_CC_list
        vdiffr::expect_doppelganger("CC_testfile_3", pl_CC)
    })
})


test_that("main_RF_multiplot", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance")
    RES <- CalibraCurve(
        D,
        verbose = FALSE,
        plot_type = "single_plots",
        RF_plot_type = "multiplot"
    )

    pl_RF <- RES$plot_RF_list[[1]]
    vdiffr::expect_doppelganger("RF_testfile_main_multiplot", pl_RF)
})


test_that("main_RF_allin1", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance")
    RES <- CalibraCurve(
        D,
        verbose = FALSE,
        plot_type = "single_plots",
        RF_plot_type = "all_in_one"
    )

    pl_RF <- RES$plot_RF_list[[1]]
    vdiffr::expect_doppelganger("RF_testfile_main_all_in_one", pl_RF)
})
