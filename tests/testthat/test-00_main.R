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
