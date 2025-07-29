test_that("prediction", {
    file <- system.file("extdata", "MSQC1/msqc1_dil_GGPFSDSYR.rds", package = "CalibraCurve")
    D <- readDataSE(file, concColName = "amount_fmol", substColName = "Substance")
    RES <- CalibraCurve(D, verbose = FALSE)
    newdata <- c(1000000, #1e6
                 10000000, #1e7
                 100000000) # 2e7

    pred <- predictConcentration(RES$RES[[4]], newdata = newdata, verbose = FALSE)

    expect_equal(nrow(pred), 3)
    expect_equal(ncol(pred), 3)
    expect_equal(colnames(pred), c("intensity", "predicted_concentrations", "linear_range"))
    expect_equal(pred$linear_range, c(TRUE, TRUE, FALSE))
})
