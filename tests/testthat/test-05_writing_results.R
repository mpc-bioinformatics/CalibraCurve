test_that("assemble result tables", {
    file <- system.file("extdata", "MSQC1_xlsx/GGPFSDSYR_QTRAP_y5.xlsx",
                        package = "CalibraCurve")
    D <- readDataTable(file, fileType = "xlsx", concCol = 16,
                       measCol = 12)
    D_cleaned <- cleanData(D)
    RES_PLR <- calculate_PLR(D_cleaned, cvThres = 10,
                             calcContinuousPrelimRanges = TRUE)
    RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

    RFs <- CalibraCurve::calcRF(D_cleaned, mod = RES_FLR$mod)


    tables <- CalibraCurve::assemble_results(
        X = D,
        dataCleaned = D_cleaned,
        cvThres = 20,
        PLR_res = RES_PLR,
        resFacDataV = RFs$RFs,
        avgResFacDataV = RFs$meanRFs,
        FLR_res = RES_FLR,
        mod = RES_FLR$mod,
        RfThresL = 80,
        RfThresU = 120,
        substance = "substance1"
    )

    expect_equal(length(tables), 2)
    expect_equal(class(tables), "list")
    expect_equal(names(tables), c("result_table_conc_levels", "result_table_obs"))

    expect_equal(nrow(tables$result_table_conc_levels), 6)
    expect_equal(ncol(tables$result_table_conc_levels), 14)
    expect_equal(
        colnames(tables$result_table_conc_levels),
        c(
            "substance", "concentration", "mean_measurement", "estimated_measurement",
            "removed_while_cleaning",
            "CV", "CV_within_thres", "preliminary_linear_range", "mean_percentage_bias",
            "SD_percentage_bias", "CV_percentage_bias", "mean_response_factor",
            "RF_within_thres", "final_linear_range"
        )
    )


    expect_equal(nrow(tables$result_table_obs), 6 * 3)
    expect_equal(ncol(tables$result_table_obs), 8)
    expect_equal(
        colnames(tables$result_table_obs),
        c(
            "substance", "concentration", "measurement", "removed_while_cleaning",
            "percentage_bias", "response_factor",
            "RF_within_thres", "final_linear_range"
        )
    )
})
