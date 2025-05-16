test_that("assemble result tables", {

  data(D_MFAP4)
  D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
  RES_PLR <- calculate_PLR(D_MFAP4_cleaned,
                           cv_thres = 20,
                           calcContinuousPrelimRanges = TRUE)
  RES_FLR <- calculate_FLR(RES_PLR$dataPrelim)

  resFacDataV <- CalibraCurve::calcRFLevels(D_MFAP4_cleaned, mod = RES_FLR$mod)
  avgResFacDataV <- CalibraCurve::calcRFMeans(resFacDataV)


  tables <- CalibraCurve::assemble_results(X = D_MFAP4,
                                           dataCleaned = D_MFAP4_cleaned,
                                           cv_thres = 20,
                                           PLR_res = RES_PLR,
                                           resFacDataV = resFacDataV,
                                           avgResFacDataV = avgResFacDataV,
                                           FLR_res = RES_FLR,
                                           mod = RES_FLR$mod,
                                           RfThresL = 80,
                                           RfThresU = 120,
                                           substance = "substance1"
  )

  expect_equal(length(tables), 2)
  expect_equal(class(tables), "list")
  expect_equal(names(tables), c("result_table_conc_levels", "result_table_obs"))

  expect_equal(nrow(tables$result_table_conc_levels), 11)
  expect_equal(ncol(tables$result_table_conc_levels), 14)
  expect_equal(colnames(tables$result_table_conc_levels),
               c("substance", "concentration", "mean_measurement", "estimated_measurement",
                 "removed_while_cleaning",
                 "CV", "CV_within_thres", "preliminary_linear_range", "mean_percentage_bias",
                 "SD_percentage_bias", "CV_percentage_bias", "mean_response_factor",
                 "RF_within_thres", "final_linear_range"))


  expect_equal(nrow(tables$result_table_obs), 11*5)
  expect_equal(ncol(tables$result_table_obs), 8)
  expect_equal(colnames(tables$result_table_obs),
               c("substance", "concentration", "measurement", "removed_while_cleaning",
                 "percentage_bias", "response_factor",
                 "RF_within_thres", "final_linear_range"))
})



