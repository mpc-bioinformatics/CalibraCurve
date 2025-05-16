

test_that("prediction", {

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

  RES <- list(mod = RES_FLR$mod,
              final_linear_range = as.numeric(names(RES_FLR$dataFinal)),
              dataCleaned = D_MFAP4_cleaned,
              weightingMethod = "1/x^2",
              result_table_conc_levels = tables$result_table_conc_levels,
              result_table_obs = tables$result_table_obs)


  RES$result_table_conc_levels



  newdata <- c(5, 0.2, 10)
  pred <- predictConcentration(CC_res = RES, newdata = newdata)

  expect_equal(nrow(pred), 3)
  expect_equal(ncol(pred), 3)
  expect_equal(colnames(pred), c("intensity", "predicted_concentrations", "linear_range"))


  ## concentrations that lead to predictions outside of the linear range
  newdata2 <- c(5, 0.02, 100)
  expect_warning(predictConcentration(CC_res = RES, newdata = newdata2))


})
