

test_that("prediction", {

  data(RES_MFAP4)

  newdata <- c(5, 0.2, 10)
  pred <- predictConcentration(CC_res = list(RES = list("MFAP4" = RES_MFAP4)), newdata = newdata)

  expect_equal(nrow(pred), 3)
  expect_equal(ncol(pred), 3)
  expect_equal(colnames(pred), c("intensity", "predicted_concentrations", "linear_range"))


  ## concentrations that lead to predictions outside of the linear range
  newdata2 <- c(5, 0.02, 100)
  expect_warning(predictConcentration(CC_res = list(RES = list("MFAP4" = RES_MFAP4)), newdata = newdata2))


})
