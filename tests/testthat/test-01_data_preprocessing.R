test_that("readData xlsx", {
  file <- system.file("extdata", "MFAP4_WTVFQK_y4.xlsx", package = "CalibraCurve")
  D <- readData(file,
                filetype = "xlsx",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_MFAP4)

  file <- system.file("extdata", "Apolipoprotein_TPAYYPNAGLIK_y7.xlsx", package = "CalibraCurve")
  D <- readData(file,
                filetype = "xlsx",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_Apolipoprotein)

  file <- system.file("extdata", "ALB_LVNEVTEFAK_y8.xlsx", package = "CalibraCurve")
  D <- readData(file,
                filetype = "xlsx",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_ALB)
})


test_that("readData csv", {
  file <- system.file("extdata", "MFAP4_WTVFQK_y4.csv", package = "CalibraCurve")
  D <- readData(file,
                filetype = "csv",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_MFAP4)

  file <- system.file("extdata", "Apolipoprotein_TPAYYPNAGLIK_y7.csv", package = "CalibraCurve")
  D <- readData(file,
                filetype = "csv",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_Apolipoprotein)

  file <- system.file("extdata", "ALB_LVNEVTEFAK_y8.csv", package = "CalibraCurve")
  D <- readData(file,
                filetype = "csv",
                conc_col = 6,
                meas_col = 7)
  expect_equal(D, CalibraCurve::D_ALB)
})




test_that("cleanData", {
  data(D_ALB)

  D_ALB_cleaned <- cleanData(D_ALB)
  ### should result in a list of length 7, each containing a data.frame with 3 rows
  expect_equal(length(D_ALB_cleaned), 7)
  expect_equal(class(D_ALB_cleaned), "list")
  expect_equal(nrow(D_ALB_cleaned[[1]]), 3)

  ### error because no concentration level is left
  expect_error(cleanData(D_ALB, min_replicates = 5))

})


test_that("cleanData with zeros and NAs", {
  data(D_MFAP4)

  ### introduce zeros and NAs to the data
  D_MFAP4[1,1] <- 0
  D_MFAP4[2,2] <- NA

  D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 3)
  ### should result in a list of length 11, each containing a data.frame with 5 rows (except the first one with 3 rows)
  expect_equal(length(D_MFAP4_cleaned), 11)
  expect_equal(class(D_MFAP4_cleaned), "list")
  expect_equal(nrow(D_MFAP4_cleaned[[1]]), 3)
  expect_equal(nrow(D_MFAP4_cleaned[[2]]), 5)

  D_MFAP4_cleaned <- cleanData(D_MFAP4, min_replicates = 5)
  ### first concentration should be removed (because only 3 replicates are left)
  expect_equal(length(D_MFAP4_cleaned), 10)
  expect_equal(class(D_MFAP4_cleaned), "list")
  expect_equal(nrow(D_MFAP4_cleaned[[1]]), 5)
  expect_equal(nrow(D_MFAP4_cleaned[[2]]), 5)

})
