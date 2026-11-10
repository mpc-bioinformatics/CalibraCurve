test_that("package exposes a Shiny app launcher", {
  skip_if_not_installed("shiny")

  app_obj <- CalibraCurveApp()
  expect_s3_class(app_obj, "shiny.appobj")
})

test_that("import errors highlight wrong columns before file parsing issues", {
  err <- simpleError("Concentration and measurement columns must be numeric.")

  msg <- .format_import_error_message(err, file_type = "csv")

  expect_match(msg, "correct concentration and measurement columns are selected", ignore.case = TRUE)
  expect_match(msg, "wrong decimal separator|wrong.*separator", ignore.case = TRUE)
  expect_match(msg, "\\.tsv, \\.txt and \\.csv|\\.csv, \\.tsv and \\.txt", ignore.case = TRUE)
  expect_false(grepl("Original error", msg, ignore.case = TRUE))
})

test_that("multi-substance calibration results can be subset by substance", {
  skip_if_not_installed("shiny")

  file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
    package = "CalibraCurve"
  )
  D_single <- readDataSE(
    dataPath = file,
    concColName = "amount_fmol",
    substColName = "Substance",
    assayNumber = 1
  )
  D_list <- list(
    first = D_single[[1]],
    second = D_single[[1]]
  )
  names(D_list) <- c("substance_A", "substance_B")

  result <- CalibraCurve(
    D_list = D_list,
    plot_type = "single_plots",
    verbose = FALSE,
    output_path = NULL
  )

  subset_result <- .select_calibration_dataset(result, "substance_A")

  expect_equal(names(subset_result$RES), "substance_A")
  expect_equal(nrow(subset_result$summary_tab), 1L)
  expect_identical(subset_result$plot_CC_list[[1]], result$plot_CC_list[[1]])
})

test_that("calibration exports to a zip bundle with all plots and tables", {
  file <- system.file("extdata", "MSQC1", "msqc1_dil_GGPFSDSYR.rds",
    package = "CalibraCurve"
  )
  D_single <- readDataSE(
    dataPath = file,
    concColName = "amount_fmol",
    substColName = "Substance",
    assayNumber = 1
  )
  result <- CalibraCurve(
    D_list = list(substance_A = D_single[[1]], substance_B = D_single[[1]]),
    plot_type = "single_plots",
    verbose = FALSE,
    output_path = NULL
  )

  zip_path <- tempfile(fileext = ".zip")
  .export_calibration_results_zip(
    result,
    zip_path = zip_path,
    table_format = "csv",
    plot_format = "png",
    plot_width = 6,
    plot_height = 4,
    plot_dpi = 72
  )

  entries <- utils::unzip(zip_path, list = TRUE)$Name

  expect_true(file.exists(zip_path))
  expect_true(file.size(zip_path) > 0)
  expect_false(any(grepl("/", entries)))
  expect_true(any(grepl("calibration_curve_.*substance_A", entries)))
  expect_true(any(grepl("response_factor_plot_.*substance_B", entries)))
})
