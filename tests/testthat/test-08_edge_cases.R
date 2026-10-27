test_that("CalibraCurve works with unequal replicate counts per concentration", {
    D_unequal <- data.frame(
        Concentration = c(rep(1, 2), rep(2, 3), rep(5, 4)),
        Measurement = c(10.0, 10.1, 20.0, 19.9, 20.1, 50.0, 49.8, 50.1, 50.2)
    )

    expect_no_error({
        res <- CalibraCurve(D_unequal, minReplicates = 2, verbose = FALSE)
    })

    expect_equal(length(res$RES), 1)
    expect_true(is.numeric(res$RES[[1]]$final_linear_range))
    expect_true(length(res$RES[[1]]$final_linear_range) >= 2)
})


test_that("CalibraCurve drops levels below minReplicates and still computes", {
    D_mixed <- data.frame(
        Concentration = c(rep(1, 1), rep(2, 2), rep(5, 3), rep(10, 3)),
        Measurement = c(10.0, 20.0, 20.2, 50.0, 50.1, 49.9, 100.0, 100.2, 99.8)
    )

    cleaned <- cleanData(D_mixed, minReplicates = 2)

    expect_false(1 %in% as.numeric(names(cleaned)))
    expect_equal(sort(as.numeric(names(cleaned))), c(2, 5, 10))

    expect_no_error({
        res <- CalibraCurve(D_mixed, minReplicates = 2, verbose = FALSE)
    })

    expect_equal(length(res$RES), 1)
    expect_true(is.numeric(res$RES[[1]]$final_linear_range))
    expect_true(length(res$RES[[1]]$final_linear_range) >= 2)
})


test_that("CalibraCurve skips failed curves and omits NA-suffixed output files", {
    good <- data.frame(
        Concentration = c(1, 1, 1, 2, 2, 2, 5, 5, 5, 10, 10, 10),
        Measurement = c(10.0, 10.5, 9.8, 20.0, 20.2, 19.9, 50.0, 50.1, 50.2,
                        100.0, 99.8, 100.3)
    )
    bad <- data.frame(
        Concentration = rep(1, 9),
        Measurement = c(1.0, 1.2, 1.1, 1.3, 1.1, 1.2, 1.4, 1.1, 1.2)
    )

    out_dir <- tempfile("cc_skip_failed_")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    expect_no_error({
        res <- CalibraCurve(list(good = good, bad = bad), output_path = out_dir,
                           verbose = FALSE)
    })

    expect_equal(names(res$RES), "good")
    files <- list.files(out_dir)
    expect_false(any(grepl("_NA\\.", files)))
    expect_true(any(grepl("result_table_conc_levels_good\\.xlsx", files)))
    expect_false(file.exists(file.path(out_dir, "result_table_conc_levels_NA.xlsx")))
})
