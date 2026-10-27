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
