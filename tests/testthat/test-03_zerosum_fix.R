test_that("basic zero-sum filter removes zero-sum columns", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result), 4)
    expect_true("CD3" %in% colnames(result))
    expect_true("CD4" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_true("CD56" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})

test_that("explicit marker_to_process works with zero-sum marker", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        marker_to_process = c("CD3", "CD19"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result), 2)
    expect_true("CD3" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})

test_that("bimodal markers with ZERO in list excludes zero-sum marker", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        bimodal_marker = c("CD3", "ZERO", "CD56"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result), 4)
    expect_true("CD3" %in% colnames(result))
    expect_true("CD4" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_true("CD56" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})

test_that("positive peak alignment completes without error", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        positive_peak = list(ADT = c("CD4", "CD56")),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result), 4)
})

test_that("all-zero input errors with descriptive message", {
    d <- make_test_data_allzero()
    expect_error(
        suppressWarnings(ADTnorm(
            d$cell_x_adt, d$cell_x_feature,
            peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
        )),
        "zero total counts"
    )
})
