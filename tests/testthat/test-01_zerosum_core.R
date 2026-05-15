test_that("zero-sum column is removed from output", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result), 4)
})

test_that("output column names exclude zero-sum marker", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_true("CD3" %in% colnames(result))
    expect_true("CD4" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_true("CD56" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})
