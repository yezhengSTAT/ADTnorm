test_that("marker_to_process works with zero-sum marker present", {
    d <- make_test_data_zerosum()
    result <- expect_no_error(suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        marker_to_process = c("CD3", "CD19"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    )))
    expect_equal(ncol(result), 2)
    expect_true("CD3" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})
