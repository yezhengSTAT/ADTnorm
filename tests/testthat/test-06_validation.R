test_that("Fix 5: all-zero markers produce clear error", {
    d <- make_test_data_allzero()
    expect_error(suppressWarnings(ADTnorm(d$cell_x_adt, d$cell_x_feature,
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE)),
        "zero total counts")
})

test_that("Fix 4: bimodal warning emitted for zero-sum marker", {
    d <- make_test_data_zerosum()
    expect_message(suppressWarnings(ADTnorm(d$cell_x_adt, d$cell_x_feature,
        bimodal_marker = c("CD3", "ZERO", "CD56"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE)),
        "bimodal")
})

test_that("Fix 10: descriptive error for all-zero marker_to_process", {
    d <- make_test_data_zerosum()
    expect_error(suppressWarnings(ADTnorm(d$cell_x_adt, d$cell_x_feature,
        marker_to_process = "ZERO",
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE)),
        "removed|found")
})

test_that("Fix 9: all-zero bimodal marker doesn't crash", {
    d <- make_test_data_zerosum()
    result <- expect_no_error(suppressWarnings(ADTnorm(d$cell_x_adt,
        d$cell_x_feature, bimodal_marker = "ZERO",
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE)))
    expect_equal(ncol(result), 4)
    expect_false("ZERO" %in% colnames(result))
})

test_that("Fix 3: positive_peak with zero-sum marker stays synced", {
    d <- make_test_data_zerosum()
    result <- expect_no_error(suppressWarnings(ADTnorm(d$cell_x_adt,
        d$cell_x_feature,
        positive_peak = list(ADT = c("CD3", "ZERO"), sample = c("S1", "S2")),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE)))
    expect_equal(ncol(result), 4)
    expect_false("ZERO" %in% colnames(result))
})
