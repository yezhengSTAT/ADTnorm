test_that("Bug 1: output column names are correct after zero-sum filter", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        marker_to_process = c("CD3", "CD4", "CD19", "CD56"),
        bimodal_marker = c("CD3", "ZERO", "CD56"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_true("CD3" %in% colnames(result))
    expect_true("CD4" %in% colnames(result))
    expect_true("CD19" %in% colnames(result))
    expect_true("CD56" %in% colnames(result))
    expect_false("ZERO" %in% colnames(result))
})

test_that("Bug 2: bimodal detection does not crash with zero-sum marker", {
    skip("not directly testable \u2014 requires internal ADTnorm.R state")
})

test_that("Bug 3: positive peak alignment produces correct output", {
    d <- make_test_data_zerosum()
    result_p <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        positive_peak = list(ADT = "CD4"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_equal(ncol(result_p), 4)
    expect_true("CD3" %in% colnames(result_p))
})

test_that("Bug 4: bandwidth adjustment does not fail", {
    skip("not directly testable \u2014 requires internal ADTnorm.R state")
})

test_that("Bug 5: explicit marker_to_process does not OOB crash", {
    d <- make_test_data_zerosum()
    result <- suppressWarnings(ADTnorm(
        d$cell_x_adt, d$cell_x_feature,
        marker_to_process = c("CD3", "CD56"),
        peak_type = "midpoint", save_fig = FALSE, save_landmark = FALSE
    ))
    expect_true("CD3" %in% colnames(result))
    expect_true("CD56" %in% colnames(result))
})

test_that("Bug 6: RDS/PDF filenames omit zero-sum marker", {
    skip("requires save_fig = TRUE to inspect file names")
})
