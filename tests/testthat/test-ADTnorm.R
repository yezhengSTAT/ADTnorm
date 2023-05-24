test_that("ADTnorm works", {
  data(cell_x_adt)
  data(cell_x_feature)
  save_outpath <- tempdir()
  study_name <- "ADTnorm_demoRun"

  suppressWarnings({
    res <- ADTnorm(
      cell_x_adt = cell_x_adt,
      cell_x_feature = cell_x_feature,
      save_outpath = save_outpath,
      study_name = study_name,
      marker_to_process = c("CD3", "CD4", "CD8")
    )
  })

  expect_type(res, "list")
  expect_equal(nrow(res), 422682)
  expect_equal(ncol(res), 3)
})
