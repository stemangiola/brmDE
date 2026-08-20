test_that("add_tidybulk_offset stores log(1/multiplier) on the full SE", {
  skip_if_not_installed("tidybulk")

  se <- airway_for_hpc()
  out <- suppressWarnings(add_tidybulk_offset(se, abundance = "counts", method = "TMM"))
  expect_s4_class(out, "SummarizedExperiment")
  expect_equal(nrow(out), nrow(se))
  expect_true("multiplier" %in% colnames(SummarizedExperiment::colData(out)))
  expect_true("offset" %in% colnames(SummarizedExperiment::colData(out)))
  expect_equal(out$offset, log(1 / out$multiplier))
  expect_true("counts_scaled" %in% SummarizedExperiment::assayNames(out))
})

test_that("add_tidybulk_offset rejects non-SE input", {
  se <- airway_one_gene()
  expect_error(
    add_tidybulk_offset(as.data.frame(SummarizedExperiment::colData(se))),
    "SummarizedExperiment"
  )
})
