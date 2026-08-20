test_that("estimate_gene requires a single offset column name", {
  dat <- airway_one_gene_tbl()
  expect_error(estimate_gene(dat, ~ dex), "offset")
  expect_error(estimate_gene(dat, ~ dex, offset = NULL), "offset")
  expect_error(estimate_gene(dat, ~ dex, offset = ""), "offset")
  expect_error(estimate_gene(dat, ~ dex, offset = c("offset", "lib")), "offset")
  expect_error(
    estimate_gene(dat, ~ dex, offset = "missing_offset"),
    "not found"
  )
  expect_error(
    estimate_gene(dat, ~ dex, offset = "offset", dispersion = "dispersion"),
    "not found"
  )
})

test_that("prepare_gene_data rejects NA counts rather than dropping them", {
  dat <- airway_one_gene_tbl()
  dat$counts[2] <- NA_integer_
  expect_error(
    brmDE:::prepare_gene_data(dat, abundance = "counts", offset = "offset"),
    "1 NA"
  )
})

test_that("prepare_gene_data keeps offset and coerces counts", {
  dat <- airway_one_gene_tbl()
  prepared <- brmDE:::prepare_gene_data(dat, abundance = "counts", offset = "offset")
  expect_equal(nrow(prepared$data), nrow(dat))
  expect_equal(prepared$offset, "offset")
  expect_type(prepared$data$counts, "integer")
})

test_that("prepare_formula adds response and offset", {
  se <- airway_one_gene()
  expect_true(all(c("dex", "cell") %in% colnames(SummarizedExperiment::colData(se))))
  f <- brmDE:::prepare_formula(~ dex + (1 | cell), abundance = "counts", offset = "offset")
  txt <- paste(deparse(f), collapse = " ")
  expect_match(txt, "counts")
  expect_match(txt, "offset\\s*\\(\\s*offset\\s*\\)")
})

test_that("prepare_formula does not duplicate an existing offset", {
  se <- airway_one_gene()
  expect_true("offset" %in% colnames(SummarizedExperiment::colData(se)))
  f <- brmDE:::prepare_formula(
    counts ~ dex + offset(offset),
    abundance = "counts",
    offset = "offset"
  )
  txt <- paste(deparse(f), collapse = " ")
  expect_equal(length(gregexpr("offset\\s*\\(", txt)[[1]]), 1)
})

test_that("sanitize_names collapses repeated underscores", {
  dat <- airway_one_gene_tbl()
  dat$assay_groups___altered <- as.character(dat$dex)
  prepared <- brmDE:::prepare_gene_data(
    dat,
    abundance = "counts",
    offset = "offset",
    sanitize_names = TRUE
  )
  expect_true("assay_groups_altered" %in% names(prepared$data))
  expect_false("assay_groups___altered" %in% names(prepared$data))
})

test_that("random_intercept_vs_rest_from_names builds level-vs-rest contrasts", {
  se <- airway_one_gene()
  cells <- levels(se$cell)
  params <- paste0("r_cell[", cells, ",Intercept]")
  eq <- brmDE:::random_intercept_vs_rest_from_names(
    params,
    grouping = "cell"
  )
  expect_named(eq, cells)
  expect_match(
    eq[[cells[[1]]]],
    sprintf("`cell\\[%s,Intercept\\]`", cells[[1]])
  )
  expect_match(eq[[cells[[1]]]], cells[[2]])
  expect_no_match(
    eq[[cells[[1]]]],
    paste0(cells[[1]], ",Intercept].*", cells[[1]], ",Intercept")
  )
})

test_that("nullify_newdata zeroes offset and NAs selected covariates", {
  dat <- airway_one_gene_tbl()
  out <- brmDE:::nullify_newdata(
    dat,
    nullify = "dex",
    offset = "offset",
    offset_value = 0
  )
  expect_equal(out$offset, rep(0, nrow(dat)))
  expect_true(all(is.na(out$dex)))
  expect_equal(out$cell, dat$cell)
})

test_that("nullify_newdata rejects a column that is not in the model data", {
  dat <- airway_one_gene_tbl()
  expect_error(
    brmDE:::nullify_newdata(dat, nullify = c("dex", "missing_col")),
    "missing_col"
  )
})

test_that("as_gene_tibble accepts a one-gene SummarizedExperiment", {
  se <- airway_one_gene()
  tbl <- brmDE:::as_gene_tibble(se, abundance = "counts")
  expect_equal(nrow(tbl), ncol(se))
  expect_equal(unique(tbl$.feature), rownames(se))
  expect_equal(
    tbl$counts,
    as.integer(SummarizedExperiment::assay(se, "counts")[1, ])
  )
})

test_that("default ZINB priors omit shape submodel terms without shape ~", {
  dat <- airway_one_gene_tbl()
  priors <- c(
    brmDE:::zinb_location_priors(dat, "counts", "offset"),
    brms::prior(student_t(3, 0, 2), class = shape)
  )
  expect_s3_class(priors, "brmsprior")
  expect_true("Intercept" %in% priors$class)
  expect_true("shape" %in% priors$class)
  expect_false(any(priors$dpar == "shape"))
})

test_that("as_gene_tibble rejects multi-gene SummarizedExperiment", {
  se <- airway_se()[seq_len(2), ]
  expect_error(brmDE:::as_gene_tibble(se), "one gene")
})
