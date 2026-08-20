test_that("estimate, hypothesis, and adjust append to one HPCell script", {
  skip_if_not_installed("HPCell")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()

  store <- tempfile("brmde_pipe_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
      unlink(paste0(store, "_input"), recursive = TRUE)
    },
    add = TRUE
  )

  pipeline <- se |>
    brmDE(store = store, features = c("ENSG00000120129")) |>
    estimate(~ dex, offset = "offset", dispersion = "dispersion") |>
    hypothesis("dextrt = 0") |>
    adjust(nullify = "dex")

  expect_s3_class(pipeline, "brmDE_hpc")
  expect_s3_class(pipeline, "HPCell")
  expect_true("se_dispersion" %in% names(pipeline))
  expect_true("brms_fit" %in% names(pipeline))
  expect_true("hypothesis_tbl" %in% names(pipeline))
  expect_true("adjust_tbl" %in% names(pipeline))

  lines <- readLines(paste0(store, ".R"))
  script <- paste(lines, collapse = "\n")
  expect_match(script, "hpc_internal")
  expect_match(script, "features\\.rds")
  expect_match(script, "estimate_dispersion_from_args")
  expect_match(script, "estimate_gene_from_se")
  expect_match(script, "hypothesis_gene_from_fit")
  expect_match(script, "adjust_gene_from_fit")
  expect_false(identical(trimws(lines[length(lines)]), "target_list"))
})

test_that("brmDE features must be a character vector", {
  se <- airway_one_gene()
  expect_error(brmDE(se, features = quote(ENSG00000120129)), "character vector")
  expect_error(brmDE(se, features = 1L), "character vector")
})

test_that("estimate requires offset and dispersion column names", {
  skip_if_not_installed("HPCell")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  store <- tempfile("brmde_offset_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
      unlink(paste0(store, "_input"), recursive = TRUE)
    },
    add = TRUE
  )
  pipeline <- brmDE(se, store = store, features = c("ENSG00000120129"))
  expect_error(estimate(pipeline, ~ dex), "offset")
  expect_error(estimate(pipeline, ~ dex, offset = "offset"), "dispersion")
})

test_that("estimate |> hypothesis |> adjust evaluate as one pipeline", {
  skip_on_cran()
  skip_if_not_installed("HPCell")
  skip_if_no_cmdstan()
  skip_if_not_installed("targets")
  skip_if_not_installed("tidybulk")

  se <- airway_for_hpc()

  store <- tempfile("brmde_pipe_eval_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
      unlink(paste0(store, "_input"), recursive = TRUE)
    },
    add = TRUE
  )

  out <- suppressWarnings(
    se |>
      brmDE(
        features = c("ENSG00000120129"),
        store = store,
        method = "TMM"
      ) |>
      estimate(
        ~ dex + (1 | cell),
        offset = "offset",
        dispersion = "dispersion",
        family = brms::negbinomial(),
        chains = 1,
        iter = 250,
        warmup = 100,
        cores = 1,
        backend = "cmdstanr",
        refresh = 0,
        silent = 2
      ) |>
      hypothesis("dextrt = 0") |>
      adjust(nullify = "dex") |>
      evaluate_hpc()
  )

  expect_equal(out$.feature, "ENSG00000120129")
  expect_s3_class(out$brms_fit[[1]], "brmsfit")
  expect_s3_class(out$hypothesis[[1]], "tbl_df")
  expect_s3_class(out$adjust[[1]], "tbl_df")
  expect_true("adjusted___Estimate" %in% names(out$adjust[[1]]) ||
    any(grepl("^adjusted___", names(out$adjust[[1]]))))
})
