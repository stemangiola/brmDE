test_that("estimate, hypothesis, and adjust append to one tidytargets script", {
  skip_if_not_installed("tidytargets")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()

  store <- tempfile("brmde_pipe_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
    },
    add = TRUE
  )

  pipeline <- se |>
    brmDE(store = store, features = c("ENSG00000120129")) |>
    estimate(~ dex, offset = "offset",
             dispersion_prior_log_mean = "dispersion_trended") |>
    hypothesis("dextrt = 0") |>
    adjust(nullify = "dex")

  expect_s3_class(pipeline, "brmDE_hpc")
  expect_s3_class(pipeline, "tidytargets")
  expect_equal(pipeline$gene_id$iterate, "map")
  expect_true("brms_fit" %in% names(pipeline))
  expect_true("formula_abundance" %in% names(pipeline))
  expect_true("formula_dispersion" %in% names(pipeline))
  expect_equal(pipeline$formula_abundance$iterate, "none")
  expect_true("hypothesis_tbl" %in% names(pipeline))
  expect_true("adjust_tbl" %in% names(pipeline))

  # The pipeline remembers what it asked for, and remembers it through the
  # later steps. Nothing was said about sampling here, so those entries are
  # estimate_gene()'s own defaults rather than gaps.
  settings <- tidytargets::tt_metadata(pipeline)
  expect_equal(settings$abundance, "counts")
  expect_equal(settings$features, "ENSG00000120129")
  expect_equal(settings$gene_ids, list("ENSG00000120129"))
  expect_equal(settings$offset, "offset")
  expect_equal(settings$dispersion_prior_log_mean, "dispersion_trended")
  expect_equal(settings$chains, 2)
  expect_equal(settings$draws_warmup, 300)
  expect_equal(settings$draws_sampling, 500)

  lines <- readLines(paste0(store, ".R"))
  script <- paste(lines, collapse = "\n")
  expect_match(script, "tt_factory")
  # Formulas and extra arguments are tt_data() snapshots; the iterate
  # command names those targets so changing them invalidates the fits.
  expect_match(script, "estimate_gene_from_se")
  expect_match(script, "hypothesis_gene_from_fit")
  expect_match(script, "adjust_gene_from_fit")
  expect_true(file.exists(file.path(store, "se_input_data.qs")))
  expect_true(file.exists(file.path(store, "gene_id_data.qs")))
  expect_true(file.exists(file.path(store, "formula_abundance_data.qs")))
  expect_true(file.exists(file.path(store, "estimate_args_data.qs")))
  expect_true(file.exists(file.path(store, "temp_computing_resources.qs")))
  expect_true("gene_id" %in% names(pipeline))
  # Snapshots run on the main process, not the crew controller.
  expect_match(script, 'deployment = "main"')
  # Dispersion is estimated upstream of the pipeline, not inside it.
  expect_no_match(script, "estimate_dispersion", fixed = TRUE)
  expect_false(identical(trimws(lines[length(lines)]), "target_list"))
})

test_that("changed arguments change the targets script, unchanged ones do not", {
  skip_if_not_installed("tidytargets")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  stores <- character(0)
  on.exit(
    for (s in stores) {
      unlink(s, recursive = TRUE)
      unlink(paste0(s, ".R"))
    },
    add = TRUE
  )

  other_feature <- setdiff(rownames(se), "ENSG00000120129")[[1]]

  # Session values are snapshotted by tt_data() / tt_data_list(), so a
  # change shows up in the qs file rather than in the generated script.
  snapshot_for <- function(formula_abundance = ~dex,
                            formula_dispersion = ~1,
                            features = "ENSG00000120129",
                            family = brms::negbinomial()) {
    store <- tempfile("brmde_hash_")
    stores <<- c(stores, store)
    pipeline <- se |>
      brmDE(store = store, features = features) |>
      estimate(
        formula_abundance,
        formula_dispersion = formula_dispersion,
        offset = "offset",
        dispersion_prior_log_mean = "dispersion_trended",
        family = family
      )
    qs_hash <- function(name) {
      tools::md5sum(file.path(store, paste0(name, "_data.qs")))
    }
    list(
      script = gsub(
        basename(store), "STORE", readLines(paste0(store, ".R")),
        fixed = TRUE
      ),
      gene_ids = tidytargets::tt_metadata(pipeline)$gene_ids,
      formula_abundance = qs_hash("formula_abundance"),
      formula_dispersion = qs_hash("formula_dispersion"),
      estimate_args = qs_hash("estimate_args")
    )
  }

  base <- snapshot_for()
  expect_identical(snapshot_for()$script, base$script)
  expect_false(identical(
    snapshot_for(formula_abundance = ~ dex + (1 | cell))$formula_abundance,
    base$formula_abundance
  ))
  expect_false(identical(
    snapshot_for(formula_dispersion = ~cell)$formula_dispersion,
    base$formula_dispersion
  ))
  expect_false(identical(snapshot_for(features = other_feature)$gene_ids, base$gene_ids))
  expect_false(
    identical(
      snapshot_for(family = brms::zero_inflated_negbinomial())$estimate_args,
      base$estimate_args
    )
  )
})

test_that("brmDE features must be a character vector", {
  se <- airway_one_gene()
  expect_error(brmDE(se, features = quote(ENSG00000120129)), "character vector")
  expect_error(brmDE(se, features = 1L), "character vector")
})

test_that("estimate requires an offset column name", {
  skip_if_not_installed("tidytargets")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  store <- tempfile("brmde_offset_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
    },
    add = TRUE
  )
  pipeline <- brmDE(se, store = store, features = c("ENSG00000120129"))
  expect_error(estimate(pipeline, ~ dex), "offset")
  expect_s3_class(
    estimate(pipeline, ~ dex, offset = "offset"),
    "brmDE_hpc"
  )
})

test_that("bundle regroups the genes into fewer fit targets", {
  skip_if_not_installed("tidytargets")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  stores <- character(0)
  on.exit(
    for (s in stores) {
      unlink(s, recursive = TRUE)
      unlink(paste0(s, ".R"))
    },
    add = TRUE
  )

  pipeline_for <- function(...) {
    store <- tempfile("brmde_bundle_")
    stores <<- c(stores, store)
    se |>
      brmDE(store = store, features = rownames(se)[1:6]) |>
      estimate(~dex, offset = "offset", ...) |>
      hypothesis("dextrt = 0")
  }
  script_of <- function(pipeline) {
    store <- pipeline$initialisation$store
    gsub(basename(store), "STORE", readLines(paste0(store, ".R")), fixed = TRUE)
  }

  unbundled <- pipeline_for(bundle = 1)
  bundled <- pipeline_for(bundle = 3)

  # bundle = 1 is one gene per branch; the graph is the same, only n_units
  # changes.
  expect_equal(unbundled$gene_bundle$n_units, 6L)
  expect_equal(bundled$gene_bundle$n_units, 2L)
  expect_true(any(grepl('other_arguments_to_map = "gene_bundle"', script_of(unbundled), fixed = TRUE)))
  expect_true(any(grepl('other_arguments_to_map = "gene_bundle"', script_of(bundled), fixed = TRUE)))

  # hypothesis() maps over the fit target, so it inherits the coarser branching.
  expect_true(any(grepl('other_arguments_to_map = "brms_fit"', script_of(bundled), fixed = TRUE)))
})

test_that("bundle_gene_ids groups genes by size and keeps their order", {
  ids <- as.list(paste0("g", 1:10))

  expect_identical(unlist(bundle_gene_ids(ids, 4)), unlist(ids))
  # Bundles are `bundle` genes each, so only the last one is a remainder.
  expect_identical(lengths(bundle_gene_ids(ids, 4)), c(4L, 4L, 2L))
  expect_identical(lengths(bundle_gene_ids(ids, 5)), c(5L, 5L))
  # bundle = 1 is one gene per target.
  expect_identical(lengths(bundle_gene_ids(ids, 1)), rep(1L, 10))
  # A bundle bigger than the gene set is one target, not an empty one.
  expect_identical(lengths(bundle_gene_ids(ids, 50)), 10L)
})

test_that("bundle must be a positive whole number", {
  skip_if_not_installed("tidytargets")
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  store <- tempfile("brmde_nbundles_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
    },
    add = TRUE
  )
  pipeline <- brmDE(se, store = store, features = c("ENSG00000120129"))

  expect_error(estimate(pipeline, ~dex, offset = "offset", bundle = 0), "bundle")
  expect_error(estimate(pipeline, ~dex, offset = "offset", bundle = 2.5), "bundle")
  expect_error(estimate(pipeline, ~dex, offset = "offset", bundle = c(1, 2)), "bundle")
  expect_error(estimate(pipeline, ~dex, offset = "offset", bundle = NULL), "bundle")
})

test_that("bundled genes give one row per gene, as unbundled ones do", {
  skip_on_cran()
  skip_if_not_installed("tidytargets")
  skip_if_no_cmdstan()
  skip_if_not_installed("targets")

  se <- airway_for_hpc()
  features <- c(
    "ENSG00000120129",
    setdiff(rownames(se), "ENSG00000120129")[1:3]
  )

  store <- tempfile("brmde_bundle_eval_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
    },
    add = TRUE
  )

  pipeline <- se |>
    brmDE(features = features, store = store) |>
    estimate(
      ~dex,
      offset = "offset",
      dispersion_prior_log_mean = "dispersion_trended",
      bundle = 2,
      family = brms::negbinomial(),
      chains = 1,
      draws_warmup = 100,
      draws_sampling = 100,
      cores = 1,
      backend = "cmdstanr",
      refresh = 0,
      silent = 2
    ) |>
    hypothesis("dextrt = 0")

  # Nothing has been fitted yet: plotting evaluates the pipeline itself, as
  # printing it would, and floors its band at 1 / 100 from the metadata alone.
  band <- suppressWarnings(
    ggplot2::ggplot_build(plot_volcano(pipeline))$data[[1]]
  )
  expect_equal(sort(10^(-c(band$ymin, band$ymax))), c(1e-3, 1e-2))

  # A second evaluation skips everything it already has.
  out <- suppressWarnings(pipeline |> tt_evaluate())

  # Four genes, but only two fit targets were built.
  branches <- targets::tar_read_raw("brms_fit", store = store)
  expect_length(branches, 2L)
  fits <- unname(unlist(branches, recursive = FALSE))
  expect_true(all(vapply(fits, inherits, logical(1), "brmsfit")))

  # One hypothesis per gene here, so unnesting still leaves one row per gene.
  expect_equal(out$.feature, features)

  # estimate() contributes the gene column and nothing else: the fits are far
  # too large to gather, and the statistics that make the run readable are the
  # hypothesis ones, inline rather than in a list column.
  expect_false("hypothesis_tbl" %in% names(out))
  expect_true(
    all(c("hypothesis", "pH0", "rhat", "ess_bulk", "mcse") %in% names(out))
  )
  expect_equal(out$hypothesis, rep("dextrt = 0", length(features)))

  # The metadata the band came from is what the fits actually kept.
  expect_true(all(vapply(fits, function(f) brms::ndraws(f) == 100, logical(1))))

  # Bundles keep gene order, so unpacking them lines the fits up with `.feature`.
  counts <- SummarizedExperiment::assay(se, "counts")
  expect_equal(
    lapply(fits, function(f) as.numeric(f$data$counts)),
    lapply(features, function(g) as.numeric(counts[g, ]))
  )
})

test_that("estimate |> hypothesis |> adjust evaluate as one pipeline", {
  skip_on_cran()
  skip_if_not_installed("tidytargets")
  skip_if_no_cmdstan()
  skip_if_not_installed("targets")

  se <- airway_for_hpc()

  store <- tempfile("brmde_pipe_eval_")
  on.exit(
    {
      unlink(store, recursive = TRUE)
      unlink(paste0(store, ".R"))
    },
    add = TRUE
  )

  out <- suppressWarnings(
    se |>
      brmDE(
        features = c("ENSG00000120129"),
        store = store
      ) |>
      estimate(
        ~ dex + (1 | cell),
        offset = "offset",
        dispersion_prior_log_mean = "dispersion_trended",
        family = brms::negbinomial(),
        chains = 1,
        draws_warmup = 100,
        draws_sampling = 150,
        cores = 1,
        backend = "cmdstanr",
        refresh = 0,
        silent = 2
      ) |>
      hypothesis("dextrt = 0") |>
      adjust(nullify = "dex") |>
      tt_evaluate()
  )

  expect_equal(out$.feature, "ENSG00000120129")
  # The fit is not in the result, but it is in the store. A branch holds a
  # bundle, so reading one gives a list of fits even at bundle = 1.
  expect_false("brms_fit" %in% names(out))
  branch <- targets::tar_read_raw("brms_fit", branches = 1L, store = store)
  expect_s3_class(unlist(branch, recursive = FALSE)[[1]], "brmsfit")

  # The hypothesis columns are inline; adjust stays a per-gene list column.
  expect_equal(out$hypothesis, "dextrt = 0")
  expect_true(is.finite(out$estimate))
  expect_s3_class(out$adjust[[1]], "tbl_df")
  expect_true("adjusted___Estimate" %in% names(out$adjust[[1]]) ||
    any(grepl("^adjusted___", names(out$adjust[[1]]))))
})
