test_that("estimate, hypothesis, and adjust append to one tidytargets script", {
  skip_if_not_installed("tidytargets")
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
    estimate(~ dex, offset = "offset", dispersion = "dispersion",
             dispersion_degrees_freedom = "dispersion_degrees_freedom") |>
    hypothesis("dextrt = 0") |>
    adjust(nullify = "dex")

  expect_s3_class(pipeline, "brmDE_hpc")
  expect_s3_class(pipeline, "tidytargets")
  expect_true("brms_fit" %in% names(pipeline))
  expect_true("hypothesis_tbl" %in% names(pipeline))
  expect_true("adjust_tbl" %in% names(pipeline))

  lines <- readLines(paste0(store, ".R"))
  script <- paste(lines, collapse = "\n")
  expect_match(script, "tt_factory")
  # Arguments are targets, so targets invalidates the fits when they change and
  # no worker reads them off disk. Only the input object itself is a file.
  expect_match(script, 'target_output = "brms_fit_args"')
  expect_match(script, 'target_output = "hypothesis_tbl_args"')
  expect_match(script, 'target_output = "adjust_tbl_args"')
  expect_match(script, 'offset = "offset"', fixed = TRUE)
  expect_setequal(
    list.files(paste0(store, "_input")),
    c("se.rds", "computing_resources.rds")
  )
  expect_match(script, 'target_output = "formula_abundance_text"')
  expect_match(script, 'target_output = "formula_dispersion_text"')
  # Setup / identity targets run on the main process, not the crew controller.
  expect_match(script, 'deployment = "main"')
  # Dispersion is estimated upstream of the pipeline, not inside it.
  expect_no_match(script, "estimate_dispersion", fixed = TRUE)
  expect_match(script, "estimate_gene_from_se")
  expect_match(script, "hypothesis_gene_from_fit")
  expect_match(script, "adjust_gene_from_fit")
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
      unlink(paste0(s, "_input"), recursive = TRUE)
    },
    add = TRUE
  )

  # targets hashes the target command, so anything that should force a refit
  # has to be visible in the generated script. Arguments passed through `...`
  # get there as the code they were written as, hence the quoted family.
  script_for <- function(formula_abundance = ~dex,
                         formula_dispersion = ~1,
                         features = "ENSG00000120129",
                         family = quote(brms::negbinomial())) {
    store <- tempfile("brmde_hash_")
    stores <<- c(stores, store)
    eval(bquote(
      se |>
        brmDE(store = store, features = features) |>
        estimate(
          formula_abundance,
          formula_dispersion = formula_dispersion,
          offset = "offset",
          dispersion = "dispersion",
          dispersion_degrees_freedom = "dispersion_degrees_freedom",
          family = .(family)
        )
    ))
    # gsub, not sub: a single line can mention the store more than once.
    gsub(basename(store), "STORE", readLines(paste0(store, ".R")), fixed = TRUE)
  }

  base <- script_for()
  expect_true(any(grepl("brms::negbinomial()", base, fixed = TRUE)))
  expect_identical(script_for(), base)
  expect_false(identical(script_for(formula_abundance = ~ dex + (1 | cell)), base))
  expect_false(identical(script_for(formula_dispersion = ~cell), base))
  expect_false(identical(script_for(features = "ENSG00000000003"), base))
  expect_false(
    identical(
      script_for(family = quote(brms::zero_inflated_negbinomial())),
      base
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
      unlink(paste0(store, "_input"), recursive = TRUE)
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
      unlink(paste0(s, "_input"), recursive = TRUE)
    },
    add = TRUE
  )

  script_for <- function(...) {
    store <- tempfile("brmde_bundle_")
    stores <<- c(stores, store)
    se |>
      brmDE(store = store, features = rownames(se)[1:6]) |>
      estimate(~dex, offset = "offset", ...) |>
      hypothesis("dextrt = 0")
    # The store path is baked into the header, so blank it out before compared
    # scripts can be expected to match.
    gsub(basename(store), "STORE", readLines(paste0(store, ".R")), fixed = TRUE)
  }

  plain <- script_for()
  bundled <- script_for(bundle = 3)

  # The default is the graph as it was before bundling existed: one fit target
  # per gene, and no extra target in between.
  expect_false(any(grepl("gene_bundle", plain, fixed = TRUE)))
  expect_true(any(grepl('other_arguments_to_map = "gene_id"', plain, fixed = TRUE)))
  expect_identical(script_for(bundle = 1), plain)

  # With bundling the fits map over the bundles instead, and `bundle` is in
  # the command so that changing it invalidates the fits.
  expect_true(any(grepl('target_output = "gene_bundle"', bundled, fixed = TRUE)))
  expect_true(any(grepl('other_arguments_to_map = "gene_bundle"', bundled, fixed = TRUE)))
  expect_true(any(grepl("bundle = 3L", bundled, fixed = TRUE)))

  # hypothesis() maps over the fit target, so it inherits the coarser branching.
  expect_true(any(grepl('other_arguments_to_map = "brms_fit"', bundled, fixed = TRUE)))
})

test_that("bundle_gene_ids groups genes by size and keeps their order", {
  ids <- as.list(paste0("g", 1:10))

  expect_identical(unlist(bundle_gene_ids(ids, 4)), unlist(ids))
  # Bundles are `bundle` genes each, so only the last one is a remainder.
  expect_identical(lengths(bundle_gene_ids(ids, 4)), c(4L, 4L, 2L))
  expect_identical(lengths(bundle_gene_ids(ids, 5)), c(5L, 5L))
  # The default is one target per gene.
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
      unlink(paste0(store, "_input"), recursive = TRUE)
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
      unlink(paste0(store, "_input"), recursive = TRUE)
    },
    add = TRUE
  )

  out <- suppressWarnings(
    se |>
      brmDE(features = features, store = store) |>
      estimate(
        ~dex,
        offset = "offset",
        dispersion = "dispersion",
        dispersion_degrees_freedom = "dispersion_degrees_freedom",
        bundle = 2,
        family = brms::negbinomial(),
        chains = 1,
        iter = 200,
        warmup = 100,
        cores = 1,
        backend = "cmdstanr",
        refresh = 0,
        silent = 2
      ) |>
      hypothesis("dextrt = 0") |>
      tt_evaluate()
  )

  # Four genes, but only two fit targets were built.
  branches <- targets::tar_read_raw("brms_fit", store = store)
  expect_length(branches, 2L)
  fits <- unlist(branches, recursive = FALSE)
  expect_true(all(vapply(fits, inherits, logical(1), "brmsfit")))

  expect_equal(out$.feature, features)
  expect_true(all(vapply(out$hypothesis, inherits, logical(1), "tbl_df")))

  # estimate() contributes the gene column and nothing else: the fits are far
  # too large to gather, and the statistics that make the run readable are the
  # hypothesis ones.
  expect_named(out, c(".feature", "hypothesis"))
  expect_true(all(vapply(
    out$hypothesis,
    function(h) all(c("rhat", "ess_bulk", "mcse") %in% names(h)),
    logical(1)
  )))

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
      unlink(paste0(store, "_input"), recursive = TRUE)
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
        dispersion = "dispersion",
        dispersion_degrees_freedom = "dispersion_degrees_freedom",
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
      tt_evaluate()
  )

  expect_equal(out$.feature, "ENSG00000120129")
  expect_false("brms_fit" %in% names(out))
  expect_s3_class(
    targets::tar_read_raw("brms_fit", branches = 1L, store = store)[[1]],
    "brmsfit"
  )
  expect_s3_class(out$hypothesis[[1]], "tbl_df")
  expect_s3_class(out$adjust[[1]], "tbl_df")
  expect_true("adjusted___Estimate" %in% names(out$adjust[[1]]) ||
    any(grepl("^adjusted___", names(out$adjust[[1]]))))
})
