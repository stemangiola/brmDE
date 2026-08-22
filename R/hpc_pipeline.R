# tidytargets grammar wrappers. The gene-wise engines stay in estimate_gene(),
# hypothesis_gene(), and adjust_gene(); these methods only append them to a
# tidytargets graph via tt_single() / tt_iterate().

as_brmde_hpc <- function(x) {
  class(x) <- unique(c("brmDE_hpc", "tidytargets", class(x)))
  x
}

brmde_input_dir <- function(store) {
  paste0(store, "_input")
}

# Formulas cross into the targets graph as text, so they carry no environment
# and compare cleanly between runs. Terms resolve against the data first and
# the global environment after, matching estimate_gene().
as_pipeline_formula <- function(text) {
  stats::as.formula(text, env = globalenv())
}

write_brmde_hpc_header <- function(script,
                                   packages,
                                   controller_rds,
                                   debug_step,
                                   update,
                                   garbage_collection,
                                   workspace_on_error) {
  expr <- substitute(
    {
      library(tidytargets)
      library(brmDE)
      library(targets)

      tar_option_set(
        memory = "transient",
        garbage_collection = GARBAGE,
        storage = "main",
        retrieval = "main",
        debug = DEBUG,
        cue = tar_cue(mode = UPDATE),
        controller = readRDS(CONTROLLER),
        packages = PACKAGES,
        workspace_on_error = WORKSPACE,
        format = "rds"
      )

      target_list <- list()
    },
    list(
      GARBAGE = garbage_collection,
      DEBUG = debug_step,
      UPDATE = update,
      WORKSPACE = workspace_on_error,
      CONTROLLER = controller_rds,
      PACKAGES = packages
    )
  )

  lines <- deparse(expr, width.cutoff = 500L)
  lines <- lines[-1]
  lines <- lines[-length(lines)]
  writeLines(lines, script)
}

check_features <- function(features) {
  if (is.null(features)) {
    return(features)
  }
  if (is.symbol(features) || is.name(features) || is.language(features)) {
    stop(
      "`features` must be a character vector of gene ids, e.g. c(\"gene1\", \"gene2\").",
      call. = FALSE
    )
  }
  if (!is.character(features) || length(features) < 1L || anyNA(features) || !all(nzchar(features))) {
    stop(
      "`features` must be a character vector of gene ids, e.g. c(\"gene1\", \"gene2\").",
      call. = FALSE
    )
  }
  as.character(features)
}

#' @keywords internal
#' @export
gene_ids_for_hpc <- function(se, features = NULL) {
  ids <- rownames(se)
  if (!is.null(features)) {
    ids <- intersect(ids, features)
  }
  if (length(ids) == 0L) {
    stop("No genes left to estimate after applying `features`.", call. = FALSE)
  }
  as.list(ids)
}

# One branch per element, so this list is what sets the number of HPC jobs.
# Splitting genes into fewer, larger elements trades scheduler overhead for
# coarser invalidation: the ids in a bundle are hashed together, so changing
# the gene set reshuffles the bundles and refits all of them.
#' @keywords internal
#' @export
bundle_gene_ids <- function(ids, bundle) {
  ids <- unlist(ids, use.names = FALSE)
  unname(split(ids, ceiling(seq_along(ids) / max(1L, as.integer(bundle)))))
}

# Every branch returns a list of gene-wise results, one element per gene it was
# given, whether or not genes were bundled. collect_branches() relies on that.
#' @keywords internal
#' @export
estimate_gene_from_se <- function(se,
                                  feature_id,
                                  formula_abundance,
                                  formula_dispersion,
                                  args) {
  lapply(unlist(feature_id), function(id) {
    do.call(
      estimate_gene,
      c(
        list(
          data = se[id, , drop = FALSE],
          formula_abundance = as_pipeline_formula(formula_abundance),
          formula_dispersion = as_pipeline_formula(formula_dispersion)
        ),
        args
      )
    )
  })
}

#' @keywords internal
#' @export
hypothesis_gene_from_fit <- function(fit, args) {
  lapply(as_branch_list(fit), function(f) {
    do.call(hypothesis_gene, c(list(fit = f), args))
  })
}

#' @keywords internal
#' @export
adjust_gene_from_fit <- function(fit, args) {
  lapply(as_branch_list(fit), function(f) {
    do.call(adjust_gene, c(list(fit = f), args))
  })
}

as_branch_list <- function(x) {
  if (inherits(x, c("brmsfit", "data.frame"))) {
    return(list(x))
  }
  if (is.list(x)) {
    return(x)
  }
  list(x)
}

# A branch holds one bundle's worth of results, so drop that level to get back
# to one element per gene. Branches come back in bundle order and bundles keep
# gene order, so this lines up with `gene_id`.
collect_branches <- function(x) {
  do.call(c, unname(lapply(as_branch_list(x), as_branch_list)))
}

collect_brmde_hpc <- function(input_hpc) {
  store <- input_hpc$initialisation$store
  gene_id <- unlist(
    targets::tar_read_raw("gene_id", store = store),
    use.names = FALSE
  )
  out <- tibble::tibble(.feature = as.character(gene_id))

  if ("brms_fit" %in% names(input_hpc)) {
    out$brms_fit <- collect_branches(
      targets::tar_read_raw("brms_fit", store = store)
    )
  }
  if ("hypothesis_tbl" %in% names(input_hpc)) {
    out$hypothesis <- collect_branches(
      targets::tar_read_raw("hypothesis_tbl", store = store)
    )
  }
  if ("adjust_tbl" %in% names(input_hpc)) {
    out$adjust <- collect_branches(
      targets::tar_read_raw("adjust_tbl", store = store)
    )
  }
  out
}

#' Start a gene-wise tidytargets pipeline
#'
#' Analogue of tidytargets' [tidytargets::tt_initialise()]: writes the targets
#' header (`target_list`) and the shared steps (load SE, gene ids). Later
#' calls to [estimate()], [hypothesis()], and [adjust()] append
#' `tt_iterate()` branches onto the same graph. Printing the object (or
#' calling [tt_evaluate()]) runs the pipeline, as in tidytargets.
#'
#' The pipeline is gene-wise only. Anything estimated across the whole matrix,
#' namely the library size offset and the dispersion, belongs upstream of
#' [brmDE()] and arrives as ordinary `colData` / `rowData` columns.
#'
#' @param .data A `SummarizedExperiment` that already carries a library size
#'   offset in `colData` and dispersion in `rowData`. Neither is
#'   computed here: both are whole-matrix quantities, so derive them before
#'   the pipeline (`tidybulk::scale_abundance()` taking `log(1 / multiplier)`,
#'   and [tidybulk::estimate_dispersion()]) and name the columns when you call
#'   [estimate()].
#' @param abundance Count assay name.
#' @param features Optional character vector of gene ids to fit, e.g.
#'   `c("gene1", "gene2")`. Subsetting here only limits which genes are
#'   fitted; the offset and dispersion already reflect all genes.
#' @param store Directory used as the targets store. The script is
#'   `{store}.R`.
#' @param computing_resources A crew controller, or `NULL` for sequential
#'   targets (default). On SLURM, use
#'   `crew.cluster::crew_controller_slurm()`; see the vignette section
#'   *SLURM* for an example.
#' @param debug_step Optional target name for `tar_option_set(debug = )`.
#' @param callr_function Passed to [targets::tar_make()] by
#'   [tt_evaluate()]. `NULL` (default) runs in the current session.
#' @param packages Character vector of packages loaded in workers.
#' @param update,garbage_collection,workspace_on_error Passed through
#'   to `tar_option_set()` / `tar_cue()` as in tidytargets.
#' @param verbosity Reporter passed to [targets::tar_make()].
#'
#' @return A `tidytargets` pipeline object (`brmDE_hpc`).
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#'
#' # Whole-matrix quantities are prepared before the gene-wise pipeline.
#' se <- airway
#' se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
#' se <- tidybulk::estimate_dispersion(se, formula_abundance = ~ dex + cell)
#'
#' se |>
#'   brmDE(features = c("ENSG00000120129")) |>
#'   estimate(
#'     ~ dex + (1 | cell),
#'     offset = "offset",
#'     dispersion = "dispersion",
#'     dispersion_degrees_freedom = "dispersion_degrees_freedom",
#'     family = brms::negbinomial()
#'   ) |>
#'   hypothesis("dextrt = 0") |>
#'   adjust(nullify = "dex")
#'
#' # SLURM (from a cluster login node; requires crew.cluster):
#' # library(crew.cluster)
#' # se |>
#' #   brmDE(
#' #     store = "/path/to/brmde_store",
#' #     computing_resources = crew_controller_slurm(
#' #       workers = 100,
#' #       options_cluster = crew_options_slurm(
#' #         memory_gigabytes_per_cpu = 5
#' #       )
#' #     )
#' #   ) |>
#' #   estimate(~ dex + (1 | cell), offset = "offset", dispersion = "dispersion",
#' #            dispersion_degrees_freedom = "dispersion_degrees_freedom")
#' }
#'
#' @export
brmDE <- function(.data,
                      abundance = "counts",
                      features = NULL,
                      store = tempfile(tmpdir = tempdir(), pattern = "brmde_"),
                      computing_resources = NULL,
                      debug_step = NULL,
                      callr_function = NULL,
                      packages = c(
                        "tidytargets",
                        "brmDE",
                        "SummarizedExperiment",
                        "brms",
                        "tibble"
                      ),
                      update = "thorough",
                      garbage_collection = TRUE,
                      workspace_on_error = FALSE,
                      verbosity = targets::tar_config_get("reporter_make")) {
  if (!inherits(.data, "SummarizedExperiment")) {
    stop("`.data` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (nrow(.data) < 1L) {
    stop("`.data` has no genes.", call. = FALSE)
  }
  features <- check_features(features)
  if (!requireNamespace("tidytargets", quietly = TRUE)) {
    stop("Install tidytargets to build the pipeline.", call. = FALSE)
  }

  dir.create(store, showWarnings = FALSE, recursive = TRUE)
  input_dir <- brmde_input_dir(store)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  script <- paste0(store, ".R")

  se_rds <- normalizePath(file.path(input_dir, "se.rds"), mustWork = FALSE)
  controller_rds <- normalizePath(
    file.path(input_dir, "computing_resources.rds"),
    mustWork = FALSE
  )
  saveRDS(.data, se_rds)
  saveRDS(computing_resources, controller_rds)

  write_brmde_hpc_header(
    script = script,
    packages = packages,
    controller_rds = controller_rds,
    debug_step = debug_step,
    update = update,
    garbage_collection = garbage_collection,
    workspace_on_error = workspace_on_error
  )

  input_hpc <- as_brmde_hpc(
    list(
      initialisation = list(
        store = store,
        computing_resources = computing_resources,
        tier = 1L,
        debug_step = debug_step,
        verbosity = verbosity,
        abundance = abundance,
        features = features,
        callr_function = callr_function
      )
    )
  )

  input_hpc |>
    # Cheap setup targets stay on the main process; only gene-wise branches
    # go through the crew controller.
    tidytargets::tt_single(
      "file_se",
      se_rds,
      format = "file",
      deployment = "main"
    ) |>
    tidytargets::tt_single(
      target_output = "se_input",
      user_function = readRDS |> quote(),
      file = "file_se" |> tidytargets::is_target(),
      deployment = "main"
    ) |>
    tidytargets::tt_single(
      target_output = "gene_id",
      user_function = gene_ids_for_hpc |> quote(),
      se = "se_input" |> tidytargets::is_target(),
      features = features,
      iterate = "map",
      deployment = "main"
    ) |>
    as_brmde_hpc()
}

#' Estimate genes on a tidytargets pipeline
#'
#' Appends a [tidytargets::tt_iterate()] step that calls [estimate_gene()] once
#' per gene.
#'
#' @param input_hpc A `tidytargets` / `brmDE_hpc` pipeline from [brmDE()].
#' @param formula_abundance Model for the mean, passed to [estimate_gene()].
#' @param formula_dispersion One-sided model for the negative binomial shape,
#'   passed to [estimate_gene()]. `~1` by default.
#' @param offset Required name of the precomputed offset column in `colData`
#'   (the same argument as [estimate_gene()]).
#' @param dispersion Optional name of the `rowData` dispersion column, as
#'   written by [tidybulk::estimate_dispersion()]. Default `NULL` puts a zero
#'   offset on the shape submodel (see [estimate_gene()]).
#' @param dispersion_degrees_freedom Optional name of the effective degrees of
#'   freedom column written alongside it by [tidybulk::estimate_dispersion()].
#'   Default `NULL` uses a default Student-t SD of 1 on the log-shape intercept.
#'
#' @details
#' Neither the offset nor the dispersion is computed here. Run
#' [tidybulk::estimate_dispersion()] on the whole object before [brmDE()],
#' exactly as you compute the offset, so that both are ordinary columns of the
#' input, then pass those column names here. Gene-wise fitting is the only
#' thing this pipeline does. The dispersion columns are a prior on the
#' negative binomial shape (see [estimate_gene()]): they inform each gene
#' without fixing the posterior to the external estimate. Omitting
#' `dispersion` and/or `dispersion_degrees_freedom` is valid (zero shape
#' offset, default prior SD) but computing and passing both is the preferred
#' starting point.
#'
#' Prior constants derived from each gene are passed to Stan as data, so every
#' gene generates identical Stan code and cmdstanr compiles it at most once per
#' worker process rather than once per gene.
#' @param bundle Number of genes to fit per target. `1` (default) gives one
#'   target per gene. A larger value fits that many genes sequentially inside
#'   one target, which is how you stop tens of thousands of genes from
#'   swamping an HPC scheduler with tiny jobs. See *Bundling* below.
#' @param target_output Name of the targets output.
#' @param ... Passed to [estimate_gene()] (e.g. `family`, `chains`, `iter`).
#'
#' @return The updated `tidytargets` pipeline.
#'
#' @section Bundling:
#' One target per gene gives the scheduler tens of thousands of jobs whose
#' queueing cost rivals the fit itself. `bundle` puts that many genes in each
#' fit target instead, and [hypothesis()] and [adjust()] inherit the coarser
#' branching because they map over the fit target. The result of
#' [tt_evaluate()] is unchanged: one row per gene, in the same order.
#'
#' Size `bundle` against one job rather than against the gene count, since
#' that is what it controls.
#'
#' * A bundle is fitted sequentially in one worker, so it needs about `bundle`
#'   times the walltime of a single fit.
#' * A bundle is one `.rds` holding `bundle` fits, and it is read back into
#'   the main process, so the size of a `brmsfit` is what limits how far you
#'   can push this.
#' * Invalidation is per bundle. One gene failing loses its bundle, and
#'   changing the gene set reshuffles bundle membership and refits everything.
#'   Leave `bundle` at `1` when you rely on incremental reruns.
#'
#' @export
estimate <- function(input_hpc,
                     formula_abundance,
                     formula_dispersion = ~1,
                     offset,
                     dispersion = NULL,
                     dispersion_degrees_freedom = NULL,
                     ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
estimate.default <- function(input_hpc,
                             formula_abundance,
                             formula_dispersion = ~1,
                             offset,
                             dispersion = NULL,
                             dispersion_degrees_freedom = NULL,
                             ...) {
  stop(
    "estimate() expects a pipeline from brmDE(). ",
    "For one gene, use estimate_gene().",
    call. = FALSE
  )
}

#' @rdname estimate
#' @export
estimate.tidytargets <- function(input_hpc,
                            formula_abundance,
                            formula_dispersion = ~1,
                            offset,
                            dispersion = NULL,
                            dispersion_degrees_freedom = NULL,
                            bundle = 1L,
                            target_output = "brms_fit",
                            ...) {
  offset <- check_offset_name(offset)
  bundle <- check_bundle(bundle)
  if (!is.null(dispersion)) {
    dispersion <- check_dispersion_name(dispersion)
  }
  if (!is.null(dispersion_degrees_freedom)) {
    dispersion_degrees_freedom <-
      check_degrees_freedom_name(dispersion_degrees_freedom)
  }
  abundance <- input_hpc$initialisation$abundance
  args_target <- paste0(target_output, "_args")

  # The arguments travel as the code they were written as, so targets hashes
  # that code and the worker rebuilds the objects rather than reading them off
  # disk. As anywhere in targets, an argument naming something only this session
  # holds cannot be rebuilt there: pass `brms::negbinomial()`, not a variable
  # standing for it.
  args <- as.call(c(
    quote(list),
    list(
      abundance = abundance,
      offset = offset,
      dispersion = dispersion,
      dispersion_degrees_freedom = dispersion_degrees_freedom
    ),
    as.list(substitute(list(...)))[-1L]
  ))

  # The formulas and the arguments are targets of their own so that editing one
  # invalidates the fits that depend on it. They only assemble small objects, so
  # they run on main.
  pipeline <- input_hpc |>
    tidytargets::tt_single(
      target_output = "formula_abundance_text",
      user_function = identity |> quote(),
      x = formula_text(formula_abundance),
      deployment = "main"
    ) |>
    tidytargets::tt_single(
      target_output = "formula_dispersion_text",
      user_function = identity |> quote(),
      x = formula_text(formula_dispersion),
      deployment = "main"
    ) |>
    tidytargets::tt_single(
      target_output = args_target,
      user_function = identity |> quote(),
      x = call("quote", args),
      deployment = "main"
    )

  # The branch count is whatever this target's list is long, so bundling is
  # just a regrouping of the gene ids upstream of the fit. At bundle = 1 that
  # regrouping is the identity, and leaving the target out keeps the graph
  # (and so the stored hashes) as it was before bundling existed.
  feature_target <- "gene_id"
  if (bundle > 1L) {
    feature_target <- "gene_bundle"
    pipeline <- pipeline |>
      tidytargets::tt_single(
        target_output = feature_target,
        user_function = bundle_gene_ids |> quote(),
        ids = "gene_id" |> tidytargets::is_target(),
        bundle = bundle,
        iterate = "map",
        deployment = "main"
      )
  }

  pipeline |>
    tidytargets::tt_iterate(
      target_output = target_output,
      user_function = estimate_gene_from_se |> quote(),
      se = "se_input" |> tidytargets::is_target(),
      feature_id = feature_target |> tidytargets::is_target(),
      formula_abundance = "formula_abundance_text" |> tidytargets::is_target(),
      formula_dispersion = "formula_dispersion_text" |> tidytargets::is_target(),
      args = args_target |> tidytargets::is_target()
    ) |>
    as_brmde_hpc()
}

#' Hypothesis tests on a tidytargets pipeline
#'
#' Appends a [tidytargets::tt_iterate()] step that calls [hypothesis_gene()]
#' on each [estimate()] fit. This is an S3 method for
#' [brms::hypothesis()].
#'
#' @param x A `tidytargets` pipeline from [brmDE()].
#' @param hypothesis Hypothesis strings passed to [hypothesis_gene()].
#' @param target_input Name of the upstream fit target.
#' @param target_output Name of the targets output.
#' @param ... Passed to [hypothesis_gene()].
#'
#' @return The updated `tidytargets` pipeline.
#'
#' @exportS3Method brms::hypothesis
hypothesis.tidytargets <- function(x,
                              hypothesis,
                              target_input = "brms_fit",
                              target_output = "hypothesis_tbl",
                              ...) {
  if (!target_input %in% names(x)) {
    stop(
      "hypothesis() needs estimate() on the pipeline first ",
      "(missing target '", target_input, "').",
      call. = FALSE
    )
  }

  args_target <- paste0(target_output, "_args")
  args <- as.call(c(
    quote(list),
    list(hypothesis = hypothesis),
    as.list(substitute(list(...)))[-1L]
  ))

  x |>
    tidytargets::tt_single(
      target_output = args_target,
      user_function = identity |> quote(),
      x = call("quote", args),
      deployment = "main"
    ) |>
    tidytargets::tt_iterate(
      target_output = target_output,
      user_function = hypothesis_gene_from_fit |> quote(),
      fit = target_input |> tidytargets::is_target(),
      args = args_target |> tidytargets::is_target()
    ) |>
    as_brmde_hpc()
}

#' Adjust genes on a tidytargets pipeline
#'
#' Appends a [tidytargets::tt_iterate()] step that calls [adjust_gene()] on
#' each [estimate()] fit.
#'
#' @param input_hpc A `tidytargets` pipeline, or a `brmsfit`.
#' @param nullify Covariates to nullify, passed to [adjust_gene()].
#' @param target_input Name of the upstream fit target.
#' @param target_output Name of the targets output.
#' @param ... Passed to [adjust_gene()].
#'
#' @return The updated `tidytargets` pipeline, or an adjustment tibble.
#'
#' @export
adjust <- function(input_hpc, ...) {
  UseMethod("adjust")
}

#' @rdname adjust
#' @export
adjust.default <- function(input_hpc, ...) {
  if (inherits(input_hpc, "brmsfit")) {
    return(adjust_gene(input_hpc, ...))
  }
  stop(
    "adjust() expects a pipeline from brmDE(), or a brmsfit.",
    call. = FALSE
  )
}

#' @rdname adjust
#' @export
adjust.tidytargets <- function(input_hpc,
                          nullify = NULL,
                          target_input = "brms_fit",
                          target_output = "adjust_tbl",
                          ...) {
  if (!target_input %in% names(input_hpc)) {
    stop(
      "adjust() needs estimate() on the pipeline first ",
      "(missing target '", target_input, "').",
      call. = FALSE
    )
  }

  args_target <- paste0(target_output, "_args")
  args <- as.call(c(
    quote(list),
    list(nullify = nullify),
    as.list(substitute(list(...)))[-1L]
  ))

  input_hpc |>
    tidytargets::tt_single(
      target_output = args_target,
      user_function = identity |> quote(),
      x = call("quote", args),
      deployment = "main"
    ) |>
    tidytargets::tt_iterate(
      target_output = target_output,
      user_function = adjust_gene_from_fit |> quote(),
      fit = target_input |> tidytargets::is_target(),
      args = args_target |> tidytargets::is_target()
    ) |>
    as_brmde_hpc()
}

#' @export
#' @importFrom tidytargets tt_evaluate
tidytargets::tt_evaluate

#' @export
#' @importFrom brms hypothesis
brms::hypothesis

localize_target_append <- function(script) {
  lines <- readLines(script)
  lines <- vapply(
    lines,
    function(line) {
      if (!grepl("target_list \\|> target_append\\(", line, perl = TRUE)) {
        return(line)
      }
      inner <- sub("^.*target_list \\|> target_append\\(", "", line, perl = TRUE)
      inner <- sub("\\)\\s*$", "", inner, perl = TRUE)
      paste0("target_list <- c(target_list, list(", inner, "))")
    },
    character(1),
    USE.NAMES = FALSE
  )
  writeLines(lines, script)
}

#' @exportS3Method tidytargets::tt_evaluate
tt_evaluate.brmDE_hpc <- function(tt_input) {
  store <- tt_input$initialisation$store
  script <- paste0(store, ".R")
  localize_target_append(script)
  cat("target_list\n", file = script, append = TRUE)

  reporter <- tt_input$initialisation$verbosity
  make_args <- list(
    script = script,
    store = store,
    callr_function = tt_input$initialisation$callr_function
  )
  if (!is.null(reporter) && !identical(reporter, "undefined")) {
    make_args$reporter <- reporter
  }
  do.call(targets::tar_make, make_args)

  collect_brmde_hpc(tt_input)
}

#' @rdname brmDE
#' @export
print.brmDE_hpc <- function(x, ...) {
  x |>
    tt_evaluate() |>
    print(...)
}
