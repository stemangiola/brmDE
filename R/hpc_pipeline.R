#' Mark a tidytargets graph as a `brmDE` pipeline.
#'
#' @param x A tidytargets pipeline list.
#'
#' @return `x` with classes `brmDE_hpc` and `tidytargets`.
#'
#' @keywords internal
#' @noRd
as_brmde_hpc <- function(x) {
  class(x) <- unique(c("brmDE_hpc", "tidytargets", class(x)))
  x
}

# Formulas cross into the targets graph as text, so they carry no environment
# and compare cleanly between runs. Terms resolve against the data first and
# the global environment after, matching estimate_gene().
as_pipeline_formula <- function(text) {
  stats::as.formula(text, env = globalenv())
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

#' Gene ids for HPC mapping
#'
#' @param se A `SummarizedExperiment`.
#' @param features Optional character vector of gene ids.
#'
#' @return A list of gene ids, one per target branch.
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
#' Bundle gene ids into HPC jobs
#'
#' @param ids Gene ids, typically a list from [gene_ids_for_hpc()].
#' @param bundle Number of genes per branch.
#'
#' @return A list of character vectors of gene ids.
#' @keywords internal
#' @export
bundle_gene_ids <- function(ids, bundle) {
  ids <- unlist(ids, use.names = FALSE)
  unname(split(ids, ceiling(seq_along(ids) / max(1L, as.integer(bundle)))))
}

# Every branch returns a list of gene-wise results, one element per gene it was
# given, whether or not genes were bundled. collect_branches() relies on that.
#' Fit one HPC branch of genes
#'
#' @param se A `SummarizedExperiment`.
#' @param feature_id Gene id or bundle of ids for this branch.
#' @param formula_abundance,formula_dispersion Formulas stored as text on the
#'   pipeline.
#' @param args Named list of extra [estimate_gene()] arguments.
#'
#' @return A list of `brmsfit` objects, one per gene in the branch.
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

#' Test hypotheses for one HPC branch
#'
#' @param fit A `brmsfit`, or a list of them from one branch.
#' @param args Named list of extra [hypothesis_gene()] arguments.
#'
#' @return A list of hypothesis tables, one per gene in the branch.
#' @keywords internal
#' @export
hypothesis_gene_from_fit <- function(fit, args) {
  lapply(as_branch_list(fit), function(f) {
    do.call(hypothesis_gene, c(list(fit = f), args))
  })
}

#' Adjust one HPC branch
#'
#' @param fit A `brmsfit`, or a list of them from one branch.
#' @param args Named list of extra [adjust_gene()] arguments.
#'
#' @return A list of adjusted tables, one per gene in the branch.
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

  # Twenty thousand brmsfit objects do not fit in one table, so estimate()
  # contributes the gene column alone. The fits stay in the store; read one
  # with targets::tar_read(brms_fit, branches = i, store = store).
  if ("hypothesis_tbl" %in% names(input_hpc)) {
    # A false discovery rate is a property of the gene set, so it is the one
    # quantity that cannot be computed inside a gene-wise target. The tables
    # are tiny, so it is added here, where every gene is in hand at once.
    out$hypothesis <- add_hypothesis_fdr(
      collect_branches(targets::tar_read_raw("hypothesis_tbl", store = store))
    )
     out <- tidyr::unnest(out, cols = "hypothesis")

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
#' Calls [tidytargets::tt_initialise()] for the script header, then
#' [tidytargets::tt_data()] to snapshot the `SummarizedExperiment` onto the
#' store and [tidytargets::tt_data_list()] for the gene ids to map over.
#' Later calls to [estimate()], [hypothesis()], and [adjust()] append
#' `tt_iterate()` branches onto the same graph. Printing the object (or
#' calling [tt_evaluate()]) runs the pipeline, as in tidytargets. Assigning
#' it does not; an interactive session then says the pipeline is ready to be
#' evaluated, rather than appearing to do nothing.
#'
#' The pipeline is gene-wise only. Anything estimated across the whole matrix,
#' namely the library size offset and the dispersion, belongs upstream of
#' [brmDE()] and arrives as ordinary `colData` / `rowData` columns.
#'
#' Running the pipeline gives one row per gene and contrast. [estimate()]
#' contributes the gene column; [hypothesis()] adds its columns inline, each
#' carrying the convergence of that contrast, and [adjust()] its own table.
#' The fits are not in that
#' result, because a transcriptome's worth of `brmsfit` objects will not fit
#' in one table: they stay in the store. Read one with
#' `targets::tar_read(brms_fit, branches = i, store = store)`, using the same
#' `store` you passed to [brmDE()], or
#' `tidytargets::tt_explore(pipeline, "brms_fit", index = i)`.
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
#' @return A `tidytargets` pipeline object (`brmDE_hpc`). The graph is not
#'   run until you print it or call [tt_evaluate()].
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#'
#' # Whole-matrix quantities are prepared before the gene-wise pipeline.
#' se <- airway
#' se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
#' se <- tidybulk::estimate_dispersion(se, formula_abundance = ~ dex + cell)
#' se <- estimate_dispersion_prior(se, formula_abundance = ~ dex + cell)
#'
#' se |>
#'   brmDE(features = c("ENSG00000120129")) |>
#'   estimate(
#'     ~ dex + (1 | cell),
#'     offset = "offset",
#'     dispersion_prior_log_mean = "dispersion_prior_log_mean",
#'     dispersion_prior_log_sd = "dispersion_prior_log_sd",
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
#' #   estimate(~ dex + (1 | cell), offset = "offset",
#' #            dispersion_prior_log_mean = "dispersion_prior_log_mean",
#' #            dispersion_prior_log_sd = "dispersion_prior_log_sd")
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

  # Mapped units are registered by tt_data_list() now, not by
  # tt_single(iterate = "map"), so the gene ids are snapshotted here.
  gene_ids <- gene_ids_for_hpc(.data, features)

  tidytargets::tt_initialise(
    store = store,
    computing_resources = computing_resources,
    debug_step = debug_step,
    verbosity = verbosity,
    update = update,
    garbage_collection = garbage_collection,
    workspace_on_error = workspace_on_error,
    packages = packages
  ) |>
    tidytargets::tt_metadata(
      abundance = abundance,
      features = features,
      callr_function = callr_function,
      gene_ids = gene_ids
    ) |>
    tidytargets::tt_data(.data, target_output = "se_input") |>
    tidytargets::tt_data_list(gene_ids, target_output = "gene_id") |>
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
#' @param dispersion_prior_log_mean Optional name of the `rowData` dispersion
#'   column. [estimate_dispersion_prior()] writes `"dispersion_prior_log_mean"`
#'   (`trended.dispersion`). You can also pass
#'   [tidybulk::estimate_dispersion()]'s `dispersion_trended` (\eqn{s_0^2}),
#'   not `dispersion_shrinked`. Default `NULL` puts a zero offset on the
#'   shape submodel (see [estimate_gene()]).
#' @param dispersion_prior_log_sd Optional name of a log-dispersion prior SD
#'   column written by [estimate_dispersion_prior()] as
#'   `"dispersion_prior_log_sd"`. That function's `method` chooses the width.
#'   Default `NULL` uses a log-scale SD of 1.
#'
#' @details
#' Neither the offset nor the dispersion is computed here. Run
#' [estimate_dispersion_prior()] on the whole object before [brmDE()],
#' exactly as you compute the offset, so that both are ordinary columns of
#' the input, then pass those column names here. Gene-wise fitting is the
#' only thing this pipeline does. The dispersion columns are a prior on the
#' negative binomial shape (see [estimate_gene()]): they inform each gene
#' without fixing the posterior to the external estimate. Omitting
#' `dispersion_prior_log_mean` and/or `dispersion_prior_log_sd` is valid
#' (zero shape offset, log-scale SD of 1). Computing both, then passing
#' both names, is the preferred starting point.
#'
#' Prior constants derived from each gene are passed to Stan as data, so every
#' gene generates identical Stan code and cmdstanr compiles it at most once per
#' worker process rather than once per gene.
#' @param bundle Number of genes to fit per target. Default `10`. `1` gives
#'   one target per gene. A larger value fits that many genes sequentially
#'   inside one target, which is how you stop tens of thousands of genes from
#'   swamping an HPC scheduler with tiny jobs. See *Bundling* below.
#' @param target_output Name of the targets output.
#' @param ... Passed to [estimate_gene()] (e.g. `family`, `chains`,
#'   `draws_sampling`).
#'
#' @return The updated `tidytargets` pipeline.
#'
#' @section Where the fits are:
#' A `brmsfit` carries its draws and its data, so at 20,000 genes the fits are
#' far too large to be gathered into one table, and gathering them would undo
#' the point of having a store. [estimate()] therefore contributes the gene
#' column alone; what makes the result worth reading is [hypothesis()], whose
#' per-contrast statistics and convergence are small enough to hold for every
#' gene at once (see [hypothesis_gene()]). The fits stay in the store. Read
#' one with `targets::tar_read(brms_fit, branches = i, store = store)`, using
#' the same `store` you passed to [brmDE()], or
#' `tidytargets::tt_explore(pipeline, "brms_fit", index = i)`. With
#' `bundle = 1` the branch index is the gene's row in the result.
#'
#' @section What the pipeline remembers:
#' What the fits were asked for goes onto the pipeline with
#' [tidytargets::tt_metadata()]: the column names, and every
#' [estimate_gene()] argument the sampling depends on, defaults included, so
#' the object carries its own specification while the fits stay in the store.
#' The formulas are targets, so changing one invalidates the fits. Read the
#' rest with `tt_metadata(pipeline)`. It is how [plot_volcano()] knows the
#' `chains * draws_sampling` draws a `pH0` was counted from, and so why it
#' takes the pipeline rather than the table.
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
#' * A bundle is one store file holding `bundle` fits, and it is read back
#'   into the main process, so the size of a `brmsfit` is what limits how
#'   far you can push this.
#' * Invalidation is per bundle. One gene failing loses its bundle, and
#'   changing the gene set reshuffles bundle membership and refits everything.
#'   Leave `bundle` at `1` when you rely on incremental reruns.
#'
#' @export
estimate <- function(input_hpc,
                     formula_abundance,
                     formula_dispersion = ~1,
                     offset,
                     dispersion_prior_log_mean = NULL,
                     dispersion_prior_log_sd = NULL,
                     ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
estimate.default <- function(input_hpc,
                             formula_abundance,
                             formula_dispersion = ~1,
                             offset,
                             dispersion_prior_log_mean = NULL,
                             dispersion_prior_log_sd = NULL,
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
                            dispersion_prior_log_mean = NULL,
                            dispersion_prior_log_sd = NULL,
                            bundle = 10L,
                            target_output = "brms_fit",
                            ...) {
  offset <- check_offset_name(offset)
  bundle <- check_bundle(bundle)
  abundance <- tidytargets::tt_metadata(input_hpc)$abundance

  estimate_args <- c(
    list(
      abundance = abundance,
      offset = offset,
      dispersion_prior_log_mean = dispersion_prior_log_mean,
      dispersion_prior_log_sd = dispersion_prior_log_sd
    ),
    list(...)
  )

  # What the fits were asked for rides along on the pipeline, since the fits
  # themselves stay in the store: estimate_gene()'s arguments, with whatever
  # `...` overrode, so nothing is repeated here and nothing is missing when it
  # was left at a default. What any entry is for is up to whoever reads it;
  # plot_volcano() multiplies chains by draws_sampling for the resolution of a
  # counted probability.
  settings <- utils::modifyList(as.list(formals(estimate_gene)), list(...))
  input_hpc <- tidytargets::tt_metadata(
    input_hpc,
    offset = offset,
    dispersion_prior_log_mean = dispersion_prior_log_mean,
    dispersion_prior_log_sd = dispersion_prior_log_sd,
    chains = settings$chains,
    draws_warmup = settings$draws_warmup,
    draws_sampling = settings$draws_sampling
  )

  pipeline <- input_hpc |>
    tidytargets::tt_data(formula_abundance <- formula_text(formula_abundance)) |>
    tidytargets::tt_data(formula_dispersion <- formula_text(formula_dispersion)) |>
    tidytargets::tt_data(estimate_args) |>
    tidytargets::tt_data_list(
      gene_bundle <- bundle_gene_ids(
        tidytargets::tt_metadata(input_hpc)$gene_ids,
        bundle
      )
    ) |>
    tidytargets::tt_iterate(
      command = estimate_gene_from_se(
        se_input,
        gene_bundle,
        formula_abundance,
        formula_dispersion,
        estimate_args
      ),
      target_output = target_output
    )

  as_brmde_hpc(pipeline)
}

#' Hypothesis tests on a tidytargets pipeline
#'
#' Appends a [tidytargets::tt_iterate()] step that calls [hypothesis_gene()]
#' on each [estimate()] fit. This is an S3 method for
#' [brms::hypothesis()].
#'
#' Every contrast gets its posterior statistics and the convergence of its own
#' draws, so this is where a run of 20,000 genes becomes something you can
#' read: small enough to return in full, while the fits stay in the store.
#' Collecting unnests these into the result, one row per gene and contrast.
#' See [hypothesis_gene()] for the columns.
#'
#' Collecting the pipeline adds an `fdr` column across genes when the tests
#' took the directional route, which is the default. That is the one quantity
#' a gene-wise target cannot compute, since it is a property of the gene set.
#' The rate is the cumulative mean of `pH0` in ascending order.
#'
#' @param x A `tidytargets` pipeline from [brmDE()].
#' @param hypothesis Effects to test, passed to [hypothesis_gene()].
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

  hypothesis_args <- c(list(hypothesis = hypothesis), list(...))

  x |>
    tidytargets::tt_data(hypothesis_args) |>
    tidytargets::tt_iterate(
      command = hypothesis_gene_from_fit(brms_fit, hypothesis_args),
      target_output = target_output
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

  adjust_args <- c(list(nullify = nullify), list(...))

  input_hpc |>
    tidytargets::tt_data(adjust_args) |>
    tidytargets::tt_iterate(
      command = adjust_gene_from_fit(brms_fit, adjust_args),
      target_output = target_output
    ) |>
    as_brmde_hpc()
}

#' @export
#' @importFrom tidytargets tt_evaluate
tidytargets::tt_evaluate

#' @export
#' @importFrom brms hypothesis
brms::hypothesis

#' @exportS3Method tidytargets::tt_evaluate
tt_evaluate.brmDE_hpc <- function(tt_input) {
  store <- tt_input$initialisation$store
  script <- paste0(store, ".R")
  # Each factory already assigns `target_list <- ...`; a trailing
  # `target_list` is enough to return it, and is stripped first so print()
  # is idempotent.
  lines <- readLines(script)
  lines <- lines[!grepl("^\\s*target_list\\s*$", lines)]
  writeLines(c(lines, "target_list"), script)

  reporter <- tt_input$initialisation$verbosity
  make_args <- list(
    script = script,
    store = store,
    callr_function = tidytargets::tt_metadata(tt_input)$callr_function
  )
  if (!is.null(reporter) && !identical(reporter, "undefined")) {
    make_args$reporter <- reporter
  }
  do.call(targets::tar_make, make_args)

  collect_brmde_hpc(tt_input)
}

#' @rdname brmDE
#' @param x A `brmDE_hpc` pipeline.
#' @param ... Passed to [print()].
#' @export
print.brmDE_hpc <- function(x, ...) {
  x |>
    tt_evaluate() |>
    print(...)
}
