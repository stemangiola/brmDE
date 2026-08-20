# HPCell grammar wrappers. The gene-wise engines stay in estimate_gene(),
# hypothesis_gene(), and adjust_gene(); these methods only append them to an
# HPCell targets graph via hpc_single() / hpc_iterate().

as_brmde_hpc <- function(x) {
  class(x) <- unique(c("brmDE_hpc", "HPCell", class(x)))
  x
}

brmde_input_dir <- function(store) {
  paste0(store, "_input")
}

# targets hashes the arguments baked into a target's command, not the contents
# of any file those arguments point to. A fixed args path would therefore never
# change, and edited arguments would be answered with a stale cached result.
# Naming the file after a hash of its contents moves the change into the
# command. Elements are deparsed before hashing because objects such as brms
# families carry environments, whose raw hash differs between R sessions and
# would force a rebuild on every run.
write_args_rds <- function(input_dir, prefix, args) {
  fingerprint <- rlang::hash(
    vapply(args, function(x) paste(deparse(x), collapse = " "), character(1))
  )
  path <- file.path(
    input_dir,
    sprintf("%s_%s.rds", prefix, substr(fingerprint, 1, 16))
  )
  saveRDS(args, path)
  normalizePath(path, mustWork = FALSE)
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
      library(HPCell)
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
gene_ids_for_hpc <- function(se, features_rds) {
  ids <- rownames(se)
  features <- readRDS(features_rds)$features
  if (!is.null(features)) {
    ids <- intersect(ids, features)
  }
  if (length(ids) == 0L) {
    stop("No genes left to estimate after applying `features`.", call. = FALSE)
  }
  as.list(ids)
}

#' @keywords internal
#' @export
estimate_gene_from_se <- function(se,
                                  feature_id,
                                  formula_abundance,
                                  formula_dispersion,
                                  args_rds) {
  args <- readRDS(args_rds)
  named <- c("abundance", "offset", "dispersion", "dispersion_degrees_freedom")
  do.call(
    estimate_gene,
    c(
      list(
        data = se[unlist(feature_id), , drop = FALSE],
        formula_abundance = as_pipeline_formula(formula_abundance),
        formula_dispersion = as_pipeline_formula(formula_dispersion),
        offset = args$offset,
        abundance = args$abundance,
        dispersion = args$dispersion,
        dispersion_degrees_freedom = args$dispersion_degrees_freedom
      ),
      args[setdiff(names(args), named)]
    )
  )
}

#' @keywords internal
#' @export
hypothesis_gene_from_fit <- function(fit, args_rds) {
  do.call(hypothesis_gene, c(list(fit = fit), readRDS(args_rds)))
}

#' @keywords internal
#' @export
adjust_gene_from_fit <- function(fit, args_rds) {
  args <- readRDS(args_rds)
  if (is.character(args$re_formula) && !is.na(args$re_formula)) {
    args$re_formula <- stats::as.formula(
      args$re_formula,
      env = new.env(parent = globalenv())
    )
  }
  do.call(adjust_gene, c(list(fit = fit), args))
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

collect_brmde_hpc <- function(input_hpc) {
  store <- input_hpc$initialisation$store
  gene_id <- unlist(
    targets::tar_read_raw("gene_id", store = store),
    use.names = FALSE
  )
  out <- tibble::tibble(.feature = as.character(gene_id))

  if ("brms_fit" %in% names(input_hpc)) {
    out$brms_fit <- as_branch_list(
      targets::tar_read_raw("brms_fit", store = store)
    )
  }
  if ("hypothesis_tbl" %in% names(input_hpc)) {
    out$hypothesis <- as_branch_list(
      targets::tar_read_raw("hypothesis_tbl", store = store)
    )
  }
  if ("adjust_tbl" %in% names(input_hpc)) {
    out$adjust <- as_branch_list(
      targets::tar_read_raw("adjust_tbl", store = store)
    )
  }
  out
}

#' Start a gene-wise HPCell pipeline
#'
#' Analogue of HPCell's [HPCell::initialise_hpc()]: writes the targets
#' header (`target_list`) and the shared steps (load SE, gene ids). Later
#' calls to [estimate()], [hypothesis()], and [adjust()] append
#' `hpc_iterate()` branches onto the same graph. Printing the object (or
#' calling [evaluate_hpc()]) runs the pipeline, as in HPCell.
#'
#' The pipeline is gene-wise only. Anything estimated across the whole matrix,
#' namely the library size offset and the edgeR dispersion, belongs upstream of
#' [brmDE()] and arrives as ordinary `colData` / `rowData` columns.
#'
#' @param .data A `SummarizedExperiment` that already carries a library size
#'   offset in `colData` and edgeR dispersion in `rowData`. Neither is
#'   computed here: both are whole-matrix quantities, so derive them before
#'   the pipeline (`tidybulk::scale_abundance()` taking `log(1 / multiplier)`,
#'   and [estimate_dispersion()]) and name the columns when you call
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
#'   [evaluate_hpc()]. `NULL` (default) runs in the current session.
#' @param packages Character vector of packages loaded in workers.
#' @param update,garbage_collection,workspace_on_error Passed through
#'   to `tar_option_set()` / `tar_cue()` as in HPCell.
#' @param verbosity Reporter passed to [targets::tar_make()].
#'
#' @return An `HPCell` pipeline object (`brmDE_hpc`).
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#'
#' # Whole-matrix quantities are prepared before the gene-wise pipeline.
#' se <- airway
#' se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
#' se <- estimate_dispersion(se, formula_abundance = ~ dex + (1 | cell))
#'
#' se |>
#'   brmDE(features = c("ENSG00000120129")) |>
#'   estimate(
#'     ~ dex + (1 | cell),
#'     offset = "offset",
#'     dispersion = "dispersion",
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
#' #       slurm_memory_gigabytes_per_cpu = 5,
#' #       workers = 100,
#' #       tasks_max = 5
#' #     )
#' #   ) |>
#' #   estimate(~ dex + (1 | cell), offset = "offset", dispersion = "dispersion")
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
                        "HPCell",
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
  if (!requireNamespace("HPCell", quietly = TRUE)) {
    stop("Install HPCell to build the pipeline.", call. = FALSE)
  }

  dir.create(store, showWarnings = FALSE, recursive = TRUE)
  input_dir <- brmde_input_dir(store)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  script <- paste0(store, ".R")

  se_rds <- normalizePath(file.path(input_dir, "se.rds"), mustWork = FALSE)
  features_rds <- write_args_rds(input_dir, "features", list(features = features))
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
    HPCell::hpc_single("file_se", se_rds, format = "file") |>
    HPCell::hpc_single(
      target_output = "se_input",
      user_function = readRDS |> quote(),
      file = "file_se" |> HPCell::is_target()
    ) |>
    HPCell::hpc_single(
      target_output = "gene_id",
      user_function = gene_ids_for_hpc |> quote(),
      se = "se_input" |> HPCell::is_target(),
      features_rds = features_rds,
      iterate = "map"
    ) |>
    as_brmde_hpc()
}

#' Estimate genes on an HPCell pipeline
#'
#' Appends an [HPCell::hpc_iterate()] step that calls [estimate_gene()] once
#' per gene.
#'
#' @param input_hpc An `HPCell` / `brmDE_hpc` pipeline from [brmDE()].
#' @param formula_abundance Model for the mean, passed to [estimate_gene()].
#' @param formula_dispersion One-sided model for the negative binomial shape,
#'   passed to [estimate_gene()]. `~1` by default.
#' @param offset Required name of the precomputed offset column in `colData`
#'   (the same argument as [estimate_gene()]).
#' @param dispersion Required name of the `rowData` dispersion column, as
#'   written by [estimate_dispersion()].
#' @param dispersion_degrees_freedom Name of the effective degrees of freedom
#'   column written alongside it by [estimate_dispersion()].
#'
#' @details
#' Neither the offset nor the dispersion is computed here. Run
#' [estimate_dispersion()] on the whole object before [brmDE()], exactly as
#' you compute the offset, so that both are ordinary columns of the input.
#' Gene-wise fitting is the only thing this pipeline does.
#' @param target_output Name of the targets output.
#' @param ... Passed to [estimate_gene()] (e.g. `family`, `chains`, `iter`).
#'
#' @return The updated `HPCell` pipeline.
#'
#' @export
estimate <- function(input_hpc,
                     formula_abundance,
                     formula_dispersion = ~1,
                     offset,
                     dispersion,
                     ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
estimate.default <- function(input_hpc,
                             formula_abundance,
                             formula_dispersion = ~1,
                             offset,
                             dispersion,
                             ...) {
  stop(
    "estimate() expects a pipeline from brmDE(). ",
    "For one gene, use estimate_gene().",
    call. = FALSE
  )
}

#' @rdname estimate
#' @export
estimate.HPCell <- function(input_hpc,
                            formula_abundance,
                            formula_dispersion = ~1,
                            offset,
                            dispersion,
                            dispersion_degrees_freedom = "dispersion_degrees_freedom",
                            target_output = "brms_fit",
                            ...) {
  offset <- check_offset_name(offset)
  dispersion <- check_dispersion_name(dispersion)
  dispersion_degrees_freedom <-
    check_degrees_freedom_name(dispersion_degrees_freedom)
  abundance <- input_hpc$initialisation$abundance

  dots <- list(...)
  input_dir <- brmde_input_dir(input_hpc$initialisation$store)
  args_rds <- write_args_rds(
    input_dir,
    "estimate_args",
    c(
      list(
        abundance = abundance,
        offset = offset,
        dispersion = dispersion,
        dispersion_degrees_freedom = dispersion_degrees_freedom
      ),
      dots
    )
  )

  # The formulas are targets of their own so that editing one invalidates the
  # fits that depend on it.
  input_hpc |>
    HPCell::hpc_single(
      target_output = "formula_abundance_text",
      user_function = identity |> quote(),
      x = formula_text(formula_abundance)
    ) |>
    HPCell::hpc_single(
      target_output = "formula_dispersion_text",
      user_function = identity |> quote(),
      x = formula_text(formula_dispersion)
    ) |>
    HPCell::hpc_iterate(
      target_output = target_output,
      user_function = estimate_gene_from_se |> quote(),
      se = "se_input" |> HPCell::is_target(),
      feature_id = "gene_id" |> HPCell::is_target(),
      formula_abundance = "formula_abundance_text" |> HPCell::is_target(),
      formula_dispersion = "formula_dispersion_text" |> HPCell::is_target(),
      args_rds = args_rds
    ) |>
    as_brmde_hpc()
}

#' Hypothesis tests on an HPCell pipeline
#'
#' Appends an [HPCell::hpc_iterate()] step that calls [hypothesis_gene()]
#' on each [estimate()] fit. This is an S3 method for
#' [brms::hypothesis()].
#'
#' @param x An `HPCell` pipeline from [brmDE()].
#' @param hypothesis Hypothesis strings passed to [hypothesis_gene()].
#' @param target_input Name of the upstream fit target.
#' @param target_output Name of the targets output.
#' @param ... Passed to [hypothesis_gene()].
#'
#' @return The updated `HPCell` pipeline.
#'
#' @exportS3Method brms::hypothesis
hypothesis.HPCell <- function(x,
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

  args_rds <- write_args_rds(
    brmde_input_dir(x$initialisation$store),
    "hypothesis_args",
    c(list(hypothesis = hypothesis), list(...))
  )

  x |>
    HPCell::hpc_iterate(
      target_output = target_output,
      user_function = hypothesis_gene_from_fit |> quote(),
      fit = target_input |> HPCell::is_target(),
      args_rds = args_rds
    ) |>
    as_brmde_hpc()
}

#' Adjust genes on an HPCell pipeline
#'
#' Appends an [HPCell::hpc_iterate()] step that calls [adjust_gene()] on
#' each [estimate()] fit.
#'
#' @param input_hpc An `HPCell` pipeline, or a `brmsfit`.
#' @param nullify Covariates to nullify, passed to [adjust_gene()].
#' @param target_input Name of the upstream fit target.
#' @param target_output Name of the targets output.
#' @param ... Passed to [adjust_gene()].
#'
#' @return The updated `HPCell` pipeline, or an adjustment tibble.
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
adjust.HPCell <- function(input_hpc,
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

  dots <- list(...)
  if (!is.null(dots$re_formula) &&
      !identical(dots$re_formula, NA) &&
      inherits(dots$re_formula, "formula")) {
    dots$re_formula <- formula_text(dots$re_formula)
  }

  args_rds <- write_args_rds(
    brmde_input_dir(input_hpc$initialisation$store),
    "adjust_args",
    c(list(nullify = nullify), dots)
  )

  input_hpc |>
    HPCell::hpc_iterate(
      target_output = target_output,
      user_function = adjust_gene_from_fit |> quote(),
      fit = target_input |> HPCell::is_target(),
      args_rds = args_rds
    ) |>
    as_brmde_hpc()
}

#' @export
#' @importFrom HPCell evaluate_hpc
HPCell::evaluate_hpc

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

#' @exportS3Method HPCell::evaluate_hpc
evaluate_hpc.brmDE_hpc <- function(input_hpc) {
  store <- input_hpc$initialisation$store
  script <- paste0(store, ".R")
  localize_target_append(script)
  cat("target_list\n", file = script, append = TRUE)

  reporter <- input_hpc$initialisation$verbosity
  make_args <- list(
    script = script,
    store = store,
    callr_function = input_hpc$initialisation$callr_function
  )
  if (!is.null(reporter) && !identical(reporter, "undefined")) {
    make_args$reporter <- reporter
  }
  do.call(targets::tar_make, make_args)

  collect_brmde_hpc(input_hpc)
}

#' @rdname brmDE
#' @export
print.brmDE_hpc <- function(x, ...) {
  x |>
    evaluate_hpc() |>
    print(...)
}
