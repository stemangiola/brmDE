#' Add edgeR dispersion estimates to rowData
#'
#' Port of HPCell's `se_estimate_dispersion()`
#' (<https://github.com/MangiolaLaboratory/HPCell/blob/ed763c6f53cfedf63d3a3353ec5ffcfd09bfe0e6/R/differential_expression.R#L36>):
#' tagwise dispersion when there are fewer than 1000 samples, otherwise
#' trended dispersion on a random subset of up to 2000 samples. Random-effect
#' terms in `formula_abundance` are turned into the fixed effects they imply
#' before the design matrix is built. [estimate_gene()] does not compute this;
#' [estimate()] calls this helper once on the full object before iterating
#' genes.
#'
#' @param se A `SummarizedExperiment`.
#' @param formula_abundance Model for the mean, used to build the edgeR
#'   design; the same formula you pass to [estimate_gene()]. edgeR has no
#'   random effects, so `(1 | donor)` enters the design as `donor` and
#'   `(1 + treatment | donor)` as `donor + treatment:donor`.
#' @param abundance Assay name (default `"counts"`).
#' @param dispersion Name of the `rowData` column for tagwise (or trended)
#'   dispersion \eqn{\phi_g} (default `"dispersion"`).
#' @param dispersion_degrees_freedom Name of the `rowData` column for the
#'   effective degrees of freedom \eqn{d_{eff}} behind that estimate (default
#'   `"dispersion_degrees_freedom"`).
#'
#' @details
#' Two columns are written, named by `dispersion` and
#' `dispersion_degrees_freedom`. The first holds \eqn{\phi_g}; the second holds
#' the effective degrees of freedom behind that estimate, which
#' [estimate_gene()] turns into a prior on the negative binomial shape.
#'
#' The effective degrees of freedom are
#' \deqn{d_{eff} = (n - \mathrm{ncol}(design)) + \mathrm{prior.df}}
#' The first term is the residual degrees of freedom of the fixed-effect
#' design; `estimateDisp()` does not return it, so it is computed here. The
#' second is the prior degrees of freedom that `estimateDisp()` does return,
#' quantifying how far each gene-wise estimate is shrunk toward the
#' mean-dispersion trend by empirical Bayes. `estimateTrendedDisp()` (the
#' branch used above 1000 samples) performs no such shrinkage and returns no
#' `prior.df`, so there \eqn{d_{eff}} is the residual degrees of freedom alone.
#'
#' Because the grouping factors enter as fixed effects, a factor with close to
#' one level per sample will exhaust the design and leave no residual degrees
#' of freedom. That is an error rather than a column of `NA`s: pass a simpler
#' `formula_abundance` for this step if it happens.
#'
#' The link from degrees of freedom to a standard deviation runs through the
#' chi-square distribution of the estimator. For \eqn{s^2 \sim \sigma^2
#' \chi^2_d / d}, the log of the estimate has variance
#' \deqn{\mathrm{Var}(\log s^2) = \psi'(d/2)}
#' with \eqn{\psi'} the trigamma function, so
#' \eqn{\mathrm{SD}(\log \hat\phi_g) \approx \sqrt{\psi'(d_{eff}/2)}}. The
#' familiar \eqn{\sqrt{2/d}} is the large-\eqn{d} approximation to this and is
#' 12% too small at \eqn{d = 4}, so the trigamma form is used directly. This
#' is the same empirical-Bayes variance model that limma and edgeR use to
#' moderate gene-wise variances; it is exact for a scaled chi-square and is
#' applied to the negative binomial dispersion by analogy.
#'
#' The same \eqn{d_{eff}} supports a second, conjugate parameterisation, which
#' [estimate_gene()] offers as `shape_prior = "gamma"`. Inverting the scaled
#' inverse chi-square gives a gamma, so the precision \eqn{1/\phi_g} -- the
#' quantity brms calls `shape` -- has prior
#' \deqn{\mathrm{Gamma}(d_{eff}/2,\ \mathrm{rate} = d_{eff}\phi_g/2)}
#' with mean \eqn{1/\phi_g}. Because \eqn{\mathrm{Var}(\log X) = \psi'(a)} for
#' a gamma of shape \eqn{a}, this reproduces \eqn{\psi'(d_{eff}/2)} exactly:
#' the two routes are the same calculation, since \eqn{\chi^2_d} is itself
#' \eqn{\mathrm{Gamma}(d/2, \mathrm{scale} = 2)}. Its coefficient of variation
#' is \eqn{\sqrt{2/d_{eff}}}, which is where that familiar approximation comes
#' from.
#'
#' The two forms are not reparameterisations of one another, and differ in two
#' ways worth knowing.
#'
#' They centre different summaries of the shape on edgeR's estimate. brms
#' gives the shape submodel a log link, so a Student-t on the intercept is
#' symmetric in \eqn{\log(\mathrm{shape})} and places the *median* of the
#' shape at \eqn{1/\phi_g}. The gamma places its *mean* there, and a gamma's
#' median lies below its mean. On the log scale the centres differ by
#' \deqn{\psi(d_{eff}/2) - \log(d_{eff}/2) \approx -1/d_{eff}}
#' with \eqn{\psi} the digamma function: about -0.11, or 0.22 prior standard
#' deviations, at \eqn{d_{eff} = 9.8}, and vanishing as samples accumulate.
#' Neither is wrong. Mean-centring is what the conjugate hierarchy dictates,
#' since a scaled inverse chi-square on \eqn{\phi_g} implies
#' \eqn{E[1/\phi_g] = 1/\phi_g}; median-centring is what an additive offset on
#' a log link implies.
#'
#' They also differ in tail weight, which is why the Student-t is the default.
#' Exponentiating a Student-t leaves a prior on the shape with polynomial
#' tails and no finite mean, proper but very permissive, so a gene whose true
#' dispersion is far from the shrunken edgeR estimate can still escape. The
#' gamma decays exponentially and holds such a gene closer to the trend. The
#' gamma is also left-skewed in \eqn{\log(\mathrm{shape})} (skewness
#' \eqn{\psi''(a) / \psi'(a)^{3/2}}, about -0.47 here) where the Student-t is
#' symmetric.
#'
#' Note that this standard deviation describes how precisely edgeR estimated
#' its own \eqn{\phi_g}. It does not account for the fact that the design here
#' fits the grouping factors as fixed rather than partially pooled, nor for a
#' zero-inflated likelihood in which the `zi` component absorbs part of the
#' overdispersion. Treating them as fixed is the conservative choice: it
#' spends the full degrees of freedom that shrinkage would have given back,
#' so \eqn{d_{eff}} understates rather than overstates the information behind
#' \eqn{\phi_g}.
#'
#' @references
#' Smyth GK (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. *Statistical
#' Applications in Genetics and Molecular Biology* 3(1).
#' [PDF](https://gksmyth.github.io/pubs/ebayes.pdf) — derives the
#' \eqn{\mathrm{Var}(\log s^2) = \psi'(d/2)} result and the prior degrees of
#' freedom used to moderate it.
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation.
#' *Nucleic Acids Research* 40(10):4288-4297.
#' [PMC3378882](https://pmc.ncbi.nlm.nih.gov/articles/PMC3378882/) —
#' Cox-Reid adjusted profile likelihood dispersion conditional on a design.
#'
#' Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package
#' for differential expression analysis of digital gene expression data.
#' *Bioinformatics* 26(1):139-140.
#' [PMC2796818](https://pmc.ncbi.nlm.nih.gov/articles/PMC2796818/)
#'
#' Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK (2016). Robust
#' hyperparameter estimation protects against hypervariable genes and improves
#' power to detect differential expression. *Annals of Applied Statistics*
#' 10(2):946-963. \doi{10.1214/16-AOAS920} — `robust = TRUE`, which makes
#' `prior.df` gene-specific and possibly infinite.
#'
#' [`estimateDisp()` reference manual](https://rdrr.io/bioc/edgeR/man/estimateDisp.html),
#' [edgeR on Bioconductor](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
#'
#' @return `se` with `rowData(se)[[dispersion]]` and
#'   `rowData(se)[[dispersion_degrees_freedom]]` filled.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- estimate_dispersion(airway[1:150, ], ~ dex + (1 | cell))
#' }
#'
#' @export
estimate_dispersion <- function(se,
                                formula_abundance,
                                abundance = "counts",
                                dispersion = "dispersion",
                                dispersion_degrees_freedom = "dispersion_degrees_freedom") {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (missing(formula_abundance) || is.null(formula_abundance)) {
    stop(
      "`formula_abundance` is required to build the edgeR design.",
      call. = FALSE
    )
  }
  dispersion <- check_dispersion_name(dispersion)
  dispersion_degrees_freedom <-
    check_degrees_freedom_name(dispersion_degrees_freedom)
  if (identical(dispersion, dispersion_degrees_freedom)) {
    stop(
      "`dispersion` and `dispersion_degrees_freedom` must name different columns.",
      call. = FALSE
    )
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Install edgeR to calculate dispersion.", call. = FALSE)
  }
  if (!abundance %in% SummarizedExperiment::assayNames(se)) {
    stop(
      "Assay '", abundance, "' not found. Available: ",
      paste(SummarizedExperiment::assayNames(se), collapse = ", "),
      call. = FALSE
    )
  }

  n_sample <- ncol(se)
  if (n_sample == 0L) {
    stop("`se` has no samples, so there is nothing to estimate.", call. = FALSE)
  }

  design_formula <- fixed_effects_formula(formula_abundance)

  if (n_sample < 1000L) {
    design <- dispersion_design(se, design_formula)
    check_residual_df(design, n_sample)
    counts <- SummarizedExperiment::assay(se, abundance)
    fit <- edgeR::estimateDisp(counts, design = design)
    disp <- fit$tagwise.dispersion
    d_eff <- (n_sample - ncol(design)) + fit$prior.df
  } else {
    sampled <- sample(seq_len(n_sample), size = min(n_sample, 2000L))
    se_sub <- se[, sampled, drop = FALSE]
    design <- dispersion_design(se_sub, design_formula)
    check_residual_df(design, ncol(se_sub))
    counts <- SummarizedExperiment::assay(se_sub, abundance)
    disp <- edgeR::estimateTrendedDisp(
      counts,
      design = design,
      subset = 1000,
      rowsum.filter = 10
    )
    d_eff <- ncol(se_sub) - ncol(design)
  }

  SummarizedExperiment::rowData(se)[[dispersion]] <- as.numeric(disp)
  SummarizedExperiment::rowData(se)[[dispersion_degrees_freedom]] <-
    rep_len(as.numeric(d_eff), nrow(se))
  se
}

# edgeR has no random effects, so `(f | group)` terms are expanded into the
# fixed effects they imply rather than discarded. Keeping the grouping factor
# in the design costs the degrees of freedom it really consumes, which is what
# d_eff is meant to report.
fixed_effects_formula <- function(formula) {
  if (inherits(formula, "brmsformula")) {
    formula <- formula$formula
  }
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  if (length(formula) == 2L) {
    formula <- stats::as.formula(
      paste(".response ~", deparse1(formula[[2]])),
      env = environment(formula)
    )
  }
  labels <- attr(stats::terms(formula), "term.labels")
  is_random <- grepl("|", labels, fixed = TRUE)
  fixed <- unique(c(
    labels[!is_random],
    unlist(lapply(labels[is_random], expand_random_term), use.names = FALSE)
  ))
  if (length(fixed) == 0L) {
    return(~1)
  }
  stats::as.formula(paste("~", paste(fixed, collapse = " + ")))
}

# `(1 | group)` becomes `1:group`, which is just `group`; `(1 + f1:f2 | group)`
# becomes `group + f1:f2:group`. Both sides go through terms() first, so `||`
# behaves like `|`, `(0 + f | group)` drops the intercept, and a nested
# `(1 | a/b)` expands to `a + a:b`. brms' `(f | ID | group)`, which correlates
# group-level effects across distributional parameters, contributes the same
# fixed effects as `(f | group)`: the ID only ties posteriors together.
expand_random_term <- function(label) {
  sides <- trimws(strsplit(label, "|", fixed = TRUE)[[1]])
  sides <- sides[nzchar(sides)]
  if (length(sides) < 2L) {
    stop(
      "Could not read the random-effect term '", label,
      "'. Expected the form (terms | group).",
      call. = FALSE
    )
  }
  effect <- formula_side_terms(sides[[1]])
  group <- formula_side_terms(sides[[length(sides)]])$labels
  out <- if (effect$intercept) group else character(0)
  if (length(effect$labels) > 0L && length(group) > 0L) {
    out <- c(out, as.vector(outer(effect$labels, group, paste, sep = ":")))
  }
  out
}

formula_side_terms <- function(text) {
  terms <- stats::terms(stats::as.formula(paste("~", text)))
  list(
    labels = attr(terms, "term.labels"),
    intercept = attr(terms, "intercept") == 1L
  )
}

# Expanding `(1 | group)` into a fixed `group` can exhaust the design when the
# grouping factor has close to one level per sample. edgeR would return NA for
# every gene, so stop instead of handing back a column of NAs.
check_residual_df <- function(design, n_sample) {
  if (n_sample - ncol(design) > 0L) {
    return(invisible(NULL))
  }
  stop(
    "The dispersion design has ", ncol(design), " coefficients for ",
    n_sample, " samples, leaving no residual degrees of freedom. ",
    "Random-effect terms are fitted as fixed effects here, so a grouping ",
    "factor with nearly one level per sample will exhaust the design. ",
    "Pass a simpler `formula_abundance` to estimate_dispersion().",
    call. = FALSE
  )
}

dispersion_design <- function(se, formula) {
  stats::model.matrix(
    formula,
    data = droplevels(as.data.frame(SummarizedExperiment::colData(se)))
  )
}
