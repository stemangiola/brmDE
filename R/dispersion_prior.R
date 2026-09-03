# Local Normal (Laplace) approximation of the weighted smoothed edgeR
# log-likelihood, as a prior SD on log(dispersion). See
# vignette("dispersion-priors").

#' Local quadratic SD of an unnormalised log-density
#'
#' Converts a discrete log-density on \eqn{\log(\phi)} into a Normal SD by
#' fitting \eqn{y = a + b z + c z^2} in a window around `x_mode`, with
#' \eqn{z = x - x_mode}. The curvature at the mode is \eqn{2c}; matching a
#' Normal \eqn{N(x_{\mathrm{mode}}, \sigma^2)} gives
#' \eqn{\sigma = \sqrt{-1/(2c)}}. That is the Laplace approximation used as
#' the default shape-prior width (`dispersion_prior_log_sd` with
#' `method = "curvature"`): see
#' `vignette("dispersion-priors")`.
#'
#' Returns `NA` when the window is too small, the quadratic coefficient is
#' missing, or the curve is not concave (`c >= 0`).
#'
#' @param x Numeric grid of \eqn{\log(\phi)} values, one per column of the
#'   shared log-likelihood.
#' @param y Numeric log-density values on the same grid, typically
#'   `prior.n * shared.loglik` for one gene.
#' @param x_mode Location at which the curvature is evaluated, typically
#'   `log(trended.dispersion)`.
#' @param n_local Number of grid points used in the local fit, centred on the
#'   point nearest `x_mode`. Default `5`.
#'
#' @return A single positive finite SD, or `NA_real_`.
#'
#' @noRd
dispersion_quadratic_log_sd <- function(x, y, x_mode, n_local = 10L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || !is.finite(x_mode)) {
    return(NA_real_)
  }
  i0 <- which.min(abs(x - x_mode))
  half <- floor((n_local - 1L) / 2L)
  idx <- seq.int(max(1L, i0 - half), min(length(x), i0 + half))
  if (length(idx) < 3L) {
    return(NA_real_)
  }
  z <- x[idx] - x_mode
  fit <- stats::lm(y[idx] ~ z + I(z^2))
  cf <- stats::coef(fit)
  if (length(cf) < 3L) {
    return(NA_real_)
  }
  c_hat <- unname(cf[[3L]])
  if (!is.finite(c_hat) || c_hat >= 0) {
    return(NA_real_)
  }
  sqrt(-1 / (2 * c_hat))
}

# SD of log(dispersion) implied by d_eff effective degrees of freedom.
# Var(log s^2) = trigamma(d/2) exactly for s^2 ~ sigma^2 chi^2_d / d;
# sqrt(2 / d) is the large-d approximation and is too small at low d.
# Non-finite or non-positive d_eff (edgeR's Inf / NA) become NA.
dispersion_log_sd_from_degrees_freedom <- function(d_eff) {
  d_eff <- as.numeric(d_eff)
  sd <- sqrt(trigamma(d_eff / 2))
  sd[!is.finite(d_eff) | d_eff <= 0 | !is.finite(sd) | sd <= 0] <- NA_real_
  sd
}

#' Per-gene Laplace SDs from an edgeR shared-likelihood grid
#'
#' Applies `dispersion_quadratic_log_sd()` to each gene: the log-density is
#' `prior_n[g] * shared_loglik[g, ]`, evaluated at `log(phi)` and centred
#' on `log(phi_trend[g])`.
#'
#' @param shared_loglik Numeric matrix, genes by dispersion-grid points, as
#'   returned by `edger_shared_loglik()` (`shared.loglik`).
#' @param phi Numeric vector of dispersion-grid values (the columns of
#'   `shared_loglik`), equal to `0.1 * 2^spline.pts` in edgeR.
#' @param prior_n Numeric vector of empirical-Bayes weights (`prior.n`), one
#'   per gene. Recycled to `nrow(shared_loglik)` if scalar.
#' @param phi_trend Numeric vector of trended dispersions, one per gene.
#' @param n_local Number of grid points in each local quadratic. Default `5`.
#'
#' @return Numeric vector of log-dispersion prior SDs, one per gene, with
#'   `NA` where the quadratic cannot be formed.
#'
#' @noRd
dispersion_quadratic_log_sd_over_genes <- function(shared_loglik,
                                               phi,
                                               prior_n,
                                               phi_trend,
                                               n_local = 10L) {
  x <- log(as.numeric(phi))
  n_gene <- nrow(shared_loglik)
  w <- rep_len(as.numeric(prior_n), n_gene)
  phi_trend <- as.numeric(phi_trend)
  vapply(seq_len(n_gene), function(g) {
    dispersion_quadratic_log_sd(
      x,
      w[[g]] * shared_loglik[g, ],
      log(phi_trend[[g]]),
      n_local = n_local
    )
  }, numeric(1))
}

#' Rebuild edgeR's abundance-smoothed shared log-likelihood
#'
#' `estimateDisp()` does not return `m0` / `shared.loglik`. This recomputes
#' the Cox--Reid adjusted profile-likelihood grid and the `WLEB()` smooth
#' that `estimateDisp()` uses, plus the empirical-Bayes weight `prior.n`
#' (`limma::squeezeVar()` of residual variances) and the trended dispersion.
#'
#' Genes with `rowSums(y) < min.row.sum` are left as `NA` in `shared.loglik`
#' and `prior.n`; their trended value is that of the lowest-abundance
#' selected gene, as in `estimateDisp()`.
#'
#' @param y Count matrix, genes by samples.
#' @param design Numeric design matrix for the fixed-effect mean model.
#' @param grid.length Number of dispersion-grid points. Default `61`.
#'   `estimateDisp()` uses `21` over `c(-10, 10)` (spacing 1 in `spline.pts`).
#'   The wider default `grid.range` would be coarser at 21 points, so this
#'   uses spacing 0.5.
#' @param grid.range Range of the grid in `spline.pts` units, where
#'   `phi = 0.1 * 2^spline.pts`. Default `c(-20, 10)`, wider on the small-`phi`
#'   side than edgeR's `c(-10, 10)`.
#' @param min.row.sum Minimum row sum to include a gene in the grid.
#'   Default `5`.
#' @param trend.method Smoothing method passed to [edgeR::WLEB()]. Default
#'   `"locfit"`.
#'
#' @return A list with
#'   \describe{
#'     \item{`phi`}{Dispersion-grid values, length `grid.length`.}
#'     \item{`shared.loglik`}{Matrix, genes by grid points; `NA` rows for
#'       genes below `min.row.sum`.}
#'     \item{`prior.n`}{Empirical-Bayes weight per gene; `NA` for unselected
#'       genes.}
#'     \item{`trended.dispersion`}{Trended `\phi` per gene.}
#'     \item{`effective_degrees_freedom`}{Moderated df
#'       \(d_{\mathrm{residual}} + d_0\). Unselected genes get the same
#'       scalar \(d_0\) as selected genes (tidybulk's `rep_len`); gene-wise
#'       \(d_0\) leaves them `NA`.}
#'   }
#'
#' @noRd
edger_shared_loglik <- function(y,
                                design,
                                grid.length = 61L,
                                grid.range = c(-20, 10),
                                min.row.sum = 5,
                                trend.method = "locfit") {
  y <- as.matrix(y)
  ntags <- nrow(y)
  nlibs <- ncol(y)
  design <- as.matrix(design)
  lib.size <- colSums(y)
  sel <- rowSums(y) >= min.row.sum
  sely <- y[sel, , drop = FALSE]
  offset <- matrix(log(lib.size), nrow = nrow(sely), ncol = nlibs, byrow = TRUE)

  # --- phi: candidate dispersions to get likelihoods for (vignette appendix 1) ---
  # edgeR's grid is in log2 units of phi / 0.1; convert to phi.
  spline.pts <- seq(from = grid.range[1], to = grid.range[2], length.out = grid.length)
  phi <- 0.1 * 2^spline.pts

  # --- gene-wise non-smoothed log-likelihood curves (vignette appendix 1) ---
  # l0[g, k] = Cox-Reid adjusted profile log-likelihood of gene g at phi[k].
  l0 <- matrix(0, nrow(sely), grid.length)
  last.beta <- NULL
  for (i in seq_len(grid.length)) {
    out <- edgeR::adjustedProfileLik(
      phi[[i]],
      y = sely,
      design = design,
      offset = offset,
      start = last.beta,
      get.coef = TRUE
    )
    l0[, i] <- out$apl
    last.beta <- out$beta
  }

  # --- shared.loglik: abundance-smoothed l0, edgeR's m0 (vignette appendix 2-3) ---
  # WLEB smooths each column of l0 against A_g, then reads the smooth back at
  # each gene's abundance: one shared curve ell_smooth,g(phi) per selected gene.

  # Abundance A_g, used to smooth those curves. The common dispersion is only
  # a plug-in for aveLogCPM, not a returned prior quantity.
  overall <- edgeR::maximizeInterpolant(spline.pts, matrix(colSums(l0), nrow = 1))
  AveLogCPM <- edgeR::aveLogCPM(
    y,
    lib.size = lib.size,
    dispersion = 0.1 * 2^overall
  )

  # Smoothing the log-likelihood curves
  out.1 <- edgeR::WLEB(
    theta = spline.pts,
    loglik = l0,
    covariate = AveLogCPM[sel],
    trend.method = trend.method,
    overall = FALSE,
    individual = FALSE,
    m0.out = TRUE
  )
  shared.loglik <- matrix(NA_real_, ntags, grid.length)
  shared.loglik[sel, ] <- out.1$shared.loglik

  # --- trended.dispersion: argmax of that shared curve (vignette appendix 4) ---
  # out.1$trend is spline.pts at the maximum; convert to phi. Unselected genes
  # take the trend of the lowest-abundance selected gene, as estimateDisp() does.
  trended.dispersion <- rep(0.1 * 2^out.1$trend[which.min(AveLogCPM[sel])], ntags)
  trended.dispersion[sel] <- 0.1 * 2^out.1$trend

  # --- prior.n: weight on the shared curve (vignette appendix 5) ---
  # squeezeVar's prior.df, divided by residual df of the design. Larger prior.n
  # sharpens w * ell_smooth; it does not move the mode.
  glmfit <- edgeR::glmFit(
    sely,
    design = design,
    offset = offset,
    dispersion = trended.dispersion[sel],
    prior.count = 0
  )
  df.residual <- glmfit$df.residual
  # Residual variance of the NB GLM at the trended dispersion: deviance / residual
  # df. squeezeVar() uses these as the gene-wise variances to shrink, and the
  # resulting prior.df becomes prior.n. Genes with no residual df are set to 0
  # so they do not contribute Inf.
  s2 <- glmfit$deviance / df.residual
  s2[df.residual == 0] <- 0
  s2 <- pmax(s2, 0)
  s2.fit <- limma::squeezeVar(s2, df = df.residual, covariate = AveLogCPM[sel])
  # prior.n = d0 / residual df of the design (NA for unselected genes; they
  # have no shared curve to weight). d_eff = residual df + d0, the moderated
  # df used by method = "degrees_freedom".
  # tidybulk writes the same quantity as dispersion_degrees_freedom.
  # When d0 is a scalar, fill every gene the same way tidybulk does.
  df.residual.design <- nlibs - ncol(design)
  df.prior <- s2.fit$df.prior
  prior.n <- rep(NA_real_, ntags)
  prior.n[sel] <- df.prior / df.residual.design
  effective_degrees_freedom <- rep(NA_real_, ntags)
  effective_degrees_freedom[sel] <- df.residual.design + df.prior
  if (length(df.prior) == 1L) {
    effective_degrees_freedom[!sel] <- df.residual.design + df.prior
  }

  # --- return the results ---
  list(
    phi = phi,
    shared.loglik = shared.loglik,
    prior.n = prior.n,
    trended.dispersion = trended.dispersion,
    effective_degrees_freedom = effective_degrees_freedom
  )
}

#' Log-scale SD and location of the edgeR smoothed-dispersion prior
#'
#' Writes two `rowData` columns from the shared log-likelihood:
#' `dispersion_prior_log_mean` is `trended.dispersion` (\eqn{\phi_{\mathrm{trend}}}),
#' and `dispersion_prior_log_sd` is a width on \eqn{\log(\phi)} from the local
#' Normal (Laplace) approximation of the weighted smoothed shared
#' log-likelihood. Pass both column names to [estimate_gene()] /
#' [estimate()] as `dispersion_prior_log_mean` and `dispersion_prior_log_sd`.
#' See `vignette("dispersion-priors")`.
#'
#' The shape intercept prior in [estimate_gene()] / [estimate()] is always
#' Student-t. That heavier-tailed prior is what lets the posterior leave a
#' misplaced edgeR mode under mixed-effect or zero-inflated models; this
#' function only supplies the location and a Laplace width.
#'
#' `estimateDisp()` does not return the shared log-likelihood, so this
#' rebuilds it (the same Cox--Reid profile likelihood and `WLEB()` smooth
#' that `estimateDisp()` uses) and fits a local quadratic around
#' `log(trended.dispersion)`.
#'
#' @param .data A `SummarizedExperiment` with counts in `abundance`.
#' @param formula_abundance Fixed-effect model for the mean, used to build
#'   the edgeR design. edgeR has no random effects, so write the analogue
#'   yourself (`~ dex + cell` for `~ dex + (1 | cell)`). A `|` in the
#'   formula is an error.
#' @param abundance Assay name holding integer counts. Default `"counts"`.
#' @param method How `log_sd_column` is filled. The documented value is
#'   `"curvature"` (default): local Normal approximation of the weighted
#'   smoothed shared log-likelihood.
#' @param log_sd_column Name of the `rowData` column to write for the
#'   log-scale SD. Default `"dispersion_prior_log_sd"`. Pass this name to
#'   [estimate_gene()] / [estimate()] as `dispersion_prior_log_sd`.
#' @param mean_column Name of the `rowData` column to write for the prior
#'   location, `trended.dispersion` (\eqn{\phi_{\mathrm{trend}}}). Default
#'   `"dispersion_prior_log_mean"`. Pass this name to [estimate_gene()] /
#'   [estimate()] as `dispersion_prior_log_mean`.
#' @param n_local Number of candidate dispersion values used in the local
#'   quadratic around each gene's trended dispersion. Default `10`.
#' @param grid.range Range of candidate dispersions in `spline.pts` units,
#'   where `phi = 0.1 * 2^spline.pts`. Default `c(-20, 10)`.
#' @param grid.length Number of candidate values in that range. Default `61`.
#'
#' @return `.data` with `mean_column` and `log_sd_column` added to `rowData`.
#'   The mean is the trended \eqn{\phi}. The SD is the Laplace width on
#'   \eqn{\log(\phi)} (`NA` where the quadratic cannot be formed).
#'
#' @docType methods
#' @rdname estimate_dispersion_prior-methods
#' @export
setGeneric(
  "estimate_dispersion_prior",
  function(.data,
           formula_abundance,
           abundance = "counts",
           method = c("curvature", "degrees_freedom"),
           log_sd_column = "dispersion_prior_log_sd",
           mean_column = "dispersion_prior_log_mean",
           n_local = 10L,
           grid.range = c(-20, 10),
           grid.length = 61L) {
    standardGeneric("estimate_dispersion_prior")
  }
)

.estimate_dispersion_prior_se <- function(.data,
                                           formula_abundance,
                                           abundance = "counts",
                                           method = c("curvature", "degrees_freedom"),
                                           log_sd_column = "dispersion_prior_log_sd",
                                           mean_column = "dispersion_prior_log_mean",
                                           n_local = 10L,
                                           grid.range = c(-20, 10),
                                           grid.length = 61L) {
  check_and_install_packages(c("edgeR", "limma"))
  method <- match.arg(method)
  log_sd_column <- check_column_name(log_sd_column, "log_sd_column")
  mean_column <- check_column_name(mean_column, "mean_column")
  if (identical(log_sd_column, mean_column)) {
    stop("`log_sd_column` and `mean_column` must name different columns.", call. = FALSE)
  }
  if (grepl("|", paste(deparse(formula_abundance), collapse = ""), fixed = TRUE)) {
    stop(
      "estimate_dispersion_prior() uses edgeR, which has no random effects. ",
      "Pass the fixed-effect analogue yourself, e.g. ~ dex + cell for ",
      "~ dex + (1 | cell).",
      call. = FALSE
    )
  }
  design <- stats::model.matrix(
    formula_abundance,
    data = as.data.frame(SummarizedExperiment::colData(.data))
  )
  shared <- edger_shared_loglik(
    SummarizedExperiment::assay(.data, abundance),
    design,
    grid.length = grid.length,
    grid.range = grid.range
  )
  SummarizedExperiment::rowData(.data)[[mean_column]] <-
    shared$trended.dispersion
  SummarizedExperiment::rowData(.data)[[log_sd_column]] <-
    if (identical(method, "curvature")) {
      dispersion_quadratic_log_sd_over_genes(
        shared$shared.loglik,
        shared$phi,
        shared$prior.n,
        shared$trended.dispersion,
        n_local = n_local
      )
    } else {
      dispersion_log_sd_from_degrees_freedom(shared$effective_degrees_freedom)
    }
  .data
}

#' @rdname estimate_dispersion_prior-methods
#' @export
setMethod(
  "estimate_dispersion_prior",
  "SummarizedExperiment",
  .estimate_dispersion_prior_se
)

#' @rdname estimate_dispersion_prior-methods
#' @export
setMethod(
  "estimate_dispersion_prior",
  "RangedSummarizedExperiment",
  .estimate_dispersion_prior_se
)
