skip_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  if (!instantiate::stan_cmdstan_exists()) {
    skip("CmdStan not found.")
  }
}

# Back-compat alias used by older tests.
skip_if_no_cmdstan <- skip_cmdstan

airway_se <- function(n_genes = NULL, feature = "ENSG00000120129") {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  utils::data("airway", package = "airway", envir = environment())
  se <- get("airway", envir = environment())
  se$dex <- stats::relevel(factor(se$dex), ref = "untrt")
  se$cell <- droplevels(se$cell)

  if (!is.null(n_genes)) {
    totals <- rowSums(SummarizedExperiment::assay(se, "counts"))
    keep <- union(feature, names(sort(totals, decreasing = TRUE))[seq_len(n_genes)])
    se <- se[keep, ]
  }
  se
}

airway_one_gene <- function(feature = "ENSG00000120129") {
  se <- airway_se()
  if (!feature %in% rownames(se)) {
    skip(paste("Feature", feature, "not in airway"))
  }
  lib_size <- colSums(SummarizedExperiment::assay(se, "counts"))
  se$offset <- log(lib_size)
  se[feature, , drop = FALSE]
}

airway_one_gene_tbl <- function(feature = "ENSG00000120129") {
  brmDE:::as_gene_tibble(airway_one_gene(feature), abundance = "counts")
}

# The pipeline is gene-wise only, so the whole-matrix quantities it consumes -
# the offset and the edgeR dispersion prior - have to be on the object already.
airway_with_dispersion <- function(n_genes = 150, feature = "ENSG00000120129") {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  se <- airway_se(n_genes = n_genes, feature = feature)
  se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
  estimate_dispersion_prior(se, formula_abundance = ~ dex + cell)
}

airway_for_hpc <- airway_with_dispersion
