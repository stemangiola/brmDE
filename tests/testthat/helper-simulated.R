# airway has only 8 samples, so `~ treatment * donor` exhausts the design.
# This synthetic negative binomial dataset is wide enough to keep residual
# degrees of freedom there, and records the dispersion it was generated with
# so estimates can be checked against a known truth.
simulated_se <- function(n_gene = 200,
                         n_donor = 6,
                         n_replicate = 2,
                         seed = 1) {
  skip_if_not_installed("SummarizedExperiment")

  col_data <- expand.grid(
    replicate = seq_len(n_replicate),
    treatment = c("untrt", "trt"),
    donor = paste0("D", seq_len(n_donor)),
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )
  n_sample <- nrow(col_data)

  sim <- withr::with_seed(seed, {
    phi <- stats::rgamma(n_gene, shape = 2, rate = 20)
    log_mu <-
      matrix(stats::rnorm(n_gene, mean = 6, sd = 1.5), n_gene, n_sample) +
      outer(
        stats::rnorm(n_gene, mean = 0, sd = 0.5),
        as.numeric(col_data$treatment == "trt")
      ) +
      matrix(
        stats::rnorm(n_donor, mean = 0, sd = 0.3)[factor(col_data$donor)],
        n_gene, n_sample,
        byrow = TRUE
      )
    list(
      phi = phi,
      # size = 1/phi recycles down the rows, matching one dispersion per gene.
      counts = matrix(
        stats::rnbinom(n_gene * n_sample, mu = exp(log_mu), size = 1 / phi),
        nrow = n_gene
      )
    )
  })

  counts <- sim$counts
  rownames(counts) <- sprintf("G%03d", seq_len(n_gene))
  colnames(counts) <- sprintf("S%02d", seq_len(n_sample))

  col_data$treatment <-
    stats::relevel(factor(col_data$treatment), ref = "untrt")
  col_data$donor <- factor(col_data$donor)
  col_data$offset <- log(colSums(counts))
  rownames(col_data) <- colnames(counts)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data,
    rowData = data.frame(
      true_dispersion = sim$phi,
      row.names = rownames(counts)
    )
  )
}
