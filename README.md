# brmDE

[![R build status](https://github.com/stemangiola/brmDE/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/stemangiola/brmDE/actions)

Mixed-effect modelling for differential gene expression analyses, based on Bayesian regression through [brms](https://paul-buerkner.github.io/brms/).
Each gene is fit with an arbitrary formula and a precomputed offset; hypothesis tests and covariate adjustment follow, including as a [tidytargets](https://github.com/stemangiola/tidytargets) `targets` pipeline.

## Install

```r
# from this directory
devtools::install()

# CmdStan backend (same three-step setup as sccomp)
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev/", getOption("repos"))
)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
```

## Functions

### Single-gene functions

| Function | Role |
|----------|------|
| `estimate_gene()` | Fit one gene (ZINB mixed model by default) |
| `hypothesis_gene()` | Posterior hypothesis tests on that fit |
| `adjust_gene()` | Residual-plus-fitted adjustment, dropping nuisance covariates |

### High-performance computing functions

| Function | Role |
|----------|------|
| `brmDE()` | Start a [tidytargets](https://github.com/stemangiola/tidytargets) pipeline (`tt_initialise()` + `tt_data()`) |
| `estimate()` / `hypothesis()` / `adjust()` | `tt_iterate()` steps on that same graph |

### Plotting functions

| Function | Role |
|----------|------|
| `plot_boxplot()` | Boxplot of one gene across a factor, with a posterior predictive overlay |
| `plot_volcano()` | Volcano of a pipeline, with `pH0` below `1 / ndraws` jittered in a shaded band |

## Example for one gene

```r
library(brmDE)
data("airway", package = "airway")

se <- airway["ENSG00000120129", ]
se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
se$dex <- relevel(factor(se$dex), ref = "untrt")

fit <- estimate_gene(
  se,
  formula_abundance = ~ dex + (1 | cell),
  offset = "offset",
  abundance = "counts"
)

hypothesis_gene(fit, "dextrt")

plot_boxplot(fit, factor = "dex", feature = "ENSG00000120129")
plot_boxplot(fit, factor = "dex", feature = "ENSG00000120129",
          remove_unwanted_effects = TRUE)

adjust_gene(
  fit,
  nullify = "dex",
  re_formula = ~(1 | cell)
)
```

Naming an effect on its own tests whether it clears `test_above_log2FC` (default `1`, a doubling of expression), reporting the effect in `log2_fold_change` and the probability that the better supported side is wrong in `pH0`. Both ways of being wrong are in that one number: too small to matter, or pointing the other way. Because it is a probability rather than a test statistic, pipelines turn it into an FDR by averaging and add that as an `fdr` column.

Neither `estimate_gene()` nor `brmDE()` normalises. Compute the offset yourself on all genes (TMM, or anything else), store it as a `colData` column, and pass that name to `offset`. Dispersion is the same kind of whole-matrix quantity: run `estimate_dispersion_prior()` before `brmDE()` and pass the two column names to `estimate()`. Those columns become a **prior** on the negative binomial shape, not a plug-in: the gene-wise likelihood can still pull the posterior away from the external estimate.

The two models are given separately. `formula_abundance` is the mean model (the mixed formula you pass to `estimate()` / `estimate_gene()`). `estimate_dispersion_prior()` is a separate, edgeR step: pass the fixed-effect analogue yourself (`~ dex + cell` for `~ dex + (1 | cell)`). `formula_dispersion` (default `~1`) is the model for the negative binomial shape. Neither should carry an offset: the library size offset and the `log(1/dispersion_trended)` offset are appended for you, and the assembled formulas are printed as a message so you can check them:

```
Abundance model (offset added by brmDE): counts ~ dex + (1 | cell) + offset(offset)
Dispersion model (offset added by brmDE): shape ~ 1 + offset(log(1/dispersion_trended))
```

## Example for many genes as parallel pipeline, easily deployable on your HPC

To run many genes as **one** tidytargets pipeline. The pipeline fits genes and nothing else, so the two whole-matrix quantities it consumes — the offset and the dispersion — are prepared first:

```r
data("airway", package = "airway")

se <- airway
se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
se <- estimate_dispersion_prior(se, formula_abundance = ~ dex + cell)

se |>
  brmDE(features = c("ENSG00000120129")) |>
  estimate(
    ~ dex + (1 | cell),
    offset = "offset",
    dispersion_prior_log_mean = "dispersion_prior_log_mean",
    dispersion_prior_log_sd = "dispersion_prior_log_sd",
    family = brms::negbinomial()
  ) |>
  hypothesis("dextrt") |>
  adjust(nullify = "dex")
```

Printing that object runs the pipeline (`print` calls `tt_evaluate()`, as in tidytargets). Assigning it does not; an interactive session then says the pipeline is ready to be evaluated, rather than appearing to do nothing.

The result is one row per gene and contrast. `estimate()` contributes the gene column, and `hypothesis()` the statistics worth returning for all of them, unnested into the table: each contrast's posterior summary together with the convergence of that contrast's own draws (`rhat`, `ess_bulk`, `mcse`). The fits are not in the table — twenty thousand `brmsfit` objects will not fit in one — but they are in the targets store: `targets::tar_read(brms_fit, branches = i, store = store)` reads branch `i` (the gene's row, when `bundle = 1`).

## Plotting

`plot_volcano(pipeline)` draws the run as a volcano, evaluating the pipeline first if printing it has not already. It takes the pipeline rather than the table because `pH0` is counted from draws and so cannot resolve below `1 / ndraws`, and only the pipeline knows that number: `estimate()` records the sampling settings on it with `tt_metadata()`, while the fits themselves stay in the store. The genes that come back as exactly 0 are jittered across a shaded band one decade below the dashed resolution line, which reads as "smaller than these draws can measure" rather than as a probability.

By default `estimate()` fits 10 genes per target. At transcriptome scale that can swamp an HPC scheduler with tiny jobs, so `estimate(bundle = 100)` fits 100 genes per target instead; the output is unchanged, one row per gene. Set `bundle = 1` for one target per gene.

See the package vignette for the full walkthrough, including both dispersion-prior width methods, and `vignette("dispersion-priors")` for the derivation that connects edgeR's \(s_0^2\) and \(d_0\) to those priors:

```r
vignette("brmDE", package = "brmDE")
vignette("dispersion-priors", package = "brmDE")
```

## Dispersion priors

`estimate_dispersion_prior()` writes `dispersion_prior_log_mean` (\(\phi_{\mathrm{trend}}\), the across-gene trend) and `dispersion_prior_log_sd` (a width on \(\log\phi\)). Pass those names to `estimate()` / `estimate_gene()`. `method = "curvature"` (default) fills the SD from the Laplace approximation of edgeR's weighted smoothed shared log-likelihood; `method = "degrees_freedom"` fills it from the trigamma SD of the moderated df \(d_{\mathrm{residual}}+d_0\). Do not pass a tagwise / shrinked dispersion, which already includes this gene's counts. The shape intercept prior is always Student-t. Default `NULL` for the SD is a log-scale SD of 1.

| `method` | Prior | Notes |
| --- | --- | --- |
| `"curvature"` | `student_t` on the intercept of `shape ~ 1 + offset(log(1/dispersion_prior_log_mean))` | Laplace \(\sigma\) from the shared log-likelihood |
| `"degrees_freedom"` | `student_t` on the same intercept | trigamma \(\sigma\) from \(d_{\mathrm{eff}}\) |
| skip the writer, pass `NULL` | `student_t` on the same intercept | log-scale SD 1 |

See `?estimate_gene` and `vignette("dispersion-priors")`.
