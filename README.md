# brmDE

[![R build status](https://github.com/MangiolaLaboratory/brmDE/workflows/rworkflows/badge.svg)](https://github.com/MangiolaLaboratory/brmDE/actions)

Mixed-effect modelling for differential gene expression analyses, based on Bayesian regression through [brms](https://paul-buerkner.github.io/brms/).
Each gene is fit with an arbitrary formula and a precomputed offset; hypothesis tests and covariate adjustment follow, including as a [tidytargets](https://github.com/stemangiola/tidytargets) `targets` pipeline.

The archived [HPCell](https://github.com/MangiolaLaboratory/HPCell) tree and the immuneBodyMap `dynamic_tar_script.R` are **not** part of this package.

## Functions

| Function | Role |
|----------|------|
| `estimate_gene()` | Fit one gene (ZINB mixed model by default) |
| `hypothesis_gene()` | Posterior hypothesis tests on that fit |
| `false_discovery_rate()` | Turn per-gene null probabilities into an FDR across genes |
| `adjust_gene()` | Residual-plus-fitted adjustment, dropping nuisance covariates |
| `plot_gene()` | Boxplot of one gene across a factor, with a posterior predictive overlay |
| `plot_volcano()` | Volcano of a pipeline, with `pH0` below `1 / ndraws` jittered in a shaded band |
| `brmDE()` | Start a [tidytargets](https://github.com/stemangiola/tidytargets) pipeline (`tt_initialise()` analogue) |
| `estimate()` / `hypothesis()` / `adjust()` | `tt_iterate()` steps on that same graph |

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

plot_gene(fit, factor = "dex", feature = "ENSG00000120129")
plot_gene(fit, factor = "dex", feature = "ENSG00000120129",
          remove_unwanted_effects = TRUE)

adjust_gene(
  fit,
  nullify = "dex",
  re_formula = ~(1 | cell)
)
```

Naming an effect on its own tests whether it clears `test_above_log2FC` (default `1`, a doubling of expression), reporting the effect in `log2_fold_change` and the probability that the better supported side is wrong in `pH0`. Both ways of being wrong are in that one number: too small to matter, or pointing the other way. Because it is a probability rather than a test statistic, `false_discovery_rate()` turns it into an FDR by averaging, and pipelines add that as an `fdr` column.

Neither `estimate_gene()` nor `brmDE()` normalises. Compute the offset yourself on all genes (TMM, or anything else), store it as a `colData` column, and pass that name to `offset`. Dispersion is the same kind of whole-matrix quantity: run `tidybulk::estimate_dispersion()` before `brmDE()` and pass the two column names to `estimate()`. Those columns become a **prior** on the negative binomial shape, not a plug-in: the gene-wise likelihood can still pull the posterior away from the external estimate.

The two models are given separately. `formula_abundance` is the mean model. For `tidybulk::estimate_dispersion()` write the fixed-effect analogue of that model (`~ dex + cell` for `~ dex + (1 | cell)`): edgeR has no random effects. `formula_dispersion` (default `~1`) is the model for the negative binomial shape. Neither should carry an offset: the library size offset and the `log(1/dispersion)` offset are appended for you, and the assembled formulas are printed as a message so you can check them:

```
Abundance model (offset added by brmDE): counts ~ dex + (1 | cell) + offset(offset)
Dispersion model (offset added by brmDE): shape ~ 1 + offset(log(1/dispersion))
```

To run many genes as **one** tidytargets pipeline. The pipeline fits genes and nothing else, so the two whole-matrix quantities it consumes — the offset and the dispersion — are prepared first:

```r
data("airway", package = "airway")

se <- airway
se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
se <- tidybulk::estimate_dispersion(se, formula_abundance = ~ dex + cell)

se |>
  brmDE(features = c("ENSG00000120129")) |>
  estimate(
    ~ dex + (1 | cell),
    offset = "offset",
    dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    family = brms::negbinomial()
  ) |>
  hypothesis("dextrt") |>
  adjust(nullify = "dex")
```

Printing that object runs the pipeline (`print` calls `tt_evaluate()`, as in tidytargets).

The result is one row per gene and contrast. `estimate()` contributes the gene column, and `hypothesis()` the statistics worth returning for all of them, unnested into the table: each contrast's posterior summary together with the convergence of that contrast's own draws (`rhat`, `ess_bulk`, `mcse`). The fits are not in the table — twenty thousand `brmsfit` objects will not fit in one — but they are in the targets store: `targets::tar_read(brms_fit, branches = i, store = store)` reads branch `i` (the gene's row, when `bundle = 1`).

`plot_volcano(pipeline)` draws the run as a volcano, evaluating the pipeline first if printing it has not already. It takes the pipeline rather than the table because `pH0` is counted from draws and so cannot resolve below `1 / ndraws`, and only the pipeline knows that number: `estimate()` records the sampling settings on it with `tt_metadata()`, while the fits themselves stay in the store. The genes that come back as exactly 0 are jittered across a shaded band one decade below the dashed resolution line, which reads as "smaller than these draws can measure" rather than as a probability.

By default `estimate()` fits 10 genes per target. At transcriptome scale that can swamp an HPC scheduler with tiny jobs, so `estimate(bundle = 100)` fits 100 genes per target instead; the output is unchanged, one row per gene. Set `bundle = 1` for one target per gene.

See the package vignette for the full walkthrough, including both shape-prior options:

```r
vignette("brmDE", package = "brmDE")
```

### Dispersion priors

`tidybulk::estimate_dispersion()` writes one `rowData` column per estimate, named by `dispersion_column` and `dispersion_degrees_freedom_column` (defaults `"dispersion"` for gene-wise \(\phi_g\) and `"dispersion_degrees_freedom"` for \(d_{\mathrm{eff}} = \mathrm{df.residual} + \mathrm{prior.df}\)). Pass those same names to `estimate()` / `estimate_gene()` as `dispersion=` and `dispersion_degrees_freedom=`. Those two default to `NULL` (`offset(0)` on the shape submodel, Student-t log-scale SD of 1); computing the columns and passing both names is the preferred starting point. `estimate_gene()` turns that pair into a prior on the negative binomial shape — a location from \(\phi_g\) and a tightness from \(d_{\mathrm{eff}}\) — so data-driven evidence for that gene can still diverge from the external estimate. Choose the form with `shape_prior`:

| `shape_prior` | Prior | Notes |
| --- | --- | --- |
| `"student_t"` (default) | `student_t` on the intercept of a `shape ~ 1 + offset(log(1/dispersion))` submodel (`offset(0)` if `dispersion` is omitted) | Heavier tails, tolerant of over-shrunk genes; leaves the submodel open to dispersion covariates |
| `"gamma"` | `gamma(d_eff/2, d_eff * dispersion/2)` on `shape` directly, no submodel | Conjugate to edgeR's scaled inverse chi-square hierarchy; lighter tails |

Both encode the same log-scale spread, `trigamma(d_eff / 2)`, because a chi-square *is* a gamma, and both put all their mass on a positive shape — the Student-t by exponentiating an unbounded parameter through brms' log link, the gamma by bounding it. They are not, however, reparameterisations of each other: the Student-t centres the *median* shape on edgeR's estimate while the gamma centres the *mean*, leaving their log-scale centres about `1/d_eff` apart. See `?estimate_gene` for the two parameterisations.

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

`estimate_gene()` defaults to `backend = "cmdstanr"`. Vignette fitting chunks run only when `instantiate::stan_cmdstan_exists()` is true.
