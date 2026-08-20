# brmDE

Gene-wise Bayesian differential expression with [brms](https://paul-buerkner.github.io/brms/).
This package is the engine for a later `targets` pipeline: one gene in, a precomputed offset, an arbitrary formula, then estimate / hypothesis / adjust.

The archived [HPCell](https://github.com/MangiolaLaboratory/HPCell) tree and the immuneBodyMap `dynamic_tar_script.R` are **not** part of this package.

## Functions

| Function | Role |
|----------|------|
| `estimate_gene()` | Fit one gene (ZINB mixed model by default) |
| `hypothesis_gene()` | Posterior hypothesis tests on that fit |
| `adjust_gene()` | Residual-plus-fitted adjustment, dropping nuisance covariates |
| `estimate_dispersion()` | edgeR tagwise/trended dispersion on the full SE (`rowData`) |
| `brmDE()` | Start an [HPCell](https://github.com/MangiolaLaboratory/HPCell) targets pipeline (`initialise_hpc()` analogue) |
| `estimate()` / `hypothesis()` / `adjust()` | `hpc_iterate()` steps on that same graph |

```r
library(brmDE)
data("airway", package = "airway")

se <- airway["ENSG00000120129", ]
se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
se$dex <- relevel(factor(se$dex), ref = "untrt")

fit <- estimate_gene(
  se,
  formula = ~ dex + (1 | cell),
  offset = "offset",
  abundance = "counts"
)

hypothesis_gene(fit, "dextrt = 0")
hypothesis_gene(
  fit,
  hypothesis = "random_vs_rest",
  grouping = "cell"
)

adjust_gene(
  fit,
  nullify = "dex",
  re_formula = ~(1 | cell)
)
```

`estimate_gene()` does not calculate an offset or dispersion. Compute TMM first, store it as a column, and pass that name to `offset`. The pipeline writes edgeR dispersion to `rowData` before fitting. To run many genes as **one** HPCell pipeline:

```r
data("airway", package = "airway")

airway |>
  brmDE(features = c("ENSG00000120129")) |>
  estimate(
    ~ dex + (1 | cell),
    offset = "offset",
    dispersion = "dispersion",
    family = brms::negbinomial()
  ) |>
  hypothesis("dextrt = 0") |>
  adjust(nullify = "dex")
```

Printing that object runs the pipeline (`print` calls `evaluate_hpc()`, as in HPCell).

See the package vignette for the full walkthrough, including both shape-prior options:

```r
vignette("brmDE", package = "brmDE")
```

### Dispersion priors

`estimate_dispersion()` writes one `rowData` column per estimate: `dispersion` holds edgeR's gene-wise \(\phi_g\) and `dispersion_degrees_freedom` holds the effective degrees of freedom behind it, `d_eff = df.residual + prior.df`. Both names are arguments (`dispersion=`, `dispersion_degrees_freedom=`), and the same two names are passed on to `estimate()` / `estimate_gene()`. `estimate_gene()` turns that pair into a prior on the negative binomial shape, in one of two equivalent parameterisations chosen with `shape_prior`:

| `shape_prior` | Prior | Notes |
| --- | --- | --- |
| `"student_t"` (default) | `student_t` on the intercept of a `shape ~ 1 + offset(log(1/dispersion))` submodel | Heavier tails, tolerant of over-shrunk genes; leaves the submodel open to dispersion covariates |
| `"gamma"` | `gamma(d_eff/2, d_eff * dispersion/2)` on `shape` directly, no submodel | Conjugate to edgeR's scaled inverse chi-square hierarchy; lighter tails |

Both encode the same log-scale spread, `trigamma(d_eff / 2)`, because a chi-square *is* a gamma, and both put all their mass on a positive shape — the Student-t by exponentiating an unbounded parameter through brms' log link, the gamma by bounding it. They are not, however, reparameterisations of each other: the Student-t centres the *median* shape on edgeR's estimate while the gamma centres the *mean*, leaving their log-scale centres about `1/d_eff` apart. See `?estimate_dispersion` for the derivation and references.

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
