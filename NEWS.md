# brmDE NEWS

## brmDE 0.1.1

* `estimate_dispersion_prior()` takes a fixed-effect formula for the
  edgeR design. Pass the analogue yourself (`~ dex + cell` for
  `~ dex + (1 | cell)`); a `|` in the formula is an error.
* The gene-wise HPC pipeline uses the tidytargets 0.0.7 grammar
  (`tt_data()`, `tt_data_list()`, `tt_iterate()`). Requires
  `tidytargets (>= 0.0.7)`.
* `brmDE()`, `estimate()`, `hypothesis()`, and `adjust()` usage is
  unchanged. Formulas and extra `estimate_gene()` arguments are
  pipeline targets, so changing them invalidates the fits.
* `estimate(bundle = )` always maps over gene bundles; `bundle = 1`
  is one gene per branch.
