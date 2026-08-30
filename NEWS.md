# brmDE NEWS

## brmDE 0.1.1

* The gene-wise HPC pipeline uses the tidytargets 0.0.7 grammar
  (`tt_data()`, `tt_data_list()`, `tt_iterate()`). Requires
  `tidytargets (>= 0.0.7)`.
* `brmDE()`, `estimate()`, `hypothesis()`, and `adjust()` usage is
  unchanged. Formulas and extra `estimate_gene()` arguments are
  pipeline targets, so changing them invalidates the fits.
* `estimate(bundle = )` always maps over gene bundles; `bundle = 1`
  is one gene per branch.
