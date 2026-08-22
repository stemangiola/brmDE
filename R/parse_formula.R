
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
