# Internal helpers for gene-wise brms DE.

collapse_repeated_underscores <- function(x) {
  stringr::str_replace_all(x, "_+", "_")
}

check_column_name <- function(x, arg) {
  if (
    missing(x) ||
      is.null(x) ||
      !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
  ) {
    stop(
      sprintf("`%s` must be a single column name.", arg),
      call. = FALSE
    )
  }
  x
}

check_offset_name <- function(offset) {
  check_column_name(offset, "offset")
}

check_dispersion_name <- function(dispersion) {
  check_column_name(dispersion, "dispersion")
}

check_degrees_freedom_name <- function(dispersion_degrees_freedom) {
  check_column_name(dispersion_degrees_freedom, "dispersion_degrees_freedom")
}

check_shape_prior <- function(shape_prior) {
  if (!is.character(shape_prior) || length(shape_prior) < 1L) {
    stop('`shape_prior` must be "student_t" or "gamma".', call. = FALSE)
  }
  match.arg(shape_prior[[1]], c("student_t", "gamma"))
}

format_prior_number <- function(x) {
  format(x, digits = 8, scientific = FALSE)
}

formula_text <- function(formula) {
  if (inherits(formula, "brmsformula")) {
    formula <- formula$formula
  }
  paste(deparse(formula), collapse = " ")
}

formula_has_offset <- function(formula) {
  grepl("offset\\s*\\(", formula_text(formula))
}

formula_has_lhs <- function(formula) {
  f <- if (inherits(formula, "brmsformula")) formula$formula else formula
  length(f) == 3L
}

add_response <- function(formula, abundance) {
  if (formula_has_lhs(formula)) {
    return(formula)
  }
  if (inherits(formula, "brmsformula")) {
    stop(
      "brmsformula must include a response (e.g. counts ~ ...).",
      call. = FALSE
    )
  }
  stats::update(formula, stats::as.formula(paste(abundance, "~ .")))
}

add_offset_term <- function(formula, offset) {
  if (formula_has_offset(formula) || is.null(offset)) {
    return(formula)
  }
  rhs <- paste("~ . + offset(", offset, ")")
  if (inherits(formula, "brmsformula")) {
    formula$formula <- stats::update(formula$formula, stats::as.formula(rhs))
    return(formula)
  }
  stats::update(formula, stats::as.formula(rhs))
}

add_dispersion_shape <- function(formula, dispersion) {
  if (is.null(dispersion) || has_shape_submodel(formula)) {
    return(formula)
  }
  shape_f <- stats::as.formula(
    sprintf("shape ~ 1 + offset(log(1/%s))", dispersion)
  )
  if (inherits(formula, "brmsformula")) {
    formula$pforms$shape <- shape_f
    return(formula)
  }
  brms::bf(formula, shape_f)
}

prepare_formula <- function(formula,
                            abundance,
                            offset,
                            dispersion = NULL,
                            shape_prior = "student_t") {
  formula <- add_response(formula, abundance)
  formula <- add_offset_term(formula, offset)
  if (identical(check_shape_prior(shape_prior), "gamma")) {
    # The gamma prior is placed on `shape` directly, so no submodel is built.
    return(formula)
  }
  add_dispersion_shape(formula, dispersion)
}

is_zinb_family <- function(family) {
  fam <- family
  if (is.function(fam)) {
    fam <- fam()
  }
  identical(fam$family, "zero_inflated_negbinomial")
}

intercept_location <- function(data, abundance, offset) {
  counts <- data[[abundance]]
  off <- if (is.null(offset)) 0 else data[[offset]]
  if (is.null(off)) {
    stop("Offset column '", offset, "' was not found in `data`.", call. = FALSE)
  }
  location <- mean(log1p(counts / exp(off)))
  if (!is.finite(location)) {
    stop(
      "The intercept prior location is ", location,
      ", computed from '", abundance, "' and '", offset,
      "'. Check those columns for missing or infinite values.",
      call. = FALSE
    )
  }
  location
}

has_shape_submodel <- function(formula) {
  inherits(formula, "brmsformula") && !is.null(formula$pforms$shape)
}

shape_submodel_has_terms <- function(formula) {
  if (!has_shape_submodel(formula)) {
    return(FALSE)
  }
  labels <- attr(stats::terms(formula$pforms$shape), "term.labels")
  length(setdiff(labels, grep("^offset\\(", labels, value = TRUE))) > 0L
}

check_student_df <- function(nu, arg = "shape_prior_df") {
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 2) {
    stop(
      sprintf(
        "`%s` must be a single finite number greater than 2; the Student-t %s",
        arg,
        "standard deviation is undefined for df <= 2."
      ),
      call. = FALSE
    )
  }
  as.numeric(nu)
}

# SD of log(dispersion) implied by d_eff effective degrees of freedom.
# Var(log s^2) = trigamma(d/2) exactly for s^2 ~ sigma^2 chi^2_d / d;
# sqrt(2 / d) is the large-d approximation and is too small at low d.
dispersion_log_sd <- function(d_eff) {
  sqrt(trigamma(d_eff / 2))
}

# Scale of a Student-t whose standard deviation equals `sd`. For df = nu > 2
# the Student-t SD is scale * sqrt(nu / (nu - 2)), so the scale that brms and
# Stan expect as the third argument is sd * sqrt((nu - 2) / nu).
student_t_scale_for_sd <- function(sd, nu) {
  sd * sqrt((nu - 2) / nu)
}

# Prior scale for the shape intercept, derived from the effective degrees of
# freedom that estimate_dispersion() recorded. `default` applies only when no
# dispersion was requested at all, which is the case of a user-supplied shape
# submodel. Once dispersion is in play the degrees of freedom must be there
# too, otherwise the scale would silently stop depending on the data.
shape_intercept_scale <- function(data,
                                  dispersion,
                                  dispersion_degrees_freedom,
                                  nu,
                                  default = 1) {
  if (is.null(dispersion)) {
    return(default)
  }
  require_degrees_freedom_column(data, dispersion, dispersion_degrees_freedom)
  d_eff <- check_degrees_freedom_value(data, dispersion_degrees_freedom)
  scale <- student_t_scale_for_sd(dispersion_log_sd(d_eff), nu)
  if (!is.finite(scale) || scale <= 0) {
    stop(
      "The Student-t scale implied by ", dispersion_degrees_freedom, " = ",
      d_eff, " is ", scale, ", which is not a usable prior scale.",
      call. = FALSE
    )
  }
  scale
}

require_degrees_freedom_column <- function(data,
                                           dispersion,
                                           dispersion_degrees_freedom) {
  if (!is.null(dispersion_degrees_freedom) &&
      dispersion_degrees_freedom %in% names(data)) {
    return(invisible(NULL))
  }
  stop(
    "Dispersion column '", dispersion, "' was supplied without its degrees ",
    "of freedom column '", dispersion_degrees_freedom,
    "'. estimate_dispersion() writes both; pass `dispersion_degrees_freedom` ",
    "if you named it something else.",
    call. = FALSE
  )
}

# edgeR returns prior.df = Inf when robust = TRUE shrinks a gene completely,
# and NA when it could not fit the design at all. Neither yields a prior, so
# report it instead of quietly substituting a default.
check_degrees_freedom_value <- function(data, dispersion_degrees_freedom) {
  d_eff <- as.numeric(data[[dispersion_degrees_freedom]])[[1]]
  if (!is.finite(d_eff) || d_eff <= 0) {
    stop(
      "Column '", dispersion_degrees_freedom, "' is ", d_eff,
      "; effective degrees of freedom must be finite and positive to build ",
      "a shape prior.",
      call. = FALSE
    )
  }
  d_eff
}

# Conjugate gamma prior on the brms `shape` parameter, i.e. on 1/phi.
#
# edgeR and limma model the gene-wise dispersion as a scaled inverse
# chi-square with d_eff degrees of freedom centred on phi_g. Inverting a
# scaled inverse chi-square gives a gamma, so the precision 1/phi -- exactly
# what brms calls `shape` -- has prior Gamma(d_eff/2, rate = d_eff phi_g / 2).
# Its mean is 1/phi_g, edgeR's point estimate, and because
# Var(log X) = trigamma(shape) for a gamma, its log-scale spread is
# trigamma(d_eff/2): identical to the Student-t route in
# shape_intercept_scale(), which is unsurprising since chi^2_d is itself
# Gamma(d/2, scale = 2). The two forms differ only in tail weight, the
# Student-t being the more robust to a badly shrunk edgeR estimate.
#
# Falls back to brms' own vague gamma(0.01, 0.01) only when there is no
# dispersion to build a prior from; values that are present but unusable are
# an error.
shape_gamma_parameters <- function(data,
                                   dispersion,
                                   dispersion_degrees_freedom,
                                   default = list(shape = 0.01, rate = 0.01)) {
  if (is.null(dispersion)) {
    return(default)
  }
  if (!dispersion %in% names(data)) {
    stop(
      "Dispersion column '", dispersion, "' was not found in `data`.",
      call. = FALSE
    )
  }
  require_degrees_freedom_column(data, dispersion, dispersion_degrees_freedom)
  d_eff <- check_degrees_freedom_value(data, dispersion_degrees_freedom)
  data <- check_dispersion_values(data, dispersion)
  phi <- as.numeric(data[[dispersion]])[[1]]
  list(shape = d_eff / 2, rate = d_eff * phi / 2)
}

shape_student_t_prior <- function(data,
                                  formula,
                                  dispersion,
                                  dispersion_degrees_freedom,
                                  shape_prior_df,
                                  shape_prior = "student_t") {
  nu <- check_student_df(shape_prior_df)
  if (identical(shape_prior, "gamma")) {
    stop(
      'shape_prior = "gamma" puts a prior on a scalar `shape`, but `formula` ',
      "already has a shape submodel, which has no such parameter. Drop the ",
      'submodel or use shape_prior = "student_t".',
      call. = FALSE
    )
  }
  scale <- shape_intercept_scale(data, dispersion, dispersion_degrees_freedom, nu)
  p <- brms::prior_string(
    sprintf(
      "student_t(%s, 0, %s)",
      format_prior_number(nu),
      format_prior_number(scale)
    ),
    class = "Intercept",
    dpar = "shape"
  )
  if (shape_submodel_has_terms(formula)) {
    p <- c(p, brms::prior(student_t(3, 0, 2), class = b, dpar = shape))
  }
  p
}

zinb_location_priors <- function(data, abundance, offset) {
  i <- intercept_location(data, abundance, offset)
  c(
    brms::prior_string(
      sprintf("student_t(3, %s, 1.5)", format_prior_number(i)),
      class = "Intercept"
    ),
    brms::prior(student_t(3, 0, 5), class = b)
  )
}

shape_gamma_prior <- function(data, dispersion, dispersion_degrees_freedom) {
  pars <- shape_gamma_parameters(data, dispersion, dispersion_degrees_freedom)
  brms::prior_string(
    sprintf(
      "gamma(%s, %s)",
      format_prior_number(pars$shape),
      format_prior_number(pars$rate)
    ),
    class = "shape"
  )
}

detect_cores <- function() {
  if (requireNamespace("parallelly", quietly = TRUE)) {
    as.numeric(parallelly::availableCores())
  } else {
    as.numeric(parallel::detectCores(logical = TRUE))
  }
}

make_gene_inits <- function(formula, data, family, prior, chains, abundance, offset) {
  sdata <- brms::make_standata(
    formula = formula,
    data = data,
    family = family,
    prior = prior
  )
  mu <- intercept_location(data, abundance, offset)
  lapply(seq_len(chains), function(i) {
    inits <- list(
      Intercept = stats::rnorm(1, mu, 1.5)
    )
    if (!is.null(sdata$Kc) && isTRUE(sdata$Kc > 0)) {
      inits$b <- stats::rnorm(sdata$Kc, 0, 5)
    }
    if (!is.null(sdata$Kc_shape)) {
      inits$b_shape <- stats::rnorm(sdata$Kc_shape, 0, 2)
      inits$Intercept_shape <- stats::rnorm(1, 0, 1)
    } else if (grepl("negbinomial", family$family, fixed = TRUE)) {
      inits$shape <- stats::rexp(1, rate = 1)
    }
    if (grepl("zero_inflated", family$family, fixed = TRUE)) {
      inits$zi <- stats::runif(1, 0.01, 0.2)
    }
    inits
  })
}

copy_rowdata_columns <- function(out, se, columns, optional = NULL) {
  if (length(columns) == 0L && length(optional) == 0L) {
    return(out)
  }
  rd <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)
  missing <- setdiff(columns, names(rd))
  if (length(missing) > 0L) {
    stop(
      "rowData column '", paste(missing, collapse = "', '"), "' was not found.",
      call. = FALSE
    )
  }
  for (col in c(columns, intersect(optional, names(rd)))) {
    out[[col]] <- rd[[col]][[1]]
  }
  out
}

as_gene_tibble <- function(data,
                           abundance = "counts",
                           rowdata_cols = NULL,
                           rowdata_cols_optional = NULL) {
  if (inherits(data, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Install SummarizedExperiment to pass a SummarizedExperiment.", call. = FALSE)
    }
    n_feature <- nrow(data)
    if (n_feature != 1L) {
      stop(
        "Gene-wise functions expect one gene at a time (nrow(data) == 1). ",
        "Got ", n_feature, " feature(s).",
        call. = FALSE
      )
    }
    assay_names <- SummarizedExperiment::assayNames(data)
    if (!abundance %in% assay_names) {
      stop(
        "Assay '", abundance, "' not found. Available: ",
        paste(assay_names, collapse = ", "),
        call. = FALSE
      )
    }
    counts <- as.vector(SummarizedExperiment::assay(data, abundance)[1, , drop = TRUE])
    out <- as.data.frame(SummarizedExperiment::colData(data), stringsAsFactors = FALSE)
    out[[abundance]] <- counts
    if (is.null(out[[".feature"]])) {
      out[[".feature"]] <- rownames(data)[[1]]
    }
    out <- copy_rowdata_columns(
      out,
      data,
      rowdata_cols,
      optional = rowdata_cols_optional
    )
    return(tibble::as_tibble(out))
  }

  if (is.matrix(data) || is.array(data)) {
    stop("`data` must be a data frame or a one-gene SummarizedExperiment.", call. = FALSE)
  }

  tibble::as_tibble(data)
}

check_dispersion_values <- function(data, dispersion) {
  if (is.null(dispersion)) {
    return(data)
  }
  if (!dispersion %in% names(data)) {
    stop("Dispersion column '", dispersion, "' was not found in `data`.", call. = FALSE)
  }
  val <- as.numeric(data[[dispersion]])
  bad <- !is.finite(val) | val <= 0
  if (any(bad)) {
    stop(
      "Dispersion column '", dispersion, "' has ", sum(bad),
      " non-finite or non-positive value(s). estimate_dispersion() returns ",
      "NA when edgeR cannot fit the design; fix that rather than passing ",
      "the result on.",
      call. = FALSE
    )
  }
  data[[dispersion]] <- val
  data
}

prepare_gene_data <- function(data,
                              abundance = "counts",
                              offset,
                              dispersion = NULL,
                              dispersion_degrees_freedom = NULL,
                              sanitize_names = FALSE) {
  offset <- check_offset_name(offset)
  if (!is.null(dispersion)) {
    dispersion <- check_dispersion_name(dispersion)
  }
  if (!is.null(dispersion_degrees_freedom)) {
    dispersion_degrees_freedom <-
      check_degrees_freedom_name(dispersion_degrees_freedom)
  }
  data <- as_gene_tibble(
    data,
    abundance = abundance,
    rowdata_cols = dispersion,
    rowdata_cols_optional = dispersion_degrees_freedom
  )

  if (!abundance %in% names(data)) {
    stop("Abundance column '", abundance, "' was not found in `data`.", call. = FALSE)
  }

  data[[abundance]] <- as.integer(data[[abundance]])

  n_na <- sum(is.na(data[[abundance]]))
  if (n_na > 0) {
    stop(
      glue::glue(
        "Column '{abundance}' has {n_na} NA(s). Dropping those samples would ",
        "change the model silently, so decide what to do with them before ",
        "calling estimate_gene()."
      ),
      call. = FALSE
    )
  }

  if (!is.null(offset) && !offset %in% names(data)) {
    stop("Offset column '", offset, "' was not found in `data`.", call. = FALSE)
  }
  data <- check_dispersion_values(data, dispersion)

  if (isTRUE(sanitize_names)) {
    names(data) <- collapse_repeated_underscores(names(data))
    abundance <- collapse_repeated_underscores(abundance)
    if (!is.null(offset)) {
      offset <- collapse_repeated_underscores(offset)
    }
    if (!is.null(dispersion)) {
      dispersion <- collapse_repeated_underscores(dispersion)
    }
    if (!is.null(dispersion_degrees_freedom)) {
      dispersion_degrees_freedom <-
        collapse_repeated_underscores(dispersion_degrees_freedom)
    }
  }

  list(
    data = droplevels(data),
    abundance = abundance,
    offset = offset,
    dispersion = dispersion,
    dispersion_degrees_freedom = dispersion_degrees_freedom
  )
}

random_intercept_parameters <- function(fit, grouping, par = "Intercept") {
  pattern <- paste0(
    "^r_", grouping, "\\[.*,", par, "\\]$"
  )
  params <- rownames(summary(fit$fit)[[1]])
  params[grepl(pattern, params)]
}

random_intercept_vs_rest_from_names <- function(params,
                                                grouping,
                                                par = "Intercept") {
  if (length(params) < 2L) {
    stop(
      "Need at least two '", grouping, "' levels to contrast each level ",
      "against the rest.",
      call. = FALSE
    )
  }

  quoted <- paste0("`", sub("^r_", "", params), "`")
  level_names <- sub(
    paste0("^`", grouping, "\\[(.*),", par, "\\]`$"),
    "\\1",
    quoted
  )

  equations <- vapply(seq_along(quoted), function(i) {
    this_param <- quoted[[i]]
    other_params <- quoted[-i]
    avg_expr <- paste0(
      "(", paste(other_params, collapse = " + "), ")/",
      length(other_params)
    )
    paste0(this_param, " - ", avg_expr, " = 0")
  }, character(1))
  stats::setNames(equations, level_names)
}

random_intercept_vs_rest_equations <- function(fit,
                                               grouping,
                                               par = "Intercept") {
  random_intercept_vs_rest_from_names(
    random_intercept_parameters(fit, grouping, par = par),
    grouping = grouping,
    par = par
  )
}

nullify_newdata <- function(data, nullify = NULL, offset = "offset", offset_value = 0) {
  newdata <- data
  if (!is.null(offset) && offset %in% names(newdata)) {
    newdata[[offset]] <- offset_value
  }
  missing <- setdiff(nullify, names(newdata))
  if (length(missing) > 0L) {
    stop(
      "Cannot nullify column(s) '", paste(missing, collapse = "', '"),
      "': not in the model data. Available: ",
      paste(names(newdata), collapse = ", "),
      call. = FALSE
    )
  }
  for (nm in nullify) {
    newdata[[nm]] <- NA
  }
  newdata
}

tidy_hypothesis <- function(hyp, component, grouping_label) {
  tbl <- tibble::as_tibble(hyp$hypothesis)
  if (!nrow(tbl)) {
    return(tbl)
  }
  dplyr::transmute(
    tbl,
    component = component,
    group = if ("Group" %in% names(tbl)) .data$Group else grouping_label,
    hypothesis = .data$Hypothesis,
    estimate = .data$Estimate,
    ci_lower = .data$CI.Lower,
    ci_upper = .data$CI.Upper,
    post_prob = if ("Post.Prob" %in% names(tbl)) .data$Post.Prob else NA_real_,
    evid_ratio = if ("Evid.Ratio" %in% names(tbl)) .data$Evid.Ratio else NA_real_
  )
}
