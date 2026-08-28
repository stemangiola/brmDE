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

check_dispersion_name <- function(dispersion_prior_log_mean) {
  check_column_name(dispersion_prior_log_mean, "dispersion_prior_log_mean")
}

check_log_sd_name <- function(dispersion_prior_log_sd) {
  check_column_name(dispersion_prior_log_sd, "dispersion_prior_log_sd")
}

# Suggested Bioconductor packages (edgeR, limma, ...). Same helper as tidybulk:
# offer to install via BiocManager rather than failing on requireNamespace().
check_and_install_packages <- function(packages) {
  rlang::check_installed(
    pkg = packages,
    action = function(...) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(..., ask = FALSE, update = FALSE)
    }
  )
}

check_bundle <- function(bundle) {
  if (!is.numeric(bundle) || length(bundle) != 1L || !is.finite(bundle) ||
    bundle < 1 || bundle != trunc(bundle)) {
    stop(
      "`bundle` must be a single positive whole number of genes per target; ",
      "1 fits one gene per target.",
      call. = FALSE
    )
  }
  as.integer(bundle)
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

formula_rhs_text <- function(formula, arg) {
  if (inherits(formula, "brmsformula")) {
    formula <- formula$formula
  }
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }
  if (!inherits(formula, "formula")) {
    stop(sprintf("`%s` must be a formula.", arg), call. = FALSE)
  }
  if (length(formula) == 3L) {
    stop(
      sprintf(
        "`%s` must be one-sided (e.g. ~1 or ~ dex); the response is set for you.",
        arg
      ),
      call. = FALSE
    )
  }
  deparse1(formula[[2]])
}

# The shape submodel is assembled here rather than being handed in whole, so
# that the edgeR dispersion always enters as an offset on brms' log link.
# No dispersion means that offset is 0: the intercept is then log(shape)
# itself, with prior median shape = 1.
add_dispersion_shape <- function(formula, formula_dispersion, dispersion) {
  rhs <- formula_rhs_text(formula_dispersion, "formula_dispersion")
  added_offset <- !formula_has_offset(formula_dispersion)
  if (added_offset) {
    offset_term <- if (is.null(dispersion)) {
      "offset(0)"
    } else {
      sprintf("offset(log(1/%s))", dispersion)
    }
    rhs <- sprintf("%s + %s", rhs, offset_term)
  }
  shape_f <- stats::as.formula(paste("shape ~", rhs))
  announce_formula("Dispersion model", shape_f, added_offset)

  if (inherits(formula, "brmsformula")) {
    formula$pforms$shape <- shape_f
    return(formula)
  }
  brms::bf(formula, shape_f)
}

announce_formula <- function(label, formula, added_offset) {
  message(
    label,
    if (added_offset) " (offset added by brmDE)" else "",
    ": ",
    formula_text(formula)
  )
}

prepare_formula <- function(formula_abundance,
                            formula_dispersion = ~1,
                            abundance,
                            offset,
                            dispersion = NULL) {
  if (has_shape_submodel(formula_abundance)) {
    stop(
      "`formula_abundance` carries a shape submodel. Model the dispersion ",
      "through `formula_dispersion` instead, so that the edgeR dispersion ",
      "offset is added to it.",
      call. = FALSE
    )
  }
  formula_abundance <- add_response(formula_abundance, abundance)
  added_offset <- !is.null(offset) && !formula_has_offset(formula_abundance)
  formula_abundance <- add_offset_term(formula_abundance, offset)
  announce_formula("Abundance model", formula_abundance, added_offset)

  # The edgeR dispersion is an offset on brms' log link, or 0 when it is
  # omitted. The Student-t prior on that intercept is set later; its width
  # comes from `dispersion_prior_log_sd`.
  out <- add_dispersion_shape(formula_abundance, formula_dispersion, dispersion)
  strip_formula_env(out)
}

# A formula holds a reference to the frame it was written in, so a model fitted
# inside a function serialises that entire frame - the SummarizedExperiment
# included - into every stored fit. Terms are resolved against the data first
# and the global environment after, which is where user-defined helpers live.
# This also makes a direct estimate_gene() call match the pipeline, which has
# always rebuilt formulas against globalenv() after a round trip through disk.
strip_formula_env <- function(x) {
  if (inherits(x, "brmsformula")) {
    x$formula <- strip_formula_env(x$formula)
    if (length(x$pforms)) {
      x$pforms <- lapply(x$pforms, strip_formula_env)
    }
    return(x)
  }
  if (inherits(x, "formula")) {
    environment(x) <- globalenv()
  }
  x
}

is_negbinomial_family <- function(family) {
  fam <- family
  if (is.function(fam)) {
    fam <- fam()
  }
  grepl("negbinomial", fam$family, fixed = TRUE)
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

# Scale of a Student-t whose standard deviation equals `sd`. For df = nu > 2
# the Student-t SD is scale * sqrt(nu / (nu - 2)), so the scale that brms and
# Stan expect as the third argument is sd * sqrt((nu - 2) / nu).
student_t_scale_for_sd <- function(sd, nu) {
  sd * sqrt((nu - 2) / nu)
}

# Log-scale SD of the Student-t shape intercept when no log-SD column is
# supplied. Independent of `dispersion`: that argument only shifts the offset.
# Converted to a Student-t scale through student_t_scale_for_sd(), so nu still
# controls both the df and the scale.
shape_prior_sd_default <- 1

check_log_sd_value <- function(data, dispersion_prior_log_sd) {
  if (!dispersion_prior_log_sd %in% names(data)) {
    stop(
      "Dispersion log-SD column '", dispersion_prior_log_sd,
      "' was not found in `data`. Run estimate_dispersion_prior() to write it, ",
      "or omit `dispersion_prior_log_sd` to use a log-scale SD of 1.",
      call. = FALSE
    )
  }
  sd <- as.numeric(data[[dispersion_prior_log_sd]])[[1]]
  if (!is.finite(sd) || sd <= 0) {
    stop(
      "Column '", dispersion_prior_log_sd, "' is ", sd,
      "; the log-dispersion prior SD must be finite and positive.",
      call. = FALSE
    )
  }
  sd
}

shape_student_t_prior <- function(data,
                                  formula,
                                  sd,
                                  shape_prior_df) {
  nu <- check_student_df(shape_prior_df)
  scale <- student_t_scale_for_sd(sd, nu)
  if (!is.finite(scale) || scale <= 0) {
    stop(
      "The Student-t scale implied by log-scale SD ", sd, " is ", scale,
      ", which is not a usable prior scale.",
      call. = FALSE
    )
  }
  p <- brms::prior_string(
    sprintf("student_t(%s, 0, brmde_shape_scale)", format_prior_number(nu)),
    class = "Intercept",
    dpar = "shape"
  )
  if (shape_submodel_has_terms(formula)) {
    p <- c(
      p,
      brms::prior_string("student_t(3, 0, 2)", class = "b", dpar = "shape")
    )
  }
  gene_prior(p, gene_prior_stanvar("brmde_shape_scale", scale))
}

# Default scale of the Student-t prior on the abundance coefficients. They sit
# behind a log link, so each is a natural-log fold change, and 0.7 is one log2
# fold change: a doubling of expression. Multiples follow in log2 units, 1.4
# being two log2 fold changes.
#
# This scale is not only a shrinkage choice. An "= 0" hypothesis ("dextrt = 0")
# is evaluated by the Savage-Dickey density ratio, the posterior density at 0
# divided by the prior density at 0, so widening this prior thins prior mass at
# 0 and shifts the evidence towards the null however clear the data are
# (Lindley's paradox). A scale on the order of the effects being looked for is
# what makes that ratio mean anything.
#
# estimate_gene() spells the same default out in its own signature, where a
# user reading ?estimate_gene can see it.
coefficient_prior_scale_default <- 0.7

location_priors <- function(data,
                            abundance,
                            offset,
                            coefficient_prior_scale =
                              coefficient_prior_scale_default,
                            coefficient_prior_df = 3) {
  i <- intercept_location(data, abundance, offset)
  gene_prior(
    prior = c(
      brms::prior_string(
        "student_t(3, brmde_intercept_location, 1.5)",
        class = "Intercept"
      ),
      brms::prior_string(
        sprintf(
          "student_t(%s, 0, %s)",
          format_prior_number(coefficient_prior_df),
          format_prior_number(coefficient_prior_scale)
        ),
        class = "b"
      )
    ),
    stanvars = gene_prior_stanvar("brmde_intercept_location", i)
  )
}

# Prior constants derived from a gene's own data are passed to Stan as data
# rather than pasted into the model code as literals. The generated code is
# then byte-identical for every gene, so cmdstanr compiles the model once per
# process and reuses it for every other gene. Constants that come from user
# arguments, such as the Student-t degrees of freedom, stay literal because
# they cannot vary from gene to gene.
gene_prior_stanvar <- function(name, value) {
  brms::stanvar(x = as.numeric(value), name = name, block = "data")
}

# A brmsprior together with the stanvars that define the symbols it refers to.
# The two always travel together: a prior naming `brmde_intercept_location`
# will not compile unless the matching stanvar reaches brms as well.
gene_prior <- function(prior, stanvars = NULL) {
  list(prior = prior, stanvars = stanvars)
}

combine_stanvars <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (length(parts) == 0L) {
    return(NULL)
  }
  Reduce(`+`, parts)
}

combine_gene_priors <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  gene_prior(
    prior = do.call(c, lapply(parts, `[[`, "prior")),
    stanvars = do.call(combine_stanvars, lapply(parts, `[[`, "stanvars"))
  )
}

# The default prior set for a negative binomial gene: one term for the shape,
# one for the location parameters. Every parameter it touches gets a proper
# prior, which is what lets "= 0" hypotheses be tested at all: brms can only
# draw from a prior that integrates to one.
default_gene_priors <- function(data,
                                formula,
                                abundance,
                                offset,
                                dispersion_prior_log_sd,
                                shape_prior_df,
                                coefficient_prior_scale =
                                  coefficient_prior_scale_default,
                                coefficient_prior_df = 3) {
  sd <- if (is.null(dispersion_prior_log_sd)) {
    shape_prior_sd_default
  } else {
    check_log_sd_value(data, dispersion_prior_log_sd)
  }
  combine_gene_priors(
    shape_student_t_prior(data, formula, sd, shape_prior_df),
    location_priors(
      data,
      abundance,
      offset,
      coefficient_prior_scale = coefficient_prior_scale,
      coefficient_prior_df = coefficient_prior_df
    )
  )
}

detect_cores <- function() {
  if (requireNamespace("parallelly", quietly = TRUE)) {
    as.numeric(parallelly::availableCores())
  } else {
    as.numeric(parallel::detectCores(logical = TRUE))
  }
}

make_gene_inits <- function(formula, data, family, prior, chains, abundance,
                            offset, stanvars = NULL) {
  sdata <- brms::make_standata(
    formula = formula,
    data = data,
    family = family,
    prior = prior,
    stanvars = stanvars
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

check_dispersion_values <- function(data, dispersion_prior_log_mean) {
  if (is.null(dispersion_prior_log_mean)) {
    return(data)
  }
  if (!dispersion_prior_log_mean %in% names(data)) {
    stop("Dispersion column '", dispersion_prior_log_mean, "' was not found in `data`.", call. = FALSE)
  }
  val <- as.numeric(data[[dispersion_prior_log_mean]])
  bad <- !is.finite(val) | val <= 0
  if (any(bad)) {
    stop(
      "Dispersion column '", dispersion_prior_log_mean, "' has ", sum(bad),
      " non-finite or non-positive value(s). estimate_dispersion_prior() ",
      "and tidybulk::estimate_dispersion() return NA when the design cannot ",
      "be fit; fix that rather than passing the result on.",
      call. = FALSE
    )
  }
  data[[dispersion_prior_log_mean]] <- val
  data
}

prepare_gene_data <- function(data,
                              abundance = "counts",
                              offset,
                              dispersion_prior_log_mean = NULL,
                              dispersion_prior_log_sd = NULL,
                              sanitize_names = FALSE) {
  offset <- check_offset_name(offset)
  if (!is.null(dispersion_prior_log_mean)) {
    dispersion_prior_log_mean <- check_dispersion_name(dispersion_prior_log_mean)
  }
  if (!is.null(dispersion_prior_log_sd)) {
    dispersion_prior_log_sd <- check_log_sd_name(dispersion_prior_log_sd)
  }
  data <- as_gene_tibble(
    data,
    abundance = abundance,
    rowdata_cols = c(dispersion_prior_log_mean, dispersion_prior_log_sd)
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
  data <- check_dispersion_values(data, dispersion_prior_log_mean)

  if (isTRUE(sanitize_names)) {
    names(data) <- collapse_repeated_underscores(names(data))
    abundance <- collapse_repeated_underscores(abundance)
    if (!is.null(offset)) {
      offset <- collapse_repeated_underscores(offset)
    }
    if (!is.null(dispersion_prior_log_mean)) {
      dispersion_prior_log_mean <-
        collapse_repeated_underscores(dispersion_prior_log_mean)
    }
    if (!is.null(dispersion_prior_log_sd)) {
      dispersion_prior_log_sd <-
        collapse_repeated_underscores(dispersion_prior_log_sd)
    }
  }

  list(
    data = droplevels(data),
    abundance = abundance,
    offset = offset,
    dispersion_prior_log_mean = dispersion_prior_log_mean,
    dispersion_prior_log_sd = dispersion_prior_log_sd
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


split_rows <- function(x, sizes) {
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_along(sizes), function(i) {
    if (sizes[[i]] == 0L) x[0, ] else x[starts[[i]]:ends[[i]], ]
  })
}
