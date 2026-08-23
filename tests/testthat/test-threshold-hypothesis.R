# The directional test is the default route through hypothesis_gene(). These
# tests pin down the arithmetic against draws counted by hand, the properties
# that make it prefer a one-sided claim over a two-sided one, and the false
# discovery rate built on top of it.

test_that("a bare contrast takes the directional route and '= 0' the Bayes factor", {
  expect_equal(
    brmDE:::hypothesis_is_equal_zero(c(
      "dextrt",
      "`cell[N1,Intercept]` - `cell[N2,Intercept]`",
      "dextrt = 0",
      "dextrt =0.0",
      "dextrt - dexb = 0"
    )),
    c(FALSE, FALSE, TRUE, TRUE, TRUE)
  )

  # Names are the caller's labels, not part of the equation.
  expect_equal(brmDE:::hypothesis_is_equal_zero(c(treated = "dextrt")), FALSE)
})

test_that("a right-hand side other than zero is refused", {
  # brms returns draws of `left - right`, so a non-zero right-hand side would
  # shift the fold change rather than change the question.
  for (bad in c("dextrt = 0.7", "dextrt > 0.7", "dextrt < -0.7", "dextrt > 0")) {
    expect_error(brmDE:::hypothesis_is_equal_zero(bad), "Unsupported hypothesis")
  }

  # Reported together, and the supported entries alongside them are not.
  expect_error(
    brmDE:::hypothesis_is_equal_zero(c("dextrt", "dexb > 0.7", "dexc = 1")),
    '"dexb > 0.7", "dexc = 1"',
    fixed = TRUE
  )
})

test_that("every route reports the same effect size and the same vocabulary", {
  skip_cmdstan()

  fit <- estimate_gene(
    airway_one_gene("ENSG00000120129"),
    formula_abundance = ~ dex,
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    sample_prior = "yes",
    chains = 2,
    draws_warmup = 500,
    draws_sampling = 500,
    cores = 2,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  routes <- list(
    posterior_draws = hypothesis_gene(fit, "dextrt"),
    bayes_factor = hypothesis_gene(fit, "dextrt = 0")
  )

  # One schema, so tables from different tests stack and rank together.
  schema <- lapply(routes, names)
  expect_true(all(vapply(schema, identical, logical(1), schema[[1]])))
  expect_equal(
    vapply(routes, function(x) x$pH0_from, character(1)),
    c(posterior_draws = "posterior_draws", bayes_factor = "bayes_factor")
  )

  # The effect size describes the contrast, so asking a different question
  # about it must not move it.
  draws <- as.data.frame(fit)$b_dextrt
  for (route in routes) {
    expect_equal(route$estimate, stats::median(draws))
    expect_equal(route$log2_fold_change, stats::median(draws) / log(2))
  }

  # pH0 is the null on both routes, so it is low for both here, and
  # evid_ratio is always the odds against that null.
  for (route in routes) {
    expect_lt(route$pH0, 0.1)
    expect_equal(route$evid_ratio, (1 - route$pH0) / route$pH0)
  }

  # A single call may mix routes, one row each, from one brms call.
  mixed <- hypothesis_gene(fit, c("dextrt", "dextrt = 0"))
  expect_equal(mixed$pH0_from, c("posterior_draws", "bayes_factor"))
  expect_equal(mixed$hypothesis, c("dextrt", "dextrt = 0"))
  expect_equal(mixed$log2_fold_change, rep(stats::median(draws) / log(2), 2))
  expect_equal(
    mixed$pH0,
    c(routes$posterior_draws$pH0, routes$bayes_factor$pH0)
  )
})

test_that("pH0 is counted from the draws", {
  skip_cmdstan()

  fit <- estimate_gene(
    airway_one_gene("ENSG00000120129"),
    formula_abundance = ~ dex,
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    chains = 2,
    draws_warmup = 500,
    draws_sampling = 500,
    cores = 2,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  draws <- as.data.frame(fit)$b_dextrt
  out <- hypothesis_gene(fit, "dextrt", test_above_log2FC = 1)

  # The argument is in log2 units, so the draws are compared against log(2).
  eps <- log(2)
  expect_equal(nrow(out), 1L)
  expect_equal(out$hypothesis, "dextrt")
  expect_equal(out$pH0, 1 - max(mean(draws > eps), mean(draws < -eps)))

  # Two log2 fold changes is 2 * log(2), not 2.
  two <- hypothesis_gene(fit, "dextrt", test_above_log2FC = 2)
  expect_equal(
    two$pH0,
    1 - max(mean(draws > 2 * log(2)), mean(draws < -2 * log(2)))
  )
  expect_equal(out$estimate, stats::median(draws))
  expect_equal(out$log2_fold_change, stats::median(draws) / log(2))

  # Raising the threshold past the whole posterior has to send pH0 to 1, and
  # dropping it to 0 leaves the local false sign rate.
  far <- hypothesis_gene(fit, "dextrt", test_above_log2FC = 50)
  expect_equal(far$pH0, 1)

  # A pH0 of 1 says no direction is supported, but the effect size is still
  # reported, so the sign stays readable off log2_fold_change alone.
  expect_false("direction" %in% names(far))
  expect_gt(far$log2_fold_change, 0)
  expect_lt(
    hypothesis_gene(fit, "-dextrt", test_above_log2FC = 50)$log2_fold_change,
    0
  )
  lfsr <- hypothesis_gene(fit, "dextrt", test_above_log2FC = 0)
  expect_equal(lfsr$pH0, min(mean(draws > 0), mean(draws < 0)))

  # The interval is on the coefficient scale and does not move with the
  # threshold, which only decides the probabilities.
  expect_equal(out$ci_lower, lfsr$ci_lower)
  expect_equal(out$ci_upper, lfsr$ci_upper)
  expect_lt(out$ci_lower, out$estimate)
  expect_gt(out$ci_upper, out$estimate)
})

test_that("a posterior straddling zero cannot be a discovery", {
  # This is the property that testing abs(effect) > threshold would lose: a
  # heavy-tailed posterior centred on zero has mass beyond the threshold on
  # both sides, but its sign is undetermined, so pH0 must stay near 1/2.
  set.seed(1)
  draws <- data.frame(effect = stats::rt(20000, df = 2) * 3)
  hyp <- brms::hypothesis(draws, "effect > 0.7")
  two_sided <- mean(abs(draws$effect) > 0.7)

  p_positive <- mean(draws$effect > 0.7)
  p_negative <- mean(draws$effect < -0.7)
  pH0 <- 1 - max(p_positive, p_negative)

  # The two-sided reading calls this a discovery; the directional one refuses.
  expect_gt(two_sided, 0.6)
  expect_gt(pH0, 0.5)

  # For a posterior symmetric about zero, pH0 = (1 + P(ROPE))/2, so it can
  # never fall below 1/2 however heavy the tails.
  rope <- mean(abs(draws$effect) <= 0.7)
  expect_equal(pH0, (1 + rope) / 2, tolerance = 0.01)
  expect_equal(p_positive, hyp$hypothesis$Post.Prob, tolerance = 1e-8)
})

test_that("false_discovery_rate averages null probabilities in place", {
  pH0 <- c(0.4, 0.001, 0.02, 0.9, 0.005)
  fdr <- false_discovery_rate(pH0)

  # Returned in the order given, not sorted.
  expect_equal(fdr[[2]], 0.001)
  expect_equal(fdr[[5]], mean(c(0.001, 0.005)))
  expect_equal(fdr[[3]], mean(c(0.001, 0.005, 0.02)))
  expect_equal(fdr[[1]], mean(c(0.001, 0.005, 0.02, 0.4)))
  expect_equal(fdr[[4]], mean(pH0))

  # A cumulative mean over an ascending sequence is already monotone, so no
  # Benjamini-Hochberg style step-up correction is needed.
  expect_true(!is.unsorted(fdr[order(pH0)]))
  expect_true(all(fdr <= max(pH0)))
  expect_equal(false_discovery_rate(rep(0, 5)), rep(0, 5))
})

test_that("the pipeline adds an fdr column across genes", {
  tables <- list(
    tibble::tibble(component = "fixed", hypothesis = "dextrt", pH0 = 0.4),
    tibble::tibble(component = "fixed", hypothesis = "dextrt", pH0 = 0.01),
    tibble::tibble(component = "fixed", hypothesis = "dextrt", pH0 = 0.2)
  )
  out <- brmDE:::add_hypothesis_fdr(tables)

  expect_length(out, 3L)
  expect_true(all(vapply(out, nrow, integer(1)) == 1L))
  expect_equal(vapply(out, function(x) x$pH0, numeric(1)), c(0.4, 0.01, 0.2))
  expect_equal(
    vapply(out, function(x) x$fdr, numeric(1)),
    c(mean(c(0.01, 0.2, 0.4)), 0.01, mean(c(0.01, 0.2)))
  )

  # Each hypothesis is its own family of tests.
  two <- lapply(c(0.4, 0.01), function(p) {
    tibble::tibble(
      component = "fixed",
      hypothesis = c("a", "b"),
      pH0 = c(p, p)
    )
  })
  out_two <- brmDE:::add_hypothesis_fdr(two)
  expect_equal(out_two[[1]]$fdr, rep(mean(c(0.01, 0.4)), 2))
  expect_equal(out_two[[2]]$fdr, c(0.01, 0.01))

  # Every route now reports pH0, so the FDR reaches all of them, including one
  # whose pH0 came from a Bayes factor rather than from the draws.
  from_bf <- lapply(c(0.3, 0.004), function(p) {
    tibble::tibble(component = "fixed", hypothesis = "dextrt = 0", pH0 = p)
  })
  expect_equal(
    vapply(brmDE:::add_hypothesis_fdr(from_bf), function(x) x$fdr, numeric(1)),
    c(mean(c(0.004, 0.3)), 0.004)
  )

  # A table with no pH0 at all is left alone rather than gaining an empty column.
  without <- list(
    tibble::tibble(hypothesis = "dextrt", estimate = 0.9),
    tibble::tibble(hypothesis = "dextrt", estimate = 0.5)
  )
  expect_identical(brmDE:::add_hypothesis_fdr(without), without)
})

test_that("an '= 0' test without prior draws leaves fdr NA without spoiling others", {
  # pH0 is NA when brms could not compute a Bayes factor. Those genes cannot be
  # ranked, but they must not corrupt the rate for the genes that can be.
  expect_equal(false_discovery_rate(c(0.1, NA, 0.3)), c(0.1, NA, 0.2))
  expect_true(all(is.na(false_discovery_rate(c(NA_real_, NA_real_)))))
})
