hypotheses <- data.frame(
  .feature = paste0("gene", 1:6),
  log2_fold_change = c(-3.1, 2.4, 0.2, -0.4, 4.5, 1.1),
  pH0 = c(0, 0.001, 0.8, 0.4, 0, 0.02)
)

volcano <- function(...) {
  ggplot2::ggplot_build(brmDE:::volcano_ggplot(hypotheses, ...))
}

test_that("plot_volcano wants the pipeline, which knows the draws", {
  expect_error(plot_volcano(hypotheses), "pipeline from brmDE")
  pipeline <- structure(list(initialisation = list()), class = "brmDE_hpc")
  expect_error(plot_volcano(pipeline), "hypothesis")
})

test_that("resolution is 1 / ndraws, or the smallest probability counted", {
  expect_equal(brmDE:::volcano_resolution(hypotheses$pH0, 1000), 0.001)
  expect_equal(brmDE:::volcano_resolution(hypotheses$pH0, NULL), 0.001)
  expect_equal(brmDE:::volcano_resolution(c(0.02, NA, 0.5), NULL), 0.02)
  expect_error(brmDE:::volcano_resolution(c(0, 0), NULL), "number_of_draws")
})

test_that("unresolved genes are jittered inside the band, the rest are not", {
  built <- volcano(number_of_draws = 1000, seed = 1)
  # The y scale is log10 then reversed, so the panel holds -log10(probability).
  probability_of <- function(layer) 10^(-built$data[[layer]]$y)

  expect_equal(sort(probability_of(3)), c(0.001, 0.02, 0.4, 0.8))
  expect_length(probability_of(4), 2L)
  expect_true(all(probability_of(4) > 1e-4 & probability_of(4) < 1e-3))
})

test_that("the band spans the decade below the resolution", {
  band <- volcano(number_of_draws = 1000)$data[[1]]
  expect_equal(sort(10^(-c(band$ymin, band$ymax))), c(1e-4, 1e-3))
  expect_equal(band$fill, "lightyellow")

  # Shorter fits resolve less, so the band and its line move up with them.
  band <- volcano(number_of_draws = 500)$data[[1]]
  expect_equal(sort(10^(-c(band$ymin, band$ymax))), c(2e-4, 2e-3))
})

test_that("significance colours and sizes follow the threshold", {
  built <- volcano(number_of_draws = 1000, significance_threshold = 0.05)
  points <- rbind(
    built$data[[3]][, c("colour", "size")],
    built$data[[4]][, c("colour", "size")]
  )
  expect_equal(sum(points$colour == "red"), 4L)
  expect_true(all(points$size[points$colour == "red"] == 1.5))
})
