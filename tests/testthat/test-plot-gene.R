test_that("boxplot stats are the quartiles and the full range", {
  x <- c(1:11, 500)
  stats <- brmDE:::calc_boxplot_stat(x)
  expect_equal(unname(stats["middle"]), stats::median(x))
  expect_equal(unname(stats["lower"]), unname(stats::quantile(x, 0.25)))
  expect_equal(unname(stats["upper"]), unname(stats::quantile(x, 0.75)))
  expect_equal(unname(stats[c("ymin", "ymax")]), range(x))
})

test_that("signed sqrt transform round-trips through zero and negatives", {
  x <- c(-4, -1, 0, 1, 9)
  trans <- brmDE:::S_sqrt_trans()
  expect_equal(trans$inverse(trans$transform(x)), x)
})
