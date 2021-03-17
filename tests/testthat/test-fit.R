test_that("Fit is class blblm", {
  fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = T, workers = 4)
  expect_identical(class(fit), "blblm")
})
