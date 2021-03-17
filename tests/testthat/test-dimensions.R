test_that("The dimensions of outputs are the same for parallelization", {
  fit1 = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = F)
  fit2 = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = T, workers = 4)
  expect_true(length(coef(fit1))==length(coef(fit2)))
  expect_true(all(dim(confint(fit1, c("wt", "hp")))==dim(confint(fit2, c("wt", "hp")))))
  expect_true(length(sigma(fit1, confidence = TRUE))==length(sigma(fit2, confidence = TRUE)))
  expect_true(all(dim(predict(fit1, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE))==dim(predict(fit2, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE))))
})
