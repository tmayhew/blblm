test_that("All example outputs are numeric", {
  fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = T, workers = 4)
  coef.fit = coef(fit)
  expect_equal(class(coef.fit), "numeric")

  confint.fit = confint(fit, c("wt", "hp"))
  expect_equal(class(confint.fit[1,]), "numeric")
  expect_equal(class(confint.fit[2,]), "numeric")

  sigma.fit = sigma(fit)
  expect_equal(class(sigma.fit), "numeric")

  predict.fit = predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
  expect_equal(class(predict.fit), "numeric")

})
