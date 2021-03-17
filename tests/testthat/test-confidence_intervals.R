test_that("Confidence intervals are working properly", {
  fit = blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = T, workers = 4)
  conf.fit = confint(fit, parm = c('wt', 'hp'))
  expect_equal(nrow(conf.fit), 2)
  conf.fit2 = confint(fit, parm = c('wt'))
  expect_equal(nrow(conf.fit2), 1)

  conf.fit90 = confint(fit, level = 0.90)
  conf.fit95 = confint(fit, level = 0.95)
  expect_true(all((conf.fit90[,2] - conf.fit90[,1]) < (conf.fit95[,2] - conf.fit95[,1])))

  sigma90 = sigma(fit, confidence = T, level = 0.90)
  sigma95 = sigma(fit, confidence = T, level = 0.95)
  expect_true(all((sigma90[3] - sigma90[2]) < (sigma95[3] - sigma95[2])))

  newdata = data.frame(wt = 3.2, hp = 146.6)
  predict90 = predict(object = fit, new_data = newdata, confidence = T, level = 0.90)
  predict95 = predict(object = fit, new_data = newdata, confidence = T, level = 0.95)
  expect_true(all((predict90[3] - predict90[2]) < (predict95[3] - predict95[2])))
})

