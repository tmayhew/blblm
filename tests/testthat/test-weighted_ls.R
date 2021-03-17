test_that("lm1 is returning the same coefficients for weighted least squares", {
  n = 1000;x = rnorm(n, mean = 5, sd = 1);y = 3*x + rnorm(n, mean = 0.5, sd = 0.4);z = (1/2)*y + rnorm(n, mean = 1, sd = 1)
  ex = data.frame(x = x,y,z)
  formula = y ~ x + z
  m = model.frame(formula, ex)
  X = model.matrix(formula, m)
  y = model.response(m)
  freqs = as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit1 = lm.wfit(X, y, freqs)
  fit2 = lm(y ~ x + z, data = ex, weights = freqs)
  expect_identical(fit1$coefficients, fit2$coefficients)
})
