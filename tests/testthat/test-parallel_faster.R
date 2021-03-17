test_that("Parallelization is faster", {
  n = 2000;x = rnorm(n, mean = 5, sd = 0.5);y = 3*x + rnorm(n, mean = 0.5, sd = 0.4);z = (1/2)*y + rnorm(n, mean = 1, sd = 0.5);m = (1/4)*y + rnorm(n, mean = 1, sd = 0.5);ex = data.frame(x,y,z,m)
  parallel = system.time(blblm(formula = y ~ x + z + m, data = ex, m = 10, B = 5000, parallel = T, workers = 4))
  no_parallel = system.time(blblm(formula = y ~ x + z + m, data = ex, m = 10, B = 5000))
  expect_true(parallel[3] < no_parallel[3])
})
