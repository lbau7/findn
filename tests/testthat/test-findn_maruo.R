test_that("findn_maruo finds approximately the correct sample size", {
  set.seed(20210330)
  fun_ttest <- function(n, k) {
    n_mat <- matrix(rnorm(n * k, mean = 5, sd = 10), ncol = k)
    pvals <- apply(n_mat, 2, function(x) t.test(x)$p.value)
    mean(pvals <= 0.05)
  }
  
  res_maruo1 <- findn_maruo(fun = fun_ttest, targ = 0.8, start = 10, k = 100)
  res_maruo2 <- findn_maruo(fun = fun_ttest, targ = 0.9, start = 10, k = 100)
  
  n_true1 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.8)$n)
  n_true2 <- ceiling(power.t.test(delta = 5, sd = 10, type = "one.sample", 
    power = 0.9)$n)
  
  expect_equal(res_maruo1$Point_Estimate_n, n_true1, tol = 2)
  expect_equal(res_maruo2$Point_Estimate_n, n_true2, tol = 2)
})

test_that("findn_maruo stop when power is neither 0.8 nor 0.9", {
  expect_error(findn_maruo(fun = function(n, k) n, targ = 0.85, start = 10,
    k = 100))
})